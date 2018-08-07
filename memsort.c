#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include "memsort.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"

#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <poll.h>
#include <sys/time.h>

#define	MEM_16G		(1ULL << 34)
KSEQ_DECLARE(gzFile)

//FILE * fbin = NULL; 
struct timeval sort_st, sort_et;

struct timeval bin_write_st, bin_write_et;

struct timeval bin_read_st, bin_read_et;
struct timeval bin_sort_st, bin_sort_et;
struct timeval sbin_write_st, sbin_write_et;

int64_t bin_write_time = 0;
int64_t bin_read_time = 0;
int64_t sbin_write_time = 0;
int64_t bin_sort_time = 0;
int64_t l_pac_global = 0;

mt_entry * mt;
int64_t mt_length = 0;

sort_ot_entry ** sort_ot;
int sort_ot_size = 1000;
int sort_ot_length = 0;

md_ot_entry ** md_ot;
int md_ot_size = 1000;
int md_ot_length = 0;



char * temp_unsorted_file_name;
char * sorted_file_name;
FILE * temp_unsorted_file;
FILE * sorted_file;

int sort_verbose = 0;

void sorting_init_1(int64_t l_pac){            // Assuming this l_pac is size of forward ref strand

    mt = (mt_entry * ) malloc((l_pac) * sizeof(mt_entry));
    sort_ot = (sort_ot_entry ** ) malloc(sort_ot_size * sizeof(sort_ot_entry *));
    md_ot = (md_ot_entry ** ) malloc(md_ot_size * sizeof(md_ot_entry *));
    memset(mt,0,((l_pac)) * sizeof(mt_entry));
    mt_length = l_pac;
    
    memset(sort_ot,0,(sort_ot_size * sizeof(sort_ot_entry *)));
    sort_ot_length = 0;

    memset(md_ot,0,(md_ot_size * sizeof(md_ot_entry *)));
    md_ot_length = 0;
}

void sorting_close(){
    fprintf(stderr,"Closing sorter\n");
    free(mt);
    int i = 0;
    for(i = md_ot_length-1;i>=0;i--){
        free(md_ot[i]);
    }
    for(i = sort_ot_length-1;i>=0;i--){
        free(sort_ot[i]);
    }
    free(md_ot);
    free(sort_ot);
}

void add_sort_ot_entry(bseq1s_t * seq,int64_t file_ptr,int sam_length){
    int64_t fp = file_ptr;
    if(sort_ot_length != 0 && sort_ot_length % sort_ot_size == 0){
        sort_ot = realloc(sort_ot,(sort_ot_length + sort_ot_size) * sizeof(sort_ot_entry *));
    }
    
    sort_ot[sort_ot_length] = (sort_ot_entry * ) malloc(sizeof(sort_ot_entry));
    sort_ot[sort_ot_length]->ref_pos = seq->abs_pos;
    memcpy(&sort_ot[sort_ot_length]->si.file_ptr,&fp,5);
    sort_ot[sort_ot_length]->si.sam_size = sam_length;
    sort_ot[sort_ot_length]->si.correction = seq->correction;
    sort_ot[sort_ot_length]->is_rev = seq->is_rev;
    sort_ot[sort_ot_length]->si.avg_qual = seq->avg_qual;
    sort_ot_length++;
}

void add_md_ot_entry(bseq1s_t * seq, half_mt_entry * md_en){

    // If this returns 0, then the current seq under observation is a duplicate and should be placed.

    /*if(md_ot_length != 0 && md_ot_length % md_ot_size == 0){
        md_ot = realloc(md_ot,(md_ot_length + md_ot_size) * sizeof(md_ot_entry *));
    }*/
    
    if(md_en->mdi.avg_qual < seq->avg_qual){
        // Previous entry was duplicate
        md_en->mdi.avg_qual = seq->avg_qual;
        md_en->mdi.correction = seq->correction;
    }

    /*md_ot[md_ot_length] = (md_ot_entry * ) malloc(sizeof(md_ot_entry));
    md_ot[md_ot_length]->ref_pos = seq->abs_pos - seq->correction;
    md_ot[md_ot_length]->mdi.avg_qual = seq->avg_qual;
    md_ot[md_ot_length]->mdi.correction = seq->correction;
    md_ot_length++;*/
}


void add_mt_entry_1(bseq1s_t *seq, int64_t file_ptr, int sam_length){


    half_mt_entry * sort_en;
    half_mt_entry * md_en;
    int64_t fp;
    assert(seq->abs_pos - seq->correction >= 0);
    if(seq->is_rev){
        sort_en = &mt[seq->abs_pos].re;
        md_en = &mt[seq->abs_pos - seq->correction].re;
    }
    else {
        sort_en = &mt[seq->abs_pos].fe;
        md_en = &mt[seq->abs_pos - seq->correction].fe;
    }

    if(sort_en->si.sam_size == 0){
        //Sort entry was free
        fp = file_ptr;
        memcpy(&(sort_en->si.file_ptr),&fp,5);
        sort_en->si.sam_size = (int16_t)sam_length;
        sort_en->si.correction = seq->correction;
        sort_en->si.avg_qual = seq->avg_qual;
    }
    else{
        // Add entry to Sort Overflow Table
        add_sort_ot_entry(seq,file_ptr,sam_length);
    }

    if(md_en->mdi.avg_qual == 0){
        // MD entry was free
        md_en->mdi.avg_qual = seq->avg_qual;
        md_en->mdi.correction = seq->correction;
    }
    else{
        // Add entry to MarkDuplicate Overflow Table
        // TODO: Add checking mate information for pair ended reads

        add_md_ot_entry(seq, md_en);
    }
    return;
}


void process_sam_record_1(bseq1s_t *seq, int64_t file_ptr, int sam_length, double *sam_rtime){

    double sam_rtime_start = realtime();
    if(seq == NULL){
        return;
    }
    if(seq->abs_pos <= 0){
        // TODO: Manage unmapped entries
        // For now remove them
        return;
    }
    else{
        if(seq->is_rev){
            add_mt_entry_1(seq, file_ptr, sam_length);
        }
        else{
            add_mt_entry_1(seq, file_ptr, sam_length);
        }
    }

    *sam_rtime += realtime() - sam_rtime_start;
}

int sort_comparator(const void **p, const void **q) 
{
    sort_ot_entry* l = *((sort_ot_entry **)p);
    sort_ot_entry* r = *((sort_ot_entry **)q); 
    
    if(l->ref_pos < r->ref_pos){
        return -1;
    }
    else{
        return 1;
    }
}


int is_duplicate(sort_info * sort_in, md_info * md_in){

    if(sort_in->avg_qual == 0 || md_in->avg_qual == 0){
        return 0;
    }

    if(sort_in->correction == md_in->correction && sort_in->avg_qual == md_in->avg_qual){
        return 0;
    } 
    else {
        return 1;
    }
}


void get_sorted_sam(){

    fclose(temp_unsorted_file);
    int64_t i = 0;
   
    struct timeval qsort_st, qsort_et;
    gettimeofday(&qsort_st,NULL);
    fprintf(stderr,"Total entries in sort_ot : %d\n",sort_ot_length);
    if(sort_ot_length > 0){
        fprintf(stderr,"Sorting the sort_ot\n");
        qsort(sort_ot,sort_ot_length,sizeof(sort_ot_entry *),sort_comparator);
        fprintf(stderr,"Finished sort\n");
    }
    gettimeofday(&qsort_et,NULL);
    double qsort_time = ((double)qsort_et.tv_sec - (double)qsort_st.tv_sec) + (double)((double)(qsort_et.tv_usec - qsort_st.tv_usec) / (double)(1000000));
    fprintf(stderr,"Qsort time : %f\n",qsort_time);
    int64_t count_duplicates = 0;
    int64_t head = 0;
    int64_t md_index = 0;
    if(sort_verbose >= 3){
        fprintf(stderr,"mt_length after : %ld\n",mt_length);
    }

    gettimeofday(&qsort_st,NULL);
    for(i = 0;i<mt_length;i++){
        if(head < sort_ot_length){
            while (sort_ot[head]->ref_pos < i){
                md_index = sort_ot[head]->ref_pos - sort_ot[head]->si.correction;
                if(sort_ot[head]->is_rev){
                        count_duplicates += is_duplicate(&sort_ot[head]->si, &mt[md_index].re.mdi);
                }
                else{
                        count_duplicates += is_duplicate(&sort_ot[head]->si, &mt[md_index].fe.mdi);
                }

                head++;
                if(head >= sort_ot_length){
                    break;
                }
            }
        }
        md_index = i - mt[i].fe.si.correction;
        count_duplicates += is_duplicate(&mt[i].fe.si, &mt[md_index].fe.mdi);
        
        md_index = i - mt[i].re.si.correction;
        count_duplicates += is_duplicate(&mt[i].re.si, &mt[md_index].re.mdi);
    }
    gettimeofday(&qsort_et,NULL);
    qsort_time = ((double)qsort_et.tv_sec - (double)qsort_st.tv_sec) + (double)((double)(qsort_et.tv_usec - qsort_st.tv_usec) / (double)(1000000));
    fprintf(stderr,"Mark duplicate time : %f\n",qsort_time);
    fprintf(stderr,"Total duplicates : %d\n",count_duplicates);

    return;

}

void usage(){
    printf("./bwa sort -I <input sam> -O <output sam> -v <verbose level> \n");
}


int get_sequence(char * line, size_t size, int64_t * chr_start, int64_t * chr_start_array, int * chr_index){
   
    if(sort_verbose >= 5){
        printf("Line : %s",line);
    } 
    if(line[1] != 'S'){
        return -1;
    }

    char * splits = NULL;
    char delim = '\t';

    chr_start_array[*chr_index] = *chr_start;
    *chr_index += 1;
    
    splits = strtok(line,&delim); 
    while(splits != NULL){
        if(splits[0] == 'L'){
            *chr_start += (int64_t)strtoll(splits + 3, NULL,10);
            break;
        }
        else{
            splits = strtok(NULL,&delim);
        }        
    }
    return 1;
}


int get_chr_num(char * chr){

    int size = strlen(chr);
    if(size < 4){
        return -1;
    }
    int i = 0;
    int chr_num = 0;
    for(i=3;i<size;i++){
        if(chr[i] == 'X' || chr[i] == 'x'){
            chr_num = 23;
            break;
        }
        else if(chr[i] == 'Y' || chr[i] == 'y'){
            chr_num = 24;
            break;
        }
        else if((int)chr[i] >= (int)'1' && (int)chr[i] <= (int)'9'){
            chr_num = (int)strtol(&chr[i],NULL,10);
            break;
        }
    }
    return chr_num;
}

int get_correction_from_CIGAR(char * cigar){
    int len = strlen(cigar);
    if(len < 4){
        return 0;
    }
    int i = 0;
    int correction = 0;
    int field_num = 0;

    for(i = 0;i<len;i++){
        if( (int)cigar[i] >= (int)'0' && (int)cigar[i] <= (int)'9'){
            field_num = field_num * 10 + (int)cigar[i] - 48;
        }
        else if( cigar[i] == 'M'){
            break;
        }
        else{
            correction += field_num;
            field_num = 0;
        }
    }
    return correction;
}

int get_avg_qual( char * qual, int len){
    //int len = strlen(qual);
    int i = 0;
    int total_qual = 0;
    for(i=0;i<len;i++){
        total_qual += (int)qual[i] - 33;
    }
    return (total_qual / len);
}


void get_record(char * line, size_t line_size, int *flags, int *chr_num, int64_t *pos, int *correction, int *avg_qual){
    int i = 0;
    int field_num = 1;

    int corr = 0;
    int stop = 0;

    int l_seq = 0;

    for(i = 0;i<line_size;i++){
        if(line[i] == '\t'){
            field_num++;
            stop = 0;
            if(field_num >= 12){
                break;
            }
        }
        else if(field_num == 2){
            *flags = ((*flags) * 10) + ((int)line[i] - 48);
        }
        else if(field_num == 3 && stop == 0){
            if(line[i] == 'C' || line[i] == 'c' || line[i] == 'H' || line[i] == 'h' || line[i] == 'R' || line[i] == 'r'){
                continue;    
            }
            else if((int)line[i] >= (int)'0' && (int)line[i] <= (int)'9'){
                *chr_num = (*chr_num * 10) + (int)line[i] - 48;
            }
            else if(line[i] == 'X' || line[i] == 'x'){
                *chr_num = 23;
            }
            else if(line[i] == 'Y' || line[i] == 'y'){
                *chr_num = 24;
            }
            else{
                *chr_num = -1;
                stop = 1;
            }
        }
        else if(field_num == 4 && stop == 0){
            if((int)line[i] >= (int)'0' && (int)line[i] <= (int)'9'){
                *pos = (*pos * 10) + (int)line[i] - 48;
            }
            else{
                *pos = -1;
                stop = 1;
            }
        }
        else if(field_num == 6 && stop == 0){
            if((int)line[i] >= (int)'0' && (int)line[i] <= (int)'9'){
                corr = corr * 10 + (int)line[i] - 48;
            }
            else if(line[i] != 'M'){
                *correction += corr;
                corr = 0;
            }
            else if(line[i] == 'M'){
                stop = 1;
            }
        }
        else if(field_num == 11){
            *avg_qual = *avg_qual + (int)line[i] - 33;
            l_seq++;
        }
    }

    *avg_qual = *avg_qual / l_seq;
    return;
}



void get_sam_record(char * line, size_t line_size, int64_t * chr_start_array, int64_t file_ptr, double * sam_rtime){
    
    char * splits = 0;
    int field_num = 0;

    int flags = 0;          // Field 2
    int chr_num = 0;            // Field 3
    int64_t pos = 0;            // Field 4
    int correction = 0;         // From CIGAR Field 6
    int l_seq = 0;
    int avg_qual = 0;           // Field 11
    char delim = '\t';


    get_record(line, line_size, &flags, &chr_num, &pos, &correction, &avg_qual);

    // TODO create the table


    if(sort_verbose >= 5){
        fprintf(stderr,"flags : %d, chr_num : %d, pos : %ld, corr : %d\n",flags,chr_num,pos,correction);
    }
    if(pos != -1 && chr_num != -1){
        pos += chr_start_array[chr_num-1];
    }

    bseq1s_t seq;
    seq.l_seq = 0;
    seq.id = 0;
    seq.read_id = 0;
    seq.abs_pos = pos;
    seq.correction = correction;
    seq.is_rev = (flags & 0x10) != 0 ? 1 : 0;
    seq.avg_qual = avg_qual;
    
    seq.name = NULL;
    seq.comment = NULL;
    seq.seq = NULL;
    seq.qual = NULL;
    seq.sam = NULL;

    if(chr_num != -1){
        process_sam_record_1(&seq,file_ptr,line_size, sam_rtime);
    }

} 


void print_seqs(int64_t * chr_start){
    int i = 0;
    for(i = 0;i<24;i++){
        fprintf(stderr,"[%d] %ld\n",i+1,chr_start[i]);
    }
}

/*void print_strings(char * in, int len){
    int i;
    for(i = 0;i<len;i++){
        printf("%s",in[i]);
    }
    printf("\n");
}*/

void sort_start(){
    char * line;
    size_t line_size;
    size_t ret;
    
    int64_t seq_start[24]; 
    // X is labelled as 23 
    // Y is labelled as 24


    line = NULL;
    line_size = 0;
    int i = 0;

    int64_t chr_start = 0;
    int chr_index = 0;
    int64_t * chr_start_array = (int64_t *) malloc(24 * sizeof(int64_t));
    int64_t file_ptr = 0;


    int rc = 0;

    double rtime;
    double prev_rtime, curr_rtime;
    int64_t read_count = 0;

    rtime = realtime();

    while(!feof(temp_unsorted_file)){
        file_ptr = (int64_t) ftell(temp_unsorted_file);
        ret = getline(&line,&line_size,temp_unsorted_file);
        if(ret != -1){
            if(line[0] == '@'){
                rc = get_sequence(line, line_size, &chr_start, chr_start_array, &chr_index);
                /*if(rc == -1){
                    break;
                }*/
            }
            else{
                break;
            }
        }
    }

    if(sort_verbose >= 3){
        fprintf(stderr,"Lpac : %ld\n", chr_start);
    }

    l_pac_global = chr_start;
    if(sort_verbose >= 3){
        fprintf(stderr,"Initializing sorter tables\n");
    }
    sorting_init_1(chr_start);

    if(sort_verbose >= 3){
        fprintf(stderr,"Mt_length after init : %ld\n",mt_length);
    }

    if(sort_verbose >= 3){
        fprintf(stderr,"Finished initialization | %f secs\n",realtime() - rtime);
        fprintf(stderr,"Starting reading input sam file\n");
    }
    double sam_process_time = 0.0;
    if(line != NULL){
        file_ptr = (int64_t) ftell(temp_unsorted_file);
        get_sam_record(line, strlen(line), chr_start_array, file_ptr, &sam_process_time);
        read_count++;
    }

    prev_rtime = realtime() - rtime;
    while(!feof(temp_unsorted_file)){
        int64_t file_ptr = (int64_t) ftell(temp_unsorted_file);
        ret = getline(&line,&line_size,temp_unsorted_file);
        if(ret != -1){
            get_sam_record(line, strlen(line), chr_start_array, file_ptr, &sam_process_time);
            read_count++;
        }

        if(read_count % 1000000 == 0 && sort_verbose >= 3){
            
            curr_rtime = realtime() - rtime;
            fprintf(stderr,"Read %ld reads in  %f secs\n",read_count,sam_process_time);
            prev_rtime = curr_rtime;
            //sam_process_time = 0.0;
        }
    }
    if(sort_verbose >= 3){
        fprintf(stderr,"Read %ld reads | %f secs\n",read_count,realtime() - rtime);
        fprintf(stderr,"Total sam process time : %f\n",sam_process_time);
    }

    if(line)
        free(line);
    free(chr_start_array);
    get_sorted_sam();

    sorting_close();
    
}

int main_memsort(int argc, char *argv[]){
    int rargc = 2;
    int i = 0;
    int input_file = 0;
    int output_file = 0; 
    if(argc < rargc){
        usage();
        return 0;
    }
   
    for(i=0;i<argc;i+=2){
        if(strcmp(argv[i],"-I") == 0) {
            temp_unsorted_file_name = argv[i+1];
            input_file = 1;
            continue;
        }
        if(strcmp(argv[i],"-O") == 0) {
            if(argc > i+1){
                sorted_file_name = argv[i+1];
                output_file = 1;
            }
            else{
                usage();
                return -1;
            }
            continue;
        }
        if(strcmp(argv[i],"-v") == 0){
            if(argc > i+1){
                sort_verbose = (int)strtol(argv[i+1],NULL,10);
            }
            else{
                usage();
                return -1;
            }
            continue;
        }
    }

    temp_unsorted_file = fopen(temp_unsorted_file_name,"r");
    if(temp_unsorted_file == NULL){
        fprintf(stderr,"File %s not found\n",temp_unsorted_file_name);
        return -1;
    }

    if(output_file == 1){
        sorted_file = fopen(sorted_file_name,"w");
        stdout = sorted_file;
    }

    sort_start();


    

    return 0;

}
