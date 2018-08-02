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
dpt_entry ** dpt;
umt_entry ** umt;

int dpt_size = 1000;
int dpt_length = 0;
int umt_size = 1000;
int umt_length = 0;
int64_t mt_length = 0;
char * temp_unsorted_file_name;
char * sorted_file_name;
FILE * temp_unsorted_file;
FILE * sorted_file;

int sort_verbose = 0;



void print_seqs(int64_t * seqs){
    int i = 0;
    for(i=0;i<24;i++){
        printf("Seq[%d] : %ld\n",i,seqs[i]);
    }
}



void sorting_init(int64_t l_pac){            // Assuming this l_pac is size of forward ref strand

    fprintf(stderr,"Initializing sorter\n");
    mt = (mt_entry * ) malloc((l_pac) * sizeof(mt_entry));
    dpt = (dpt_entry ** ) malloc(dpt_size * sizeof(dpt_entry *));
    umt = (umt_entry ** ) malloc(umt_size * sizeof(umt_entry *));
    memset(mt,0,((l_pac) + 1) * sizeof(mt_entry));
    memset(dpt,0,(dpt_size * sizeof(dpt_entry *)));
    dpt_length = 0;
    mt_length = l_pac;
    //TODO: Define temp_unsorted_file;
    temp_unsorted_file_name = "/home/centos/extra_ssd/tmp_unsorted.sam";
    temp_unsorted_file = fopen(temp_unsorted_file_name,"w");
}

void sorting_init_1(int64_t l_pac){            // Assuming this l_pac is size of forward ref strand

    fprintf(stderr,"Initializing sorter\n");
    mt = (mt_entry * ) malloc((l_pac) * sizeof(mt_entry));
    dpt = (dpt_entry ** ) malloc(dpt_size * sizeof(dpt_entry *));
    umt = (umt_entry ** ) malloc(umt_size * sizeof(umt_entry *));
    memset(mt,0,((l_pac) + 1) * sizeof(mt_entry));
    memset(dpt,0,(dpt_size * sizeof(dpt_entry *)));
    dpt_length = 0;
    umt_length = 0;
    mt_length = l_pac;
}

void sorting_close(){
    fprintf(stderr,"Closing sorter\n");
    free(mt);
    int i = 0;
    for(i = dpt_length-1;i>=0;i--){
        free(dpt[i]);
    }
    for(i = umt_length-1;i>=0;i--){
        free(umt[i]);
    }
    free(dpt);
    free(umt);
}

void add_dpt_entry_mt(int64_t ref_pos,half_mt_entry * in_mt_entry){
    if(dpt_length != 0 && dpt_length % dpt_size == 0){
        dpt = realloc(dpt,(dpt_length + dpt_size) * sizeof(dpt_entry *));
    }
    dpt[dpt_length] = (dpt_entry *) malloc(sizeof(dpt_entry));
    dpt[dpt_length]->ref_pos = ref_pos;
    memcpy(&dpt[dpt_length]->file_ptr,&in_mt_entry->file_ptr,5);
    dpt[dpt_length]->sam_size = in_mt_entry->sam_size;
    dpt_length++;
}

void add_dpt_entry_seq(int64_t file_ptr, bseq1_t *seq){
    int fp = file_ptr;
    if(dpt_length != 0 && dpt_length % dpt_size == 0){
        dpt = realloc(dpt,(dpt_length + dpt_size) * sizeof(dpt_entry *));
    }
    dpt[dpt_length] = (dpt_entry *) malloc(sizeof(dpt_entry));
    dpt[dpt_length]->ref_pos = seq->abs_pos;
    memcpy(&dpt[dpt_length]->file_ptr,&fp,5);
    dpt[dpt_length]->sam_size = strlen(seq->sam);
    dpt_length++;
}


void add_dpt_entry_seq_1(int64_t file_ptr, bseq1_t *seq, int sam_length){
    int fp = file_ptr;
    if(dpt_length != 0 && dpt_length % dpt_size == 0){
        dpt = realloc(dpt,(dpt_length + dpt_size) * sizeof(dpt_entry *));
    }
    dpt[dpt_length] = (dpt_entry *) malloc(sizeof(dpt_entry));
    dpt[dpt_length]->ref_pos = seq->abs_pos;
    memcpy(&dpt[dpt_length]->file_ptr,&fp,5);
    dpt[dpt_length]->sam_size = sam_length;
    dpt_length++;
}


void add_mt_entry(bseq1_t *seq, half_mt_entry * in_en){
    if(in_en->avg_qual == 0){
        int64_t fp = (int64_t)ftell(temp_unsorted_file);
        memcpy(&(in_en->file_ptr),&fp,5);
        in_en->avg_qual = seq->avg_qual;
        in_en->sam_size = strlen(seq->sam);
    } 
    else{
        int64_t fp = (int64_t)ftell(temp_unsorted_file);
        if(in_en->avg_qual < seq->avg_qual) {
            add_dpt_entry_mt(seq->abs_pos,in_en);
            memcpy(&(in_en->file_ptr),&fp,5);
            in_en->avg_qual = seq->avg_qual;
            in_en->sam_size = strlen(seq->sam);
        }
        else {
            add_dpt_entry_seq(fp,seq);
        }
    }
    return;
}

void add_mt_entry_1(bseq1_t *seq, half_mt_entry * in_en, int64_t file_ptr, int sam_length){
    if(in_en->avg_qual == 0){
        int64_t fp = file_ptr;
        memcpy(&(in_en->file_ptr),&fp,5);
        in_en->avg_qual = seq->avg_qual;
        in_en->sam_size = (int16_t)sam_length;
    } 
    else{
        int64_t fp = file_ptr;
        if(in_en->avg_qual < seq->avg_qual) {
            add_dpt_entry_mt(seq->abs_pos,in_en);
            memcpy(&(in_en->file_ptr),&fp,5);
            in_en->avg_qual = seq->avg_qual;
            in_en->sam_size = (int16_t)sam_length;
        }
        else {
            add_dpt_entry_seq_1(fp,seq, sam_length);
        }
    }
    return;
}

void add_umt_entry(bseq1_t *seq){
    int64_t fp = (int64_t)ftell(temp_unsorted_file);
    if(umt_length != 0 && umt_length % dpt_size == 0){
        umt = realloc(umt,(umt_length + umt_size) * sizeof(umt_entry *));
    }
    umt[umt_length] = (umt_entry *) malloc(sizeof(umt_entry));
    memcpy(&umt[umt_length]->file_ptr,&fp,5);
    umt[umt_length]->sam_size = strlen(seq->sam);
    umt_length++;
}

void add_umt_entry_1(bseq1_t *seq, int64_t file_ptr, int sam_length){
    int64_t fp = file_ptr;
    if(umt_length != 0 && umt_length % dpt_size == 0){
        umt = realloc(umt,(umt_length + umt_size) * sizeof(umt_entry *));
    }
    umt[umt_length] = (umt_entry *) malloc(sizeof(umt_entry));
    memcpy(&umt[umt_length]->file_ptr,&fp,5);
    umt[umt_length]->sam_size = sam_length;
    umt_length++;
}

void process_sam_record(bseq1_t *seq){
    if(seq == NULL){
        return;
    }
    int i = 0;
    if(seq->abs_pos >= mt_length){
        add_umt_entry(seq); 
    }
    else{
        if(seq->is_rev){
            add_mt_entry(seq,&(mt[seq->abs_pos].re));
        }
        else{
            add_mt_entry(seq,&(mt[seq->abs_pos].fe));
        }
    }
    err_fputs(seq->sam, temp_unsorted_file);
}

void process_sam_record_1(bseq1_t *seq, int64_t file_ptr, int sam_length){
    if(seq == NULL){
        return;
    }
    int i = 0;
    if(seq->abs_pos >= mt_length){
        add_umt_entry_1(seq,file_ptr,sam_length); 
    }
    else{
        if(seq->is_rev){
            add_mt_entry_1(seq,&(mt[seq->abs_pos].re), file_ptr, sam_length);
        }
        else{
            add_mt_entry_1(seq,&(mt[seq->abs_pos].fe), file_ptr, sam_length);
        }
    }
}

int comparator(const void **p, const void **q) 
{
    dpt_entry* l = *((dpt_entry **)p);
    dpt_entry* r = *((dpt_entry **)q); 
    
    if(l->ref_pos < r->ref_pos){
        return -1;
    }
    else{
        return 1;
    }
}

void print_mt_entry(half_mt_entry * in,FILE * fin,FILE * fout){
    if(in->avg_qual != 0){
        //TODO: Read from temp_file
        int64_t fp = 0;
        memcpy(&fp,&in->file_ptr,5);
        fseek(fin,fp,SEEK_SET);
        char * in_sam = malloc(in->sam_size + 1);
        memset(in_sam,0,in->sam_size+1);
        fread(in_sam,sizeof(char),in->sam_size,fin);
        in_sam[in->sam_size] = '\0';
        if(sort_verbose >= 5){
            fprintf(fout,"File ptr : %ld, sam_size : %d\n",fp,in->sam_size);
        }
        err_fputs(in_sam, fout);
        free(in_sam);
    }
} 

void print_dpt_entry(dpt_entry * in,FILE * fin,FILE * fout){
    int64_t fp = 0;
    memcpy(&fp,&in->file_ptr,5);
    fseek(fin,fp,SEEK_SET);
    char * in_sam = malloc(in->sam_size + 1);
    memset(in_sam,0,in->sam_size+1);
    fread(in_sam,sizeof(char),in->sam_size,fin);
    in_sam[in->sam_size] = '\0';

    //TODO: Add extra flag to sam for duplicates
    //fprintf(fout,"*****");
    err_fputs(in_sam, fout);
    free(in_sam);
}

void print_umt_entry(umt_entry * in,FILE * fin, FILE * fout){
    int64_t fp = 0;
    memcpy(&fp,&in->file_ptr,5);
    fseek(fin,fp,SEEK_SET);
    char * in_sam = malloc(in->sam_size + 1);
    memset(in_sam,0,in->sam_size+1);
    fread(in_sam,sizeof(char),in->sam_size,fin);
    in_sam[in->sam_size] = '\0';
    err_fputs(in_sam, fout);
    free(in_sam);
}

void print_full_dpt(){
    int i = 0;
    fprintf(stderr,"Sorted dpt : \n");
    for(i = 0;i<dpt_length;i++){
        fprintf(stderr,"Ref pos : %ld\n",dpt[i]->ref_pos);
    }
}

void get_sorted_sam(){

    fclose(temp_unsorted_file);
    temp_unsorted_file = fopen(temp_unsorted_file_name,"r");
    int64_t i = 0;
   
    struct timeval qsort_st, qsort_et;
    gettimeofday(&qsort_st,NULL);
    fprintf(stderr,"Total duplicates detected : %d\n",dpt_length);
    if(dpt_length > 0){
        fprintf(stderr,"Sorting the duplicate table\n");
        qsort(dpt,dpt_length,sizeof(dpt_entry *),comparator);
        fprintf(stderr,"Finished duplicate table sort\n");
    }
    gettimeofday(&qsort_et,NULL);
    double qsort_time = ((double)qsort_et.tv_sec - (double)qsort_st.tv_sec) + (double)((double)(qsort_et.tv_usec - qsort_st.tv_usec) / (double)(1000000));
    fprintf(stderr,"Qsort time : %f\n",qsort_time);

    fprintf(stderr,"Starting sorted prints\n");
    int head = 0;
    for(i=0;i<mt_length;i++){
        while(head < dpt_length && dpt[head]->ref_pos < i){
            print_dpt_entry(dpt[head],temp_unsorted_file,stdout); 
            head++;
        }
        if(sort_verbose >= 5 && (mt[i].fe.avg_qual != 0 || mt[i].re.avg_qual != 0)){
            fprintf(stdout,"Index : %ld, ",i);
        }
        print_mt_entry(&mt[i].fe,temp_unsorted_file,stdout);
        print_mt_entry(&mt[i].re,temp_unsorted_file,stdout);
    } 

    for(i=0;i<umt_length;i++){
        print_umt_entry(umt[i],temp_unsorted_file,stdout);
    }
    fclose(temp_unsorted_file);
}

void usage(){
    printf("./bwa sort -I <input sam> -O <output sam> -v <verbose level> \n");
}


int get_sequence(char * line, size_t size, int64_t * chr_start, int64_t * chr_start_array, int * chr_index){
    
    if(line[1] != 'S'){
        return -1;
    }

    int i = 0;
    char * splits;
    char delim = '\t';
    int64_t len = 0;

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

int get_avg_qual( char * qual){
    int len = strlen(qual);
    int i = 0;
    int total_qual = 0;
    for(i=0;i<len;i++){
        total_qual += (int)qual[i] - 33;
    }
    return (total_qual / len);
}


void get_sam_record(char * line, size_t line_size, int64_t * chr_start_array, int64_t file_ptr, double * sam_rtime){
    
    double sam_rtime_start = realtime();
    char * splits = 0;
    int field_num = 0;

    int16_t flags = 0;          // Field 2
    int chr_num = 0;            // Field 3
    int64_t pos = 0;            // Field 4
    int correction = 0;         // From CIGAR Field 6
    int avg_qual = 0;           // Field 11
    char delim = '\t';

    splits = strtok(line,&delim);
    while(splits != NULL){
        field_num++;
        if(field_num == 2){
            flags = (int)strtol(splits,NULL,10);
        }
        else if(field_num == 3){
            chr_num = get_chr_num(splits); 
        }
        else if(field_num == 4){
            pos = (int64_t)strtoll(splits, NULL, 10);
            if(chr_num == -1){
                pos = l_pac_global;
            }
            pos += chr_start_array[chr_num-1];
        }
        else if(field_num == 6){
            correction = get_correction_from_CIGAR(splits);
        }
        else if(field_num == 11){
            avg_qual = get_avg_qual(splits);
            break;
        }
        splits = strtok(NULL,&delim);
    }

    // TODO create the table

    bseq1_t seq;
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

    if(sort_verbose >= 5){
        printf("pos : %ld,correction : %d\n",pos,correction);
    }
    process_sam_record_1(&seq,file_ptr,line_size);

    *sam_rtime += realtime() - sam_rtime_start;
} 


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
    int64_t read_count;

    rtime = realtime();

    while(!feof(temp_unsorted_file)){
        int64_t file_ptr = (int64_t) ftell(temp_unsorted_file);
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

    if(sort_verbose >= 5){
        printf("Lpac : %ld\n", chr_start);
    }

    l_pac_global = chr_start;
    sorting_init_1(chr_start);

    if(sort_verbose >= 3){
        fprintf(stderr,"Finished initialization | %f secs\n",realtime() - rtime);
        fprintf(stderr,"Starting reading input sam file\n");
    }
    double sam_process_time = 0.0;
    if(line != NULL){
        int64_t file_ptr = (int64_t) ftell(temp_unsorted_file);
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
            fprintf(stderr,"Read %ld reads in %f secs, read processing time : %f secs \t|\t Total : %f secs \n",read_count,curr_rtime - prev_rtime,sam_process_time,curr_rtime);
            fprintf(stderr,"[DATA],%ld,%f,%f,%f\n",read_count,curr_rtime - prev_rtime,sam_process_time,curr_rtime);
            prev_rtime = curr_rtime;
            sam_process_time = 0.0;
        }
    }
    if(sort_verbose >= 3){
        fprintf(stderr,"Read %ld reads | %f secs\n",read_count,realtime() - rtime);
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
