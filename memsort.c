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
#define QUEUE_SIZE 100
#define INIT_OVERFLOW_TABLE_SIZE 1000000
KSEQ_DECLARE(gzFile)

//FILE * fbin = NULL; 

char * sorted_file_name;
FILE * sorted_file;

int sort_verbose = 0;

struct timespec ini_start, ini_end;
double ini_time = 0.0;

struct timespec read_start, read_end;
double read_time = 0.0;

struct timespec finish_start, finish_end;
double finish_time = 0.0;


void struct_init(sort_struct_t * sld, int64_t length){            

    sld->mt_length = length;

    sld->sort_ot_size = INIT_OVERFLOW_TABLE_SIZE;
    sld->sort_ot_length = 0;

    sld->md_ot_size = INIT_OVERFLOW_TABLE_SIZE;
    sld->md_ot_length = 0;


    sld->mt = (mt_entry * ) malloc(length * sizeof(mt_entry));
    memset(sld->mt,0,(length * sizeof(mt_entry)));
    
    sld->sort_ot = (sort_ot_entry ** ) malloc(INIT_OVERFLOW_TABLE_SIZE * sizeof(sort_ot_entry *));
    memset(sld->sort_ot,0,(INIT_OVERFLOW_TABLE_SIZE * sizeof(sort_ot_entry *)));

    sld->md_ot = (md_ot_entry ** ) malloc(INIT_OVERFLOW_TABLE_SIZE * sizeof(md_ot_entry *));
    memset(sld->md_ot,0,(INIT_OVERFLOW_TABLE_SIZE * sizeof(md_ot_entry *)));
}

void struct_delete(sort_struct_t * sld){

    free(sld->mt);
    int i = 0;
    for(i = sld->md_ot_length-1;i>=0;i--){
        free(sld->md_ot[i]);
    }
    for(i = sld->sort_ot_length-1;i>=0;i--){
        free(sld->sort_ot[i]);
    }
    free(sld->md_ot);
    free(sld->sort_ot);
}

void add_sort_ot_entry(bseq1s_t * seq,sort_struct_t * sld){

    if(sld->sort_ot_length != 0 && sld->sort_ot_length % sld->sort_ot_size == 0){
        sld->sort_ot = realloc(sld->sort_ot,(sld->sort_ot_length + sld->sort_ot_size) * sizeof(sort_ot_entry *));
    }
   
    int index = sld->sort_ot_length;

    sld->sort_ot[index] = (sort_ot_entry * ) malloc(sizeof(sort_ot_entry));
    sld->sort_ot[index]->ref_pos = seq->abs_pos;
    sld->sort_ot[index]->is_rev = seq->is_rev;

    sld->sort_ot[index]->si.read_id = seq->read_id;
    sld->sort_ot[index]->si.correction = seq->correction;
    sld->sort_ot[index]->si.avg_qual = seq->avg_qual;
    sld->sort_ot[index]->si.flags = seq->flags;
    sld->sort_ot[index]->si.mate_diff = seq->mate_diff;
    sld->sort_ot_length++;
}


void add_md_ot_entry(bseq1s_t * seq, half_mt_entry * md_en, sort_struct_t * sld){


    // Md_en is not empty when it enters this method
    // Check if the mate_diff is the same in the seq and md_en


    if(md_en->mdi.mate_diff == seq->mate_diff){
        // The inputs are truly duplicates
        if(md_en->mdi.avg_qual >= seq->avg_qual){
            // Seq is the duplicates
            // Do nothing
        }
        else{
            // md_en is the duplicate
            md_en->mdi.avg_qual = seq->avg_qual;
            md_en->mdi.correction = seq->correction;
            md_en->mdi.mate_diff = seq->mate_diff;
        }
    }
    else {
        // None of the entries are duplicates and add seq to the ot
        if(sld->md_ot_length != 0 && sld->md_ot_length % sld->md_ot_size == 0){
            sld->md_ot = realloc(sld->md_ot,(sld->md_ot_length + sld->md_ot_size) * sizeof(md_ot_entry *));
        }

        md_ot_entry ** tmp_md_ot = sld->md_ot;
        int index = sld->md_ot_length;

        tmp_md_ot[index] = (md_ot_entry * ) malloc(sizeof(md_ot_entry));
        tmp_md_ot[index]->ref_pos = seq->abs_pos - seq->correction;
        tmp_md_ot[index]->is_rev = seq->is_rev;
        tmp_md_ot[index]->mdi.avg_qual = seq->avg_qual;
        tmp_md_ot[index]->mdi.correction = seq->correction;
        tmp_md_ot[index]->mdi.mate_diff = seq->mate_diff;
        sld->md_ot_length++;
    }

    return;
}


void add_mt_entry_1(bseq1s_t *seq, sort_struct_t * sld){


    half_mt_entry * sort_en;
    half_mt_entry * md_en;
    assert(seq->abs_pos - seq->correction >= 0);
    if(seq->is_rev){
        sort_en = &sld->mt[seq->abs_pos].re;
        md_en = &sld->mt[seq->abs_pos - seq->correction].re;
    }
    else {
        sort_en = &sld->mt[seq->abs_pos].fe;
        md_en = &sld->mt[seq->abs_pos - seq->correction].fe;
    }

    if(sort_en->si.read_id == 0){
        //Sort entry was free
        sort_en->si.read_id = seq->read_id;
        sort_en->si.correction = seq->correction;
        sort_en->si.avg_qual = seq->avg_qual;
        sort_en->si.flags = seq->flags;
        sort_en->si.mate_diff = seq->mate_diff;
    }
    else{
        // Add entry to Sort Overflow Table
        add_sort_ot_entry(seq,sld);
    }

    if(md_en->mdi.avg_qual == 0){
        // MD entry was free
        md_en->mdi.avg_qual = seq->avg_qual;
        md_en->mdi.correction = seq->correction;
        md_en->mdi.mate_diff = seq->mate_diff;
    }
    else{
        // Add entry to MarkDuplicate Overflow Table
        // TODO: Add checking mate information for pair ended reads

        add_md_ot_entry(seq, md_en, sld);
    }
    return;
}


void process_sam_record_1(bseq1s_t *seq, sort_struct_t * sld){

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
            add_mt_entry_1(seq, sld);
        }
        else{
            add_mt_entry_1(seq, sld);
        }
    }

}

int sort_ot_comparator(const void **p, const void **q) 
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

int md_ot_comparator(const void **p, const void **q) 
{
    md_ot_entry* l = *((md_ot_entry **)p);
    md_ot_entry* r = *((md_ot_entry **)q); 
    
    if(l->ref_pos < r->ref_pos){
        return -1;
    }
    else{
        return 1;
    }
}


int is_duplicate(sort_info * sort_in, md_info * md_in){

    if(sort_in->read_id == 0 || md_in->avg_qual == 0){
        return 0;
    }
    if((sort_in->flags & 0x02) != 0){
        // Properly paired, check mate as well
        if(sort_in->mate_diff == md_in->mate_diff){
            // Could be the same read or a different read

            if(sort_in->correction == md_in->correction && sort_in->avg_qual == md_in->avg_qual){
                // Same read
                return 0;
            }
            else{
                return 1;
            }

        }
        else{
            return 0;
        }
    }
    else{
        // Not properly paired, only check correction and avg_qual scores
        if(sort_in->correction == md_in->correction && sort_in->avg_qual == md_in->avg_qual){
            // Same read
            return 0;
        }
        else{
            return 1;
        }

    }
}

int check_duplicate(sort_info * sort_in, int64_t md_in_ref_pos, int64_t md_ot_index, int is_rev, sort_struct_t * sld){
    //assert(md_ot_index >= 0 && md_ot_index < md_ot_length);


    // Check main table entry first

    int64_t mt_ref_pos = md_in_ref_pos + sort_in->correction;
    if(sort_verbose >= 5){
        fprintf(stderr,"\t\tMt_ref_pos : %ld, md_ot_index : %ld\n",mt_ref_pos,md_ot_index);
    }
    int is_dup = 0;
    if(is_rev){
        is_dup = is_duplicate(sort_in, &sld->mt[mt_ref_pos].re.mdi);
    }
    else{
        is_dup = is_duplicate(sort_in, &sld->mt[mt_ref_pos].fe.mdi);
    }

    if(is_dup){
        return 1;
    }

    if(md_ot_index < sld->md_ot_length){
        int64_t head_ref_pos = sld->md_ot[md_ot_index]->ref_pos;
        if(sort_verbose >= 5){
            fprintf(stderr,"\t\tmd_in_ref_pos : %ld\n",md_in_ref_pos);
        }
        while(head_ref_pos <= md_in_ref_pos){
            if(sort_verbose >= 5){
                fprintf(stderr,"\t\t\tHead ref_pos : %ld\n",head_ref_pos);
            }
            if(head_ref_pos == md_in_ref_pos && is_rev == sld->md_ot[md_ot_index]->is_rev){
                is_dup = is_duplicate(sort_in, &sld->md_ot[md_ot_index]->mdi);
                if(is_dup){
                    return 1;
                }
            }
            md_ot_index++;
            if(md_ot_index >= sld->md_ot_length){
                return 0;
            }
            head_ref_pos = sld->md_ot[md_ot_index]->ref_pos;
        }
    }
    return 0;
}

void check_md_ot_head(int64_t * md_ot_head, int64_t main_ref_pos, sort_struct_t * sld){
            char in_str[2];

    if(sld->md_ot_length == 0 || *md_ot_head >= sld->md_ot_length){
        return;
    }
    if((main_ref_pos - sld->md_ot[*md_ot_head]->ref_pos < 0)){
        return;
    }
    while(main_ref_pos - sld->md_ot[*md_ot_head]->ref_pos > 101){
        *md_ot_head += 1;
        if(*md_ot_head >= sld->md_ot_length){
            if(sort_verbose >= 5){
                fprintf(stderr,"main_ref_pos : %ld, md_ot_head : %ld\n",main_ref_pos, *md_ot_head);
                fprintf(stderr,"Breaking out now\n");
                scanf("%1s",in_str);
            }
            *md_ot_head = sld->md_ot_length;
            break;
        }
        else{
            if(sort_verbose >= 5){
                fprintf(stderr,"main_ref_pos : %ld, md_ot_head : %ld, ref_pos : %ld\n",main_ref_pos, *md_ot_head, sld->md_ot[*md_ot_head]->ref_pos);
                scanf("%1s",in_str);
            }
        }
    }

    if(*md_ot_head >= sld->md_ot_length){
        // Do not allow md_ot_head to go beyond md_ot_length
        *md_ot_head = sld->md_ot_length;
    }
    return;
}

void usage(){
    printf("./bwa sort -I <input sam> -O <output sam> -v <verbose level> \n");
}


int get_sequence(char * line, size_t size, int64_t * chr_start_array, int * chr_index){
   
    if(line[1] != 'S'){
        return -1;
    }

    char * splits = NULL;
    char delim = '\t';

    chr_start_array[*chr_index] = 0;
    
    splits = strtok(line,&delim); 
    while(splits != NULL){
        if(splits[0] == 'L'){
            chr_start_array[*chr_index] = (int64_t)strtoll(splits + 3, NULL,10);
            break;
        }
        else{
            splits = strtok(NULL,&delim);
        }        
    }
    *chr_index += 1;
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


void get_record(char * line, size_t line_size, int *flags, int *chr_num, int64_t *pos, int *correction, int *avg_qual, int *mate_diff){
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
        else if(field_num == 9 && stop == 0){
            if(line[i] == '-'){
                continue;
            }
            else if((int)line[i] >= (int)'0' && (int)line[i] <= (int)'9'){
                *mate_diff = (*mate_diff * 10) + (int)line[i] - 48; 
            }
            else{
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



bseq1s_t * get_sam_record(char * line, size_t line_size,uint32_t read_id){
    
    int flags = 0;          // Field 2
    int chr_num = 0;            // Field 3
    int64_t pos = 0;            // Field 4
    int correction = 0;         // From CIGAR Field 6
    int avg_qual = 0;           // Field 11
    int mate_diff = 0;

    get_record(line, line_size, &flags, &chr_num, &pos, &correction, &avg_qual, &mate_diff);

    // TODO create the table


    if(sort_verbose >= 10){
        fprintf(stderr,"flags : %d, chr_num : %d, pos : %ld, corr : %d, mate_diff : %d\n",flags,chr_num,pos,correction, mate_diff);
    }
    if(pos == -1 || chr_num == -1){
        return NULL;
    }

    bseq1s_t * seq = (bseq1s_t *)malloc(sizeof(bseq1s_t));
    seq->read_id = read_id;
    seq->chr_num = chr_num;
    seq->abs_pos = pos;
    seq->mate_diff = mate_diff;
    seq->correction = correction;
    seq->is_rev = (flags & 0x10) != 0 ? 1 : 0;
    seq->flags = flags;
    seq->avg_qual = avg_qual;
    seq->last_entry = 0;
   
    return seq;
} 


int get_sequences(FILE * temp_unsorted_file,int64_t ** chr_start_array){
    char * line = NULL;
    size_t line_size = 0;
    size_t ret = 0;

    int chr_index = 0;
    int num_chromosomes = 0;
    int rc = 0;
    *chr_start_array = (int64_t *) malloc(24 * sizeof(int64_t));

    while(!feof(temp_unsorted_file)){
        ret = getline(&line,&line_size,temp_unsorted_file);
        if(ret != -1){
            if(line[0] == '@'){
                rc = get_sequence(line, line_size, *chr_start_array, &chr_index);
                if(rc == -1){
                    break;
                }
                else{
                    num_chromosomes++;
                }
            }
            else{
                break;
            }
        }
    }
    if(line){
        free(line);
    }
    return num_chromosomes;
}

void sort_ots(sort_struct_t * sld){
    // Sort the sort_ot 
    if(sld->sort_ot_length > 0){
        qsort(sld->sort_ot,sld->sort_ot_length,sizeof(sort_ot_entry *),sort_ot_comparator);
    }

    // Sort the md_ot 
    if(sld->md_ot_length > 0){
        qsort(sld->md_ot,sld->md_ot_length,sizeof(md_ot_entry *),md_ot_comparator);
    }
    return;
}

int mark_duplicates(sort_struct_t * sld){
    int64_t i = 0;
    int64_t head_sort_ot = 0;
    int64_t head_md_ot = 0;
    int64_t md_index = 0;
    int count_duplicates = 0; 
    for(i = 0;i<sld->mt_length;i++){
        
            char in_str[2];
        // i is the head of the main table
        // Check the head of the sort_ot and the head of the main table
        // Select one entry from (merge sort of sorts)

        /*if(head_sort_ot < sld->sort_ot_length){
            while (sld->sort_ot[head_sort_ot]->ref_pos < i){
                md_index = sld->sort_ot[head_sort_ot]->ref_pos - sld->sort_ot[head_sort_ot]->si.correction;
                if(sld->sort_ot[head_sort_ot]->is_rev){
                        check_md_ot_head(&head_md_ot, sld->sort_ot[head_sort_ot]->ref_pos, sld);
                        count_duplicates += check_duplicate(&sld->sort_ot[head_sort_ot]->si, md_index, head_md_ot, 1, sld);
                }
                else{
                        check_md_ot_head(&head_md_ot, sld->sort_ot[head_sort_ot]->ref_pos, sld);
                        count_duplicates += check_duplicate(&sld->sort_ot[head_sort_ot]->si, md_index, head_md_ot, 0,sld);
                }

                head_sort_ot++;
                if(head_sort_ot >= sld->sort_ot_length){
                    break;
                }
            }
        }*/
        check_md_ot_head(&head_md_ot, i,sld);
        /*
        if(sort_verbose >= 5){
            fprintf(stderr,"\tMd_ot_head : %ld\n",head_md_ot);
            fprintf(stderr,"\ti : %ld, Md_ot_head ref_pos: %ld\n",i,sld->md_ot[head_md_ot]->ref_pos);
        }
        */

        if(sld->mt[i].fe.si.read_id != 0){
            /*md_index = i - sld->mt[i].fe.si.correction;
            if(sort_verbose >= 5){
                fprintf(stderr,"\tForward md_index : %ld\n",md_index);
            }
            count_duplicates += check_duplicate(&sld->mt[i].fe.si, md_index, head_md_ot, 0,sld);*/
            scanf("%1s",in_str);
        }
        if(sld->mt[i].re.si.read_id != 0){
            /*md_index = i - sld->mt[i].re.si.correction;
            if(sort_verbose >= 5){
                fprintf(stderr,"\tReverse md_index : %ld\n",md_index);
            }
            count_duplicates += check_duplicate(&sld->mt[i].re.si, md_index, head_md_ot, 1,sld);*/
            scanf("%1s",in_str);
        }
    }
    return count_duplicates;
}

void sort_ST(void * data){
    struct timespec proc_start, proc_end;
    sort_slave_t * s = (sort_slave_t *) data;
    s->time = 0.0;
    s->num_reads = 0;
    sort_struct_t * sld = (sort_struct_t *) malloc(sizeof(sort_struct_t));

    queue * q = s->q;
    struct_init(sld,s->length);
    void ** queue_out = (void **) malloc(sizeof(void *)) ;
    bseq1s_t * seq;

    while(1){
        getElement(q,queue_out); 
        seq = (bseq1s_t *) *queue_out;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID,&proc_start);
        if(seq->last_entry){
            break;
        }
        else{
            process_sam_record_1(seq, sld);
        }
        if(seq){
            free(seq);
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID,&proc_end);
        s->time += (double)(proc_end.tv_sec - proc_start.tv_sec) + ((double)(proc_end.tv_nsec - proc_start.tv_nsec)/(double)(1000000000));
        s->num_reads++;
    }

    sort_ots(sld);
    if(sort_verbose >= 15){
        fprintf(stderr,"md_ot_length : %d\n",sld->md_ot_length);
    }
    if(s->thread_id == 0){
        mark_duplicates(sld);
    }

    struct_delete(sld);
    free(sld);
    free(queue_out);
    pthread_exit(0);
}

void sort_MT(FILE * in_sam, int num_threads){
    int64_t * chr_len;
    char * line = NULL;
    size_t line_size = 0;
    size_t ret = 0;

    int num_chromosomes = 0;
    int i = 0;
    uint32_t read_id = 1;

    if(sort_verbose >= 3){
        clock_gettime(CLOCK_REALTIME,&ini_start);
    }
    num_chromosomes = get_sequences(in_sam, &chr_len);

    assert(num_chromosomes == 24);

    // Create threads here

    int factor = (int)ceil((double)num_chromosomes / (double)num_threads);
    int * chr_thread_id = (int *)malloc(num_chromosomes * sizeof(int));
    int64_t * lens = (int64_t *) malloc(num_threads * sizeof(int64_t));
    memset(lens,0,num_threads * sizeof(int64_t));


    int tid = 0;
    int64_t running_length = 0;
    int64_t temp = 0;
    int prev_tid = -1;
    for(i = 0;i<num_chromosomes;i++){
        tid = (i / factor);
        chr_thread_id[i] = tid;
        if(prev_tid != tid){
            if(prev_tid != -1){
                lens[prev_tid] = running_length;
            }
            running_length = 0;
        }
        temp = chr_len[i];
        chr_len[i] = running_length;
        running_length += temp;
        prev_tid = tid;
    }

    lens[prev_tid] = running_length;

    if(sort_verbose >= 5){
        fprintf(stderr,"Num chromosomes : %d\n",num_chromosomes);
        fprintf(stderr,"Factor : %d\n",factor);
        int j = 0;
        for(j=0;j<num_chromosomes;j++){
            fprintf(stderr,"tid[%d] : %d , %ld\n",j+1,chr_thread_id[j],chr_len[j]);
        }
        for(j=0;j<num_threads;j++){
            fprintf(stderr,"lens[%d] : %ld\n",j+1,lens[j]);
        }
    }

    
    queue ** qs = (queue **) malloc(num_threads * sizeof(queue *));
    pthread_t * pts = (pthread_t *) malloc(num_threads * sizeof(pthread_t));
    sort_slave_t * slaves = (sort_slave_t *) malloc(num_threads * sizeof(sort_slave_t));
    bseq1s_t * s;

    // Launching Threads
    if(sort_verbose >= 3){
        fprintf(stderr,"Launching %d threads\n",num_threads);
    }
    for(i=0;i<num_threads;i++){
        qs[i] = queueInit(QUEUE_SIZE);
        slaves[i].thread_id = i;
        slaves[i].length = lens[i];
        slaves[i].q = qs[i];
        slaves[i].time = 0;
        slaves[i].num_reads = 0;
        pthread_create(&pts[i],NULL,sort_ST,&slaves[i]);
    } 

    if(sort_verbose >= 3){
        clock_gettime(CLOCK_REALTIME,&ini_end);
        ini_time = (double)(ini_end.tv_sec - ini_start.tv_sec) + ((double)(ini_end.tv_nsec - ini_start.tv_nsec)/(double)(1000000000));
        fprintf(stderr,"Initialization done in %f secs\n",ini_time);
    }


    queue * sending_queue = NULL;
    if(sort_verbose >= 3){
        clock_gettime(CLOCK_REALTIME,&read_start);
    }
    while(!feof(in_sam)){
        ret = getline(&line,&line_size,in_sam);
        if(ret != -1 && line[0] != '@'){
            s = get_sam_record(line, strlen(line), read_id);
            read_id++;
            if(s != NULL){
                s->abs_pos += s->chr_num;
                tid = chr_thread_id[s->chr_num - 1];
                sending_queue = qs[tid];
                addElement(sending_queue, (void *)s);
            }
            if(sort_verbose >= 3 && read_id % 1000000 == 0){
                clock_gettime(CLOCK_REALTIME,&read_end);
                read_time = (double)(read_end.tv_sec - read_start.tv_sec) + ((double)(read_end.tv_nsec - read_start.tv_nsec)/(double)(1000000000));
                fprintf(stderr,"Read %d reads in %f secs\n",read_id,read_time);
            }
        }
    }


    for(i=0;i<num_threads;i++){
        sending_queue = qs[i];
        s = (bseq1s_t *) malloc( sizeof(bseq1s_t));
        s->read_id = 0;
        s->chr_num = 0;
        s->abs_pos = 0;
        s->mate_diff = 0;
        s->correction = 0;
        s->is_rev = 0;
        s->flags = 0;
        s->avg_qual = 0;
        s->last_entry = 1;
        addElement(sending_queue, (void *)s);
    }

    if(sort_verbose >= 3){
        clock_gettime(CLOCK_REALTIME,&read_end);
        read_time = (double)(read_end.tv_sec - read_start.tv_sec) + ((double)(read_end.tv_nsec - read_start.tv_nsec)/(double)(1000000000));
        fprintf(stderr,"Read %d reads in %f secs\n",read_id-1,read_time);
    }


    if(sort_verbose >= 3){
        clock_gettime(CLOCK_REALTIME,&finish_start);
        fprintf(stderr,"Waiting for threads\n");
    }


    for(i=0;i<num_threads;i++){
        pthread_join(pts[i],NULL);
    }

    if(sort_verbose >= 3){
        clock_gettime(CLOCK_REALTIME,&finish_end);
        finish_time = (double)(finish_end.tv_sec - finish_start.tv_sec) + ((double)(finish_end.tv_nsec - finish_start.tv_nsec)/(double)(1000000000));
        fprintf(stderr,"Finished all threads in %f secs\n",finish_time);
    }


    for(i=0;i<num_threads;i++){
        queueDelete(qs[i]);
    }
    
    for(i=0;i<num_threads;i++){
        fprintf(stderr,"%ld Reads processed in tid %ld in %f secs\n",slaves[i].num_reads,i,slaves[i].time);
    }

    free(slaves);
    free(pts);
    free(qs);

    free(lens);
    free(chr_thread_id);
    free(chr_len);
    return;
}

int main_memsort(int argc, char *argv[]){
    
    int rargc = 2;
    int i = 0;
    int output_file = 0; 
    int num_threads = 1;
    FILE * temp_unsorted_file;
    char * temp_unsorted_file_name;
    if(argc < rargc){
        usage();
        return 0;
    }
   
    for(i=0;i<argc;i+=2){
        if(strcmp(argv[i],"-I") == 0) {
            if(argc > i+1){
                temp_unsorted_file_name = argv[i+1];
            }
            else{
                usage();
                return -1;
            }
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
        if(strcmp(argv[i],"-t") == 0){
            if(argc > i+1){
                num_threads = (int)strtol(argv[i+1],NULL,10);
                if(num_threads > 24){
                    fprintf(stderr,"Max number of threads allowed is 24, shifting to 24 threads");
                    num_threads = 24;
                }
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

    sort_MT(temp_unsorted_file, num_threads);

    fclose(temp_unsorted_file);
    return 0;

}
