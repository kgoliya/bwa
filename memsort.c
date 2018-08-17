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

    sld->mt = (mt_entry * ) malloc(length * sizeof(mt_entry));
    memset(sld->mt,0,(length * sizeof(mt_entry)));
    
    sld->sort_ot = (ot_entry ** ) malloc(INIT_OVERFLOW_TABLE_SIZE * sizeof(ot_entry *));
    memset(sld->sort_ot,0,(INIT_OVERFLOW_TABLE_SIZE * sizeof(ot_entry *)));

}

void struct_delete(sort_struct_t * sld){

    free(sld->mt);
    int i = 0;
    for(i = sld->sort_ot_length-1;i>=0;i--){
        free(sld->sort_ot[i]);
    }
    free(sld->sort_ot);
}


void write_entry(bseq1s_t * seq, half_mt_entry * in_en){
    in_en->read_id = seq->read_id;
    in_en->correction = seq->correction;
    in_en->avg_qual = seq->avg_qual;
    in_en->flags = seq->flags;
    in_en->mate_diff = seq->mate_diff;
    return;
}

void copy_entry(half_mt_entry * src_en, half_mt_entry * dest_en){
    dest_en->read_id = src_en->read_id;
    dest_en->correction = src_en->correction;
    dest_en->avg_qual = src_en->avg_qual;
    dest_en->flags = src_en->flags;
    dest_en->mate_diff = src_en->mate_diff;
    return; 
}

void add_sort_ot_entry(bseq1s_t * seq,sort_struct_t * sld){

    if(sld->sort_ot_length != 0 && sld->sort_ot_length % sld->sort_ot_size == 0){
        sld->sort_ot = realloc(sld->sort_ot,(sld->sort_ot_length + sld->sort_ot_size) * sizeof(ot_entry *));
    }
   
    int index = sld->sort_ot_length;

    sld->sort_ot[index] = (ot_entry * ) malloc(sizeof(ot_entry));
    sld->sort_ot[index]->ref_pos = seq->abs_pos;

    write_entry(seq,&sld->sort_ot[index]->ote);
    sld->sort_ot_length++;
}

void add_sort_ot_entry_from_mt(int64_t ref_pos,half_mt_entry * in_en, sort_struct_t * sld){
    if(sld->sort_ot_length != 0 && sld->sort_ot_length % sld->sort_ot_size == 0){
        sld->sort_ot = realloc(sld->sort_ot,(sld->sort_ot_length + sld->sort_ot_size) * sizeof(ot_entry *));
    }
   
    int index = sld->sort_ot_length;

    sld->sort_ot[index] = (ot_entry * ) malloc(sizeof(ot_entry));
    sld->sort_ot[index]->ref_pos = ref_pos;
    
    copy_entry(in_en, &sld->sort_ot[index]->ote);
    
    sld->sort_ot_length++;
}


void add_mt_entry_1(bseq1s_t *seq, sort_struct_t * sld){


    half_mt_entry * sort_en;
    if(sort_verbose >= 10){
        fprintf(stderr,"Ref pos : %ld\n",seq->abs_pos);
        fprintf(stderr,"Flags : %d\n",seq->flags);
        fprintf(stderr,"Mate pos : %d\n",seq->mate_diff);
    }
    if(seq->is_rev){
        sort_en = &sld->mt[seq->abs_pos].re;
    }
    else {
        sort_en = &sld->mt[seq->abs_pos].fe;
    }

    if(sort_en->flags == 0 && ((seq->flags & 0x100) == 0)){
        //Sort entry was free
        write_entry(seq,sort_en);
    }
    else if((seq->flags & 0x100) != 0){
        add_sort_ot_entry(seq,sld);
    }
    else{
        // Add entry to Sort Overflow Table
        // Check if both the entries are properly paired
        if(((sort_en->flags & 0x02) != 0) && ((seq->flags & 0x02) != 0)){
            // Both entries are properly paired
            if(sort_en->mate_diff == seq->mate_diff){
                if(sort_en->avg_qual < seq->avg_qual){
                    sort_en->flags |= 0x400;            // Duplicate marked
                    sort_en->flags |= 0x8000;           // Entry processed
                    add_sort_ot_entry_from_mt(seq->abs_pos,sort_en, sld);
                    write_entry(seq,sort_en);
                }
                else{
                    // Seq is the duplicate entry
                    seq->flags |= 0x400;
                    seq->flags |= 0x8000;
                    add_sort_ot_entry(seq,sld);
                }
            }
            else{
                add_sort_ot_entry(seq,sld);
            }
        }
        else if(((sort_en->flags & 0x02) == 0) && ((seq->flags & 0x02) != 0)){
            // If sort_en is not properly paired and seq is properly paired
            // Put sort_en in overflow table
            if(sort_en->avg_qual <= seq->avg_qual){
                sort_en->flags |= 0x400; 
                sort_en->flags |= 0x8000;
            }
            add_sort_ot_entry_from_mt(seq->abs_pos,sort_en, sld);
            write_entry(seq,sort_en);
        }
        else if(((seq->flags & 0x02) == 0) && ((sort_en->flags & 0x02) != 0)){
             if(seq->avg_qual <= sort_en->avg_qual){
                 seq->flags |= 0x400;
                 seq->flags |= 0x8000;
             }
             add_sort_ot_entry(seq,sld);
        }
        else{
            // Both are not properly paired
            if(seq->avg_qual <= sort_en->avg_qual){
                 seq->flags |= 0x400;
                 seq->flags |= 0x8000;
                 add_sort_ot_entry(seq,sld);
            }
            else{
                sort_en->flags |= 0x400;            // Duplicate marked
                sort_en->flags |= 0x8000;           // Entry processed
                add_sort_ot_entry_from_mt(seq->abs_pos,sort_en, sld);
                write_entry(seq,sort_en);
            }
        }
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
    ot_entry* l = *((ot_entry **)p);
    ot_entry* r = *((ot_entry **)q); 
    
    if(l->ref_pos < r->ref_pos){
        return -1;
    }
    else if(l->ref_pos == r->ref_pos){
        if((((l->ote.flags & 0x100) != 0) || ((l->ote.flags & 0x400) != 0)) && ((r->ote.flags & 0x100) == 0) && ((r->ote.flags & 0x400) == 0)){
            // r should go before l
            return 1;
        }
        else if((((r->ote.flags & 0x100) != 0) || ((r->ote.flags & 0x400) != 0)) && ((l->ote.flags & 0x100) == 0) && ((l->ote.flags & 0x400) == 0)){
            // l should go before r
            return -1;
        }
        else{
            return 0;
        }
    }
    else{
        return 1;
    }
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


    if(pos == -1 || chr_num == -1){
        return NULL;
    }

    bseq1s_t * seq = (bseq1s_t *)malloc(sizeof(bseq1s_t));
    seq->read_id = read_id;
    seq->chr_num = chr_num;
    seq->abs_pos = pos - correction;
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
        qsort(sld->sort_ot,sld->sort_ot_length,sizeof(ot_entry *),sort_ot_comparator);
    }

    return;
}




// 0 : No decision taken
// 1 : en1 is duplicate
// -1 : en2 is duplicate
int is_duplicate(half_mt_entry * en1, half_mt_entry * en2){
   
    if(((en1->flags & 0x400) != 0) || ((en2->flags & 0x400) != 0) || ((en1->flags & 0x100) != 0) || ((en1->flags & 0x100) != 0)){
        // either entry is a duplicate or a secondary alignment
        return 0;
    }
    if((en1->flags & 0x10) == (en2->flags & 0x10)){
        // Both entries have the same orientation

        if((en1->flags & 0x02) == 0){
            // Not properly paired
            if(en1->avg_qual <= en2->avg_qual){
                return 1;
            }
            else{
                if((en2->flags & 0x02) == 0){
                    // if en2 is also not paired properly
                    return -1;
                }
            }
        }
        else{
            if((en2->flags & 0x02) != 0){
                if(en1->mate_diff == en2->mate_diff){
                    if(en1->avg_qual <= en2->avg_qual){
                        return 1;
                    }
                    else{
                        return -1;
                    }
                }
            }
        }
    }
    return 0;

}



void check_for_duplicates(int64_t ref_pos,half_mt_entry * in_en, int64_t md_head, sort_struct_t * sld, int *num_duplicates){
    assert (in_en != NULL);
    if(md_head >= sld->sort_ot_length || ((in_en->flags & 0x8000) != 0) || ((in_en->flags & 0x100) != 0)){
        // Ignore input if input has been processed already or its a secondary alignment
        if((in_en->flags & 0x400) != 0){
            *num_duplicates = *num_duplicates + 1;
        }
        return;
    }
    int dup = 0;
    while(sld->sort_ot[md_head]->ref_pos <= ref_pos){
        if(sld->sort_ot[md_head]->ref_pos == ref_pos){
            // The entry has the same ref pos
            if(((sld->sort_ot[md_head]->ote.flags & 0x100) != 0) || ((sld->sort_ot[md_head]->ote.flags & 0x400) != 0)){
                break;
            }
            if(sld->sort_ot[md_head]->ote.read_id != in_en->read_id){
                // Not the same read
                dup = is_duplicate(in_en,&sld->sort_ot[md_head]->ote);
                if(dup == 1){
                    in_en->flags |= 0x400;
                    *num_duplicates = *num_duplicates + 1;
                    break;
                }
                else if(dup == -1){
                    *num_duplicates = *num_duplicates + 1;
                    sld->sort_ot[md_head]->ote.flags |= 0x400;
                    sld->sort_ot[md_head]->ote.flags |= 0x8000;
                }
            }
        }
        md_head++;
        if(md_head >= sld->sort_ot_length){
            break;
        }
    }

    in_en->flags |= 0x8000;
    return;
}

/*
void mark_duplicate_mt_entry(int64_t ref_pos,int64_t md_head,half_mt_entry * in_en, sort_struct_t * sld, int check_for_duplicate){
    int64_t i = 0;
    int64_t ref_pos_without_corr = (ref_pos + in_en->correction);
    half_mt_entry * single_en;
    int64_t single_en_uncorr_ref_pos = 0;
    assert(ref_pos_without_corr < sld->mt_length);
    int64_t tmp_md_head = md_head;
   
    if(in_en->read_id == 0 || ((in_en->flags & 0x8000) != 0)){
        // If in_en is invalid or it has already been seen before
        return;
    }

    for(i = ref_pos+1;i<ref_pos_without_corr;i++){
        if(md_head < sld->sort_ot_length){
            // Get entry in the overflow table
            single_en = &sld->sort_ot[md_head]->ote;
            single_en_uncorr_ref_pos = sld->sort_ot[md_head]->ref_pos + single_en->correction;

            while(ref_pos_without_corr > single_en_uncorr_ref_pos){
                // The entry in the overflow table has a sort position smaller 
                // the current entry under consideration
                mark_duplicate_mt_entry(sld->sort_ot[md_head]->ref_pos, md_head, single_en, sld,1);
                md_head++;
                if(md_head < sld->sort_ot_length){
                    single_en = &sld->sort_ot[md_head]->ote;
                    single_en_uncorr_ref_pos = sld->sort_ot[md_head]->ref_pos + single_en->correction;
                }
                else{
                    break;
                }
            }
        }
        mark_duplicate_mt_entry(sld->sort_ot[md_head]->ref_pos, md_head, single_en, sld, 0);
    }
    if(check_for_duplicate != 0){
        check_for_duplicates(ref_pos,in_en, tmp_md_head,sld);
    }
    return;
}

*/

void update_ot_head(int64_t * ot_head, int64_t current_ref,sort_struct_t * sld){
    if(sld->sort_ot_length == 0 || sld->sort_ot_length <= *ot_head){
        
        // If sort_ot_length is zero or the head is already past the full table
        // do nothing
        return;
    }

    if(current_ref < sld->sort_ot[*ot_head]->ref_pos){
        
        // If current ref positon is less than the ref pos
        // at the head of the table, do nothing
        return;
    }

    while(current_ref > sld->sort_ot[*ot_head]->ref_pos){
        *ot_head = *ot_head + 1;
        if(*ot_head >= sld->sort_ot_length){
            *ot_head = sld->sort_ot_length;
            break;
        }
    }
    return;
}


void mark_duplicates(sort_struct_t * sld, int *num_duplicates){
    int64_t i = 0;
    int64_t md_head = 0;
    int64_t prev_ref_pos = 0;

    half_mt_entry * single_en;
    for(i=0;i<sld->sort_ot_length;i++){

        single_en = &sld->sort_ot[i]->ote;
        if(prev_ref_pos != sld->sort_ot[i]->ref_pos){
            md_head = i;
            prev_ref_pos = sld->sort_ot[i]->ref_pos;
        }

        if(sort_verbose >= 5){
            fprintf(stderr,"i : %ld,md_head : %ld, ref_pos : %ld\n",i,md_head,sld->sort_ot[md_head]->ref_pos);
        }
        check_for_duplicates(i,single_en, md_head, sld, num_duplicates);
    }

    return;
}

void sort_ST(void * data){
    struct timespec proc_start, proc_end;
    sort_slave_t * s = (sort_slave_t *) data;
    s->write_time = 0.0;
    s->num_reads = 0;
    sort_struct_t * sld = (sort_struct_t *) malloc(sizeof(sort_struct_t));

    queue * q = s->q;
    struct_init(sld,s->length);
    void ** queue_out = (void **) malloc(sizeof(void *)) ;
    bseq1s_t * seq;

    // Write Portion
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
        s->write_time += (double)(proc_end.tv_sec - proc_start.tv_sec) + ((double)(proc_end.tv_nsec - proc_start.tv_nsec)/(double)(1000000000));
        s->num_reads++;
    }
    // Write Portion End

    // Sorting the overflow table
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&proc_start);

    sort_ots(sld);
    
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&proc_end);
    s->sort_time = (double)(proc_end.tv_sec - proc_start.tv_sec) + ((double)(proc_end.tv_nsec - proc_start.tv_nsec)/(double)(1000000000));
    // Sorting end

    // Mark duplicates in the overflow table
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&proc_start);
    s->num_duplicates = 0;
    mark_duplicates(sld,&s->num_duplicates);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&proc_end);
    s->md_time = (double)(proc_end.tv_sec - proc_start.tv_sec) + ((double)(proc_end.tv_nsec - proc_start.tv_nsec)/(double)(1000000000));
    
    // Mark duplicates end

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
        slaves[i].num_duplicates = 0;

        slaves[i].write_time = 0.0;
        slaves[i].sort_time = 0.0;
        slaves[i].md_time = 0.0;
        slaves[i].num_reads = 0;

        slaves[i].q = qs[i];
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
                s->abs_pos += chr_len[s->chr_num - 1];
                tid = chr_thread_id[s->chr_num - 1];
                sending_queue = qs[tid];
                addElement(sending_queue, (void *)s);
            }
            if(sort_verbose >= 3 && read_id % 10000000 == 0){
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
    int total_duplicates = 0;
    if(sort_verbose >= 3){
        for(i=0;i<num_threads;i++){
            fprintf(stderr,"For tid : %d\n",i+1);
            fprintf(stderr,"\t%ld Reads processed in %f secs\n",slaves[i].num_reads,slaves[i].write_time);
            fprintf(stderr,"\t%d Duplicates detected\n",slaves[i].num_duplicates);
            total_duplicates += slaves[i].num_duplicates;

            fprintf(stderr,"\tSorting time : %f\n",slaves[i].sort_time);
            fprintf(stderr,"\tMark duplicates time : %f\n",slaves[i].md_time);
        }
    }

    fprintf(stderr,"Total duplicates detected : %d\n",total_duplicates);

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
