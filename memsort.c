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

int sort_verbose = 0;

struct timespec ini_start, ini_end;
double ini_time = 0.0;

struct timespec read_start, read_end;
double read_time = 0.0;

struct timespec finish_start, finish_end;
double finish_time = 0.0;


void struct_init(sort_struct_t * sld, int64_t length){            

    sld->mt_length = length;

    sld->fot_size = INIT_OVERFLOW_TABLE_SIZE;
    sld->fot_length = 0;

    sld->rot_size = INIT_OVERFLOW_TABLE_SIZE;
    sld->rot_length = 0;

    sld->fmt = (half_mt_entry ** ) malloc(length * sizeof(half_mt_entry *));
    memset(sld->fmt,0,(length * sizeof(half_mt_entry *)));
   
    sld->rmt = (half_mt_entry ** ) malloc(length * sizeof(half_mt_entry *));
    memset(sld->rmt,0,(length * sizeof(half_mt_entry *)));


    sld->fot = (ot_entry ** ) malloc(INIT_OVERFLOW_TABLE_SIZE * sizeof(ot_entry *));
    memset(sld->fot,0,(INIT_OVERFLOW_TABLE_SIZE * sizeof(ot_entry *)));

    sld->rot = (ot_entry ** ) malloc(INIT_OVERFLOW_TABLE_SIZE * sizeof(ot_entry *));
    memset(sld->rot,0,(INIT_OVERFLOW_TABLE_SIZE * sizeof(ot_entry *)));

}

void struct_delete(sort_struct_t * sld){

    
    int64_t i = 0;
    for(i = 0;i<sld->mt_length;i++){
        if(sld->fmt[i]){
            free(sld->fmt[i]);
        }
        if(sld->rmt[i]){
            free(sld->rmt[i]);
        }
    }
    free(sld->fmt);
    free(sld->rmt);
    for(i = sld->fot_length-1;i>=0;i--){
        free(sld->fot[i]);
    }
    free(sld->fot);

    for(i = sld->rot_length-1;i>=0;i--){
        free(sld->rot[i]);
    }
    free(sld->rot);
}

void copy_fileptr(uint8_t * src_fileptr, uint8_t * dest_fileptr){
    dest_fileptr[0] = src_fileptr[0];
    dest_fileptr[1] = src_fileptr[1];
    dest_fileptr[2] = src_fileptr[2];
    dest_fileptr[3] = src_fileptr[3];
    dest_fileptr[4] = src_fileptr[4];
}

void write_entry(bseq1s_t * seq, half_mt_entry * in_en){
    copy_fileptr(seq->fileptr, in_en->fileptr);
    in_en->sam_size = seq->sam_size;
    in_en->correction = seq->correction;
    in_en->avg_qual = seq->avg_qual;
    in_en->flags = seq->flags;
    in_en->mate_diff = seq->mate_diff;
    in_en->ref_pos = seq->abs_pos;
    in_en->mate_chr_num = seq->mate_chr_num;
    in_en->mate_pos = seq->mate_pos;
    return;
}

void copy_entry(half_mt_entry * src_en, half_mt_entry * dest_en){

    copy_fileptr(src_en->fileptr, dest_en->fileptr);
    dest_en->sam_size = src_en->sam_size;
    dest_en->correction = src_en->correction;
    dest_en->avg_qual = src_en->avg_qual;
    dest_en->flags = src_en->flags;
    dest_en->mate_diff = src_en->mate_diff;
    dest_en->ref_pos = src_en->ref_pos;
    dest_en->mate = src_en->mate;
    dest_en->mate_chr_num = src_en->mate_chr_num;
    dest_en->mate_pos = src_en->mate_pos;

    // Also go back to the mates position and change it to reflect dest pointer and not source ptr
    return; 
}

half_mt_entry * add_ot_entry(bseq1s_t * seq,sort_struct_t * sld){

    int is_rev = (seq->flags & 0x10);
    int index = 0;

    // Add ot entry : Return the address of sld->fot[index]->ote / sld->rot[index]->ote


    if(is_rev == 0){
        // seq belongs to forward strand
        if(sld->fot_length != 0 && sld->fot_length % sld->fot_size == 0){
            sld->fot = realloc(sld->fot,(sld->fot_length + sld->fot_size) * sizeof(ot_entry *));
        }

        index = sld->fot_length;

        sld->fot[index] = (ot_entry * ) malloc(sizeof(ot_entry));
        //sld->fot[index]->ote = (half_mt_entry *) malloc(sizeof(half_mt_entry));
        sld->fot[index]->ref_pos = seq->abs_pos;

        write_entry(seq,&sld->fot[index]->ote);
        sld->fot_length = sld->fot_length + 1;
        return &sld->fot[index]->ote;
    }
    else{
        // seq belongs to reverse strand
        if(sld->rot_length != 0 && sld->rot_length % sld->rot_size == 0){
            sld->rot = realloc(sld->rot,(sld->rot_length + sld->rot_size) * sizeof(ot_entry *));
        }

        index = sld->rot_length;

        sld->rot[index] = (ot_entry * ) malloc(sizeof(ot_entry));
        sld->rot[index]->ref_pos = seq->abs_pos;

        write_entry(seq,&sld->rot[index]->ote);
        sld->rot_length = sld->rot_length + 1;
        return &sld->rot[index]->ote;
    }

}

void add_ot_entry_from_mt(int64_t ref_pos,half_mt_entry * in_en, sort_struct_t * sld){
    int is_rev = (in_en->flags & 0x10);
    int index = 0;

    // If we ever enter this function, that means an entry from the main table is being transferred to the 
    // overflow table. Thus we will need to change the location of the mate's mate (i.e. the current position if it
    // is proerly paired and if it is not a secondary alignment

    half_mt_entry * new_en = NULL;
    new_en = (half_mt_entry *)in_en->mate;
    new_en = NULL;

    if(is_rev == 0){
        if(sld->fot_length != 0 && sld->fot_length % sld->fot_size == 0){
            sld->fot = realloc(sld->fot,(sld->fot_length + sld->fot_size) * sizeof(ot_entry *));
        }

        index = sld->fot_length;

        sld->fot[index] = (ot_entry * ) malloc(sizeof(ot_entry));
        sld->fot[index]->ref_pos = ref_pos;

        copy_entry(in_en, &sld->fot[index]->ote);
        new_en = &sld->fot[index]->ote;

        sld->fot_length = sld->fot_length + 1;
    }
    else{
        if(sld->rot_length != 0 && sld->rot_length % sld->rot_size == 0){
            sld->rot = realloc(sld->rot,(sld->rot_length + sld->rot_size) * sizeof(ot_entry *));
        }

        index = sld->rot_length;

        sld->rot[index] = (ot_entry * ) malloc(sizeof(ot_entry));
        sld->rot[index]->ref_pos = ref_pos;

        copy_entry(in_en, &sld->rot[index]->ote);
        new_en = &sld->rot[index]->ote;

        sld->rot_length = sld->rot_length + 1;
    }

    //check if new_en is properly paired, if yes then update its mate's mate 
    if(((new_en->flags & 0x100) == 0) && ((new_en->flags & 0x2) != 0)){
        // new_en is not a secondary alignment and is properly paired
        half_mt_entry * new_en_m = (half_mt_entry *)new_en->mate;
        new_en_m->mate = (void *)new_en;
    }
}


half_mt_entry * add_mt_entry(bseq1s_t *seq, sort_struct_t * sld, half_mt_entry * mate){


    half_mt_entry ** sort_en_add;
    half_mt_entry * sort_en;
    
    
    int is_rev = seq->is_rev;
    int dup = 0;
    
    if(seq->is_rev){
        sort_en_add = &sld->rmt[seq->abs_pos];
    }
    else {
        sort_en_add = &sld->fmt[seq->abs_pos];
    }

    

    if(*sort_en_add == NULL){
        *sort_en_add = (half_mt_entry *) malloc(sizeof(half_mt_entry));
        memset(*sort_en_add,0,sizeof(half_mt_entry));
        (*sort_en_add)->flags = 0; 
    }


    sort_en = *sort_en_add;


    // Process seq first then worry about seq->mate

    if(sort_en->flags == 0){
        //Sort entry was free
        //mate value will be NULL if mate is not present
        write_entry(seq,sort_en);
    }
    else{
        // Add entry to Sort Overflow Table
        // Check if both the entries are properly paired

        if(((sort_en->flags & 0x100) != 0) || ((seq->flags & 0x100) != 0)){
            //Sort_en or seq is a secondary alignment 

            if(((seq->flags & 0x100) != 0)){
                // Seq is are secondary alignments
                // Place seq in overflow table
                sort_en = add_ot_entry(seq,sld);

                // If seq is the entry being places in the overflow table, and sort_en is already valid, 
                // this means that the mate is not going to look at old sort_en anymore, thus change sort_en to 
                // the entry just written into the overflow table
            }
            else {
                // Seq is not a secondary alignment, which means sort_en is a secondary alignment
                add_ot_entry_from_mt(seq->abs_pos,sort_en, sld);
                write_entry(seq,sort_en);
            }

        }
        else if(((sort_en->flags & 0x4) != 0) || ((seq->flags & 0x4) != 0)){
            //Sort_en or seq is an unmapped alignment 
                
            if(((seq->flags & 0x4) != 0)){
                // Seq is an unmapped alignments
                // Place seq in overflow table
                sort_en = add_ot_entry(seq,sld);
            }
            else {
                // Seq is not an unmapped alignment, which means sort_en is an unmapped alignment
                add_ot_entry_from_mt(seq->abs_pos,sort_en, sld);
                write_entry(seq,sort_en);
            }

        }
        else {

            // Both entries are not secondary or unmapped alignments

            // Check if any is a unpaired read, 

            if(((sort_en->flags & 0x2) != 0) && ((seq->flags & 0x2) != 0)){
                // both are either properly paired thus they can be compared
                if(sort_en == mate){
                    dup = 0;
                }
                else{
                    dup = is_duplicate_seq_en(seq, sort_en,is_rev, seq->abs_pos);
                }

                if(dup == 1){
                    // seq is the duplicate
                    // place seq in the overflow table
                    sort_en = add_ot_entry(seq,sld);
                }
                else if(dup == -1){
                    // sort_en is the duplicate, 
                    // place sort_en in overflow table
                    add_ot_entry_from_mt(seq->abs_pos,sort_en, sld);
                    write_entry(seq,sort_en);
                }
                else if(dup == 2){
                    // seq or sort_en was a duplicate previously
                    if((sort_en->flags & 0x400) != 0){
                        // sort_en was the duplicate, 
                        // place sort_en in overflow table
                        add_ot_entry_from_mt(seq->abs_pos,sort_en, sld);
                        write_entry(seq,sort_en);
                    }
                    else{
                        // seq was marked as duplicate
                        // place seq in overflow table
                        sort_en = add_ot_entry(seq,sld);
                    }
                }
                else{
                    // No decision taken, place seq in overflow table
                    sort_en = add_ot_entry(seq,sld);
                }

            }
            else{
                // one is properly paired and the other is not properly paired
                if((sort_en->flags & 0x2) == 0){
                    // sort_en is the unpaired entry
                    // Place seq in the ot table
                    sort_en = add_ot_entry(seq,sld);
                }
                else{
                    // seq is the unpaired entry
                    add_ot_entry_from_mt(seq->abs_pos,sort_en, sld);
                    write_entry(seq,sort_en);
                }
            }





        }
    }

    return sort_en;
}


void process_sam_record_1(bseq1s_t *seq, sort_struct_t * sld){

    half_mt_entry * seq_en;
    half_mt_entry * mate_en;
    if(seq == NULL){
        return;
    }
    seq_en = add_mt_entry(seq, sld,NULL);

    // Check if the seq has a mate
    if(seq->mate != NULL){
        mate_en = add_mt_entry((bseq1s_t *)seq->mate, sld,seq_en);
        seq_en->mate = (void *) mate_en;
        mate_en->mate = (void *) seq_en;
        
    }
    else{
        seq_en->mate = NULL;
    }

}

int lexicographic_comparator(char * str0, int len0, char * str1, int len1){
    int smaller_len_index = len0 < len1 ? 0 : 1;
    int smaller_len = (smaller_len_index == 0) ? len0 : len1;

    int i = 0;
    int decision = 0;
    for(i=0;i<smaller_len;i++){
        if(str0[i] != str1[i]){
            if(str0[i] < str1[i]){
                decision = -1;
            }
            else{
                decision = 1;
            }
            break;
        }
    }

    if(decision == 0){
        // No decision could be made
        if(len0 != len1){
            // Lengths are not equal
            // The one with the smaller length should be first
            if(len0 < len1){
                // len0 should go before len1
                decision = -1;
            }
            else{
                decision = 1;
            }
        }
        // else no decision can be reached
    }

    return decision;
}

void print_mt_entry(half_mt_entry * en1, half_mt_entry * en2){
    int64_t f1 = 0,f2 = 0;
    memcpy(&f1,en1->fileptr,5*sizeof(uint8_t));
    memcpy(&f2,en2->fileptr,5*sizeof(uint8_t));
    fprintf(stderr,"Ref_pos : %ld\t|\tRef_pos : %ld\n",en1->ref_pos,en2->ref_pos);
    fprintf(stderr,"Fileptr : %ld\t\t|\tFileptr : %ld\n",f1,f2);
    fprintf(stderr,"Corr : %d\t\t|\tCorr : %d\n",en1->correction,en2->correction);
    fprintf(stderr,"Avg_qual : %d\t\t|\tAvg_qual : %d\n",en1->avg_qual,en2->avg_qual);
    fprintf(stderr,"Flags : %d\t\t|\tFlags : %d\n",(en1->flags & 0xfff),(en2->flags & 0xfff));
    fprintf(stderr,"Mate diff : %d\t\t|\tMate diff : %d\n",en1->mate_diff,en2->mate_diff);
    fprintf(stderr,"Mate chr num : %d\t\t|\tMate chr num : %d\n",en1->mate_chr_num,en2->mate_chr_num);
    fprintf(stderr,"Mate pos : %d\t\t|\tMate pos : %d\n",en1->mate_pos,en2->mate_pos);

}

void print_seq_mt_entry(bseq1s_t * en1, half_mt_entry * en2){
    int64_t f1 = 0,f2 = 0;
    memcpy(&f1,en1->fileptr,5*sizeof(uint8_t));
    memcpy(&f2,en2->fileptr,5*sizeof(uint8_t));
    fprintf(stderr,"Ref_pos : %ld\t|\tRef_pos : %ld\n",en1->abs_pos,en2->ref_pos);
    fprintf(stderr,"Fileptr : %ld\t\t|\tFileptr : %ld\n",f1,f2);
    fprintf(stderr,"Corr : %d\t\t|\tCorr : %d\n",en1->correction,en2->correction);
    fprintf(stderr,"Avg_qual : %d\t\t|\tAvg_qual : %d\n",en1->avg_qual,en2->avg_qual);
    fprintf(stderr,"Flags : %d\t\t|\tFlags : %d\n",(en1->flags & 0xfff),(en2->flags & 0xfff));
    fprintf(stderr,"Mate diff : %d\t\t|\tMate diff : %d\n",en1->mate_diff,en2->mate_diff);
    fprintf(stderr,"Mate chr num : %d\t\t|\tMate chr num : %d\n",en1->mate_chr_num,en2->mate_chr_num);
    fprintf(stderr,"Mate pos : %d\t\t|\tMate pos : %d\n",en1->mate_pos,en2->mate_pos);
}



// Return value :
// 0 : no decision for marking duplicates
// 1 : l_en1 is marked as duplicate
// -1: l_en2 is marked as duplicate
// 2 : Exit due to inputs being duplicates or secondary alignments


int is_duplicate_ens(half_mt_entry * en1, half_mt_entry * en2, int rev, int64_t ref_pos){
   
    int64_t tmp1 = 0;
    int64_t tmp2 = 0;

    if(((en1->flags & 0x400) != 0) || ((en2->flags & 0x400) != 0) || ((en1->flags & 0x100) != 0) || ((en2->flags & 0x100) != 0) || ((en1->flags & 0x4) != 0) || ((en2->flags & 0x4) != 0)){
        // either entry is a duplicate or a secondary alignment or a unmapped read
        return 2;
    }


    if(((en1->flags & 0x02) != 0) && ((en2->flags & 0x02) != 0)){
        // Both reads are properly paired
        // Should have mates
        // Check mate as well
        
        half_mt_entry * en1_m = (half_mt_entry *) en1->mate;
        half_mt_entry * en2_m = (half_mt_entry *) en2->mate;
        if(sort_verbose >= 5){
            fprintf(stderr,"------------------------------------------\n");
            fprintf(stderr,"%s\n",__func__);
            print_mt_entry(en1,en2);
            fprintf(stderr,"===========================================\n");
            print_mt_entry(en1_m,en2_m);
        }

        if((en1->mate_diff == en2->mate_diff) && ((en1->flags & 0x10) == (en2->flags & 0x10)) && (en1->mate_chr_num = en2->mate_chr_num) && (en1->mate_pos == en2->mate_pos)){
            // Mate diff is the same
            // Check the mates position now
            if((en1_m->ref_pos == en2_m->ref_pos) && ((en1_m->flags & 0x10) == (en2_m->flags & 0x10))){
                // mates position match and they have the same orientation

                if(en1->avg_qual < en2->avg_qual && en1_m->avg_qual < en2_m->avg_qual){
                    // en1 is a duplicate
                    if(sort_verbose >= 5){
                        fprintf(stderr,"1 was marked duplicate along with mate due to low avg_qual\n");
                    }
                    en1_m->flags |= 0x400;
                    en1_m->flags |= 0x8000;
                    return 1;
                }
                else if(en2->avg_qual < en1->avg_qual && en2_m->avg_qual < en1_m->avg_qual){
                    if(sort_verbose >= 5){
                        fprintf(stderr,"2 was marked duplicate along with mate due to low avg_qual\n");
                    }
                    // en2 is the duplicate
                    en2_m->flags |= 0x400;
                    en2_m->flags |= 0x8000;
                    return -1;
                }
                else{
                    // Check the fileptrs, the larger one is marked as the duplicate
                    memcpy(&tmp1,en1->fileptr,5*sizeof(uint8_t));
                    memcpy(&tmp2,en2->fileptr,5*sizeof(uint8_t));

                    if(tmp1 < tmp2){
                        if(sort_verbose >= 5){
                            fprintf(stderr,"2 was marked duplicate along with mate due to higher fileptr\n");
                        }
                        en2_m->flags |= 0x400;
                        en2_m->flags |= 0x8000;
                        return -1;
                    }
                    else {
                        if(sort_verbose >= 5){
                            fprintf(stderr,"1 was marked duplicate along with mate due to higher fileptr\n");
                        }
                        en1_m->flags |= 0x400;
                        en1_m->flags |= 0x8000;
                        return -1;
                    }
                }
            }
            else{
                return 0;
            }
        }
        else{
            // Mate diff are different, report no decision
            return 0;
        }

    }
    else {
        if(sort_verbose >= 5){
            fprintf(stderr,"------------------------------------------\n");
            fprintf(stderr,"%s\n",__func__);
            print_mt_entry(en1,en2);
        }
        // en1 is not properly paired or en2 is properly paired
        // Check orientation of the read
        if(((en1->flags & 0x10) == (en2->flags & 0x10))){

            //check if mate is mapped and compare mate positions only if mate is mapped
            
            if((((en1->flags & 0x8) == 0) && ((en2->flags & 0x8) == 0) && ((en1->flags & 0x20) == (en2->flags & 0x20)) && (en1->mate_chr_num == en2->mate_chr_num) && (en1->mate_pos == en2->mate_pos)) || ((en1->flags & 0x8) != 0) || ((en2->flags & 0x8) != 0)){
                if(en1->avg_qual < en2->avg_qual){
                    // l_en1 is the duplicate
                    if(sort_verbose >= 5){
                        fprintf(stderr,"1 was marked duplicate due to lower avg_qual\n");
                    }

                    return 1;
                }
                else if(en1->avg_qual == en2->avg_qual){
                    // Avg quality scores are the same

                    // Check relative positions now

                    if(rev == 1){
                        tmp1 = ref_pos - en1->correction;
                        tmp2 = ref_pos - en2->correction; 
                    }
                    else{
                        tmp1 = ref_pos + en1->correction;
                        tmp2 = ref_pos + en2->correction; 
                    }

                    if(tmp1 < tmp2){
                        //en2 is a duplicate
                        if(sort_verbose >= 5){
                            fprintf(stderr,"2 was marked duplicate due to higher sorted_position\n");
                        }
                        return -1;
                    }
                    else if(tmp1 > tmp2){
                        if(sort_verbose >= 5){
                            fprintf(stderr,"1 was marked duplicate due to higher sorted_position\n");
                        }
                        return 1;
                    }
                    else{
                        memcpy(&tmp1,en1->fileptr,5*sizeof(uint8_t));
                        memcpy(&tmp2,en2->fileptr,5*sizeof(uint8_t));

                        if(tmp1 < tmp2){
                            if(sort_verbose >= 5){
                                fprintf(stderr,"2 was marked duplicate due to higher fileptr\n");
                            }
                            return -1;
                        }
                        else {
                            if(sort_verbose >= 5){
                                fprintf(stderr,"1 was marked duplicate due to higher fileptr\n");
                            }
                            return 1;
                        }
                    }
                }
                else{
                    if(sort_verbose >= 5){
                        fprintf(stderr,"2 was marked duplicate due to lower avg_qual\n");
                    }
                    // l_en2 is duplicate
                    return -1;
                }
            }
            else{
                return 0;
            }
        }
        else{
            return 0;
        }

    }

    return 0;
}

// Return value :
// 0 : no decision for marking duplicates
// 1 : l_en1 is marked as duplicate
// -1: l_en2 is marked as duplicate
// 2 : Exit due to inputs being duplicates or secondary alignments

int is_duplicate_seq_en(bseq1s_t * en1, half_mt_entry * en2, int rev, int64_t ref_pos){
   
    int64_t tmp1 = 0;
    int64_t tmp2 = 0;

    if(((en1->flags & 0x400) != 0) || ((en2->flags & 0x400) != 0) || ((en1->flags & 0x100) != 0) || ((en2->flags & 0x100) != 0)){
        // either entry is a duplicate or a secondary alignment
        return 2;
    }


    if(((en1->flags & 0x02) != 0) && ((en2->flags & 0x02) != 0)){
        // Both reads are properly paired
        // Should have mates
        // Check mate as well
        
        bseq1s_t * en1_m = (bseq1s_t *) en1->mate;
        half_mt_entry * en2_m = (half_mt_entry *) en2->mate;
        if(sort_verbose >= 5){
            fprintf(stderr,"------------------------------------------\n");
            fprintf(stderr,"%s\n",__func__);
            print_seq_mt_entry(en1,en2);
            fprintf(stderr,"===========================================\n");
            print_seq_mt_entry(en1_m,en2_m);
        }

        //if(en1->mate_diff == en2->mate_diff && ((en1->flags & 0x10) == (en2->flags & 0x10)) && ((en1->flags & 0x40) == (en2->flags & 0x40))){
        if((en1->mate_diff == en2->mate_diff) && ((en1->flags & 0x10) == (en2->flags & 0x10)) && (en1->mate_chr_num = en2->mate_chr_num) && (en1->mate_pos == en2->mate_pos)){
            // Mate diff is the same
            // Mate diff is the same
            // Check the mates position now
            if((en1_m->abs_pos == en2_m->ref_pos) && ((en1_m->flags & 0x10) == (en2_m->flags & 0x10))){
                // mates position match and they have the same orientation

                if(en1->avg_qual < en2->avg_qual && en1_m->avg_qual < en2_m->avg_qual){
                    // en1 is a duplicate
                    if(sort_verbose >= 5){
                        fprintf(stderr,"1 was marked duplicate along with mate due to low avg_qual\n");
                    }
                    en1->flags |= 0x400;
                    en1->flags |= 0x8000;
                    en1_m->flags |= 0x400;
                    en1_m->flags |= 0x8000;
                    return 1;
                }
                else if(en2->avg_qual < en1->avg_qual && en2_m->avg_qual < en1_m->avg_qual){
                    // en2 is the duplicate
                    if(sort_verbose >= 5){
                        fprintf(stderr,"2 was marked duplicate along with mate due to low avg_qual\n");
                    }
                    en2->flags |= 0x400;
                    en2->flags |= 0x8000;
                    en2_m->flags |= 0x400;
                    en2_m->flags |= 0x8000;
                    return -1;
                }
                else{
                    // Check the fileptrs, the larger one is marked as the duplicate
                    memcpy(&tmp1,en1->fileptr,5*sizeof(uint8_t));
                    memcpy(&tmp2,en2->fileptr,5*sizeof(uint8_t));

                    if(tmp1 < tmp2){
                        if(sort_verbose >= 5){
                            fprintf(stderr,"2 was marked duplicate along with mate due to high fileptr\n");
                        }
                        en2->flags |= 0x400;
                        en2->flags |= 0x8000;
                        en2_m->flags |= 0x400;
                        en2_m->flags |= 0x8000;
                        return -1;
                    }
                    else {
                        if(sort_verbose >= 5){
                            fprintf(stderr,"1 was marked duplicate along with mate due to high fileptr\n");
                        }
                        en1->flags |= 0x400;
                        en1->flags |= 0x8000;
                        en1_m->flags |= 0x400;
                        en1_m->flags |= 0x8000;
                        return 1;
                    }
                }
            }
            else{
                return 0;
            }
        }
        else{
            // Mate diff are different or orientation is different, report no decision
            return 0;
        }

    }
    else {
        if(sort_verbose >= 5){
            fprintf(stderr,"------------------------------------------\n");
            fprintf(stderr,"%s\n",__func__);
            print_seq_mt_entry(en1,en2);
        }

        // en1 is not properly paired or en2 is notproperly paired
        if(((en1->flags & 0x10) == (en2->flags & 0x10))){

            //check if mate is mapped and compare mate positions only if mate is unmapped
            
            if((((en1->flags & 0x8) == 0) && ((en2->flags & 0x8) == 0) && ((en1->flags & 0x20) == (en2->flags & 0x20)) && (en1->mate_chr_num == en2->mate_chr_num) && (en1->mate_pos == en2->mate_pos)) || ((en1->flags & 0x8) != 0) || ((en2->flags & 0x8) != 0)){


                // if both entries have mapped mates, then check the orientations of the mate
                // otherwise ignore mate orientation

                if(en1->avg_qual < en2->avg_qual){
                    // l_en1 is the duplicate
                    if(sort_verbose >= 5){
                        fprintf(stderr,"1 was marked duplicate due to low avg_qual\n");
                    }
                    en1->flags |= 0x400;
                    en1->flags |= 0x8000;

                    return 1;
                }
                else if(en1->avg_qual == en2->avg_qual){
                    // Avg quality scores are the same

                    // Check relative positions now

                    if(rev == 1){
                        tmp1 = ref_pos - en1->correction;
                        tmp2 = ref_pos - en2->correction; 
                    }
                    else{
                        tmp1 = ref_pos + en1->correction;
                        tmp2 = ref_pos + en2->correction; 
                    }

                    if(tmp1 < tmp2){
                        //en2 is a duplicate
                        if(sort_verbose >= 5){
                            fprintf(stderr,"2 was marked duplicate due to higher sorted_position\n");
                        }
                        en2->flags |= 0x400;
                        en2->flags |= 0x8000;
                        return -1;
                    }
                    else if(tmp1 > tmp2){
                        if(sort_verbose >= 5){
                            fprintf(stderr,"1 was marked duplicate due to higher sorted_position\n");
                        }
                        en1->flags |= 0x400;
                        en1->flags |= 0x8000;
                        return 1;
                    }
                    else{

                        // Sorted positions are also the same now compare position in pair
                        memcpy(&tmp1,en1->fileptr,5*sizeof(uint8_t));
                        memcpy(&tmp2,en2->fileptr,5*sizeof(uint8_t));

                        if(tmp1 < tmp2){
                            if(sort_verbose >= 5){
                                fprintf(stderr,"2 was marked duplicate due to higher fileptr\n");
                            }
                            en2->flags |= 0x400;
                            en2->flags |= 0x8000;
                            return -1;
                        }
                        else {
                            if(sort_verbose >= 5){
                                fprintf(stderr,"1 was marked duplicate due to higher fileptr\n");
                            }
                            en1->flags |= 0x400;
                            en1->flags |= 0x8000;
                            return 1;
                        }
                    }
                }
                else{
                    if(sort_verbose >= 5){
                        fprintf(stderr,"2 was marked duplicate due to lower avg_qual\n");
                    }
                    en2->flags |= 0x400;
                    en2->flags |= 0x8000;
                    // l_en2 is duplicate
                    return -1;
                }
            }
            else{
                return 0;
            }
        }
        else{
            return 0;
        }

    }

    return 0;
}


int md_comparator(const void **p, const void **q) 
{
    ot_entry* l = *((ot_entry **)p);
    ot_entry* r = *((ot_entry **)q); 
    
    if(l->ref_pos < r->ref_pos){
        return -1;
    }
    else if(l->ref_pos == r->ref_pos){

        // Put properly paired reads before improperly paired reads
        if(((l->ote.flags & 0x2) == 0) && ((r->ote.flags & 0x2) != 0)){
            // l is not properly paired and r is,
            // l goes before r
            return -1; 
        }
        else if(((l->ote.flags & 0x2) != 0) && ((r->ote.flags & 0x2) == 0)){
            // l is properly paired and r is not,
            // r goes before l
            return 1;
        }
        else {
            if(l->ote.avg_qual > r->ote.avg_qual){
                return -1;
            }
            else if(l->ote.avg_qual < r->ote.avg_qual){
                return 1;
            }
            else{
                return 0;
            }
        }
    }
    else{
        return 1;
    }
}




void usage(){
    printf("./bwa sort -I <input sam> -O <output sam> -v <verbose level> -t <no. of threads> -pe \n");
    printf("-pe Run for pair ended reads\n");
}


int get_sequence(char * line, size_t size, int64_t * chr_start_array, int * chr_index){
   
    if(line[1] != 'S'){
        return -1;
    }


    chr_start_array[*chr_index] = 0;
    

    int i = 0;
    for(i=0;i<size;i++){
        if(line[i] == 'L'){
            break;
        }
    }

    chr_start_array[*chr_index] = (int64_t)strtoll(line + i + 3, NULL,10);
    *chr_index += 1;
    return 1;
}


void get_record(char * line, size_t line_size, int *flags, int *chr_num, int64_t *pos, uint8_t *fcorr, uint8_t *rcorr,int *avg_qual, int *mate_diff, int *mate_chr_num, int *mate_pos){
    int i = 0;
    int field_num = 1;

    int8_t corr = 0;
    int fcorr_done = 0;
    int stop = 0;


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
            // CIGAR Field
            if((int)line[i] >= (int)'0' && (int)line[i] <= (int)'9'){
                corr = corr * 10 + (int)line[i] - 48;
            }
            else if(line[i] != 'M'){
                if(fcorr_done == 0){
                    *fcorr = *fcorr + corr;
                }
                else{
                    if(line[i] != 'I'){
                        *rcorr = *rcorr + corr;
                    }
                }
                corr = 0;
            }
            else if(line[i] == 'M'){
                //stop = 1;
                *rcorr = *rcorr + corr; 
                fcorr_done = 1;
                corr = 0;
            }
        }
        else if(field_num == 7 && stop == 0){
            if(line[i] == 'C' || line[i] == 'c' || line[i] == 'H' || line[i] == 'h' || line[i] == 'R' || line[i] == 'r'){
                continue;    
            }
            else if((int)line[i] >= (int)'0' && (int)line[i] <= (int)'9'){
                *mate_chr_num = (*chr_num * 10) + (int)line[i] - 48;
            }
            else if(line[i] == 'X' || line[i] == 'x'){
                *mate_chr_num = 23;
            }
            else if(line[i] == 'Y' || line[i] == 'y'){
                *mate_chr_num = 24;
            }
            else if(line[i] == '='){
                *mate_chr_num = *chr_num;
            }
            else{
                *mate_chr_num = 0;
                stop = 1;
            }
        }
        else if(field_num == 8 && stop == 0){
            if((int)line[i] >= (int)'0' && (int)line[i] <= (int)'9'){
                *mate_pos = (*mate_pos * 10) + (int)line[i] - 48;
            }
            else{
                *mate_pos = -1;
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
            if(((int)line[i] - 33) >= 15){
                *avg_qual = *avg_qual + (int)line[i] - 33;
            }
        }

    }
    
    return;
}



bseq1s_t * get_sam_record(char * line, size_t line_size,uint64_t fileptr){
    
    int flags = 0;          // Field 2
    int chr_num = 0;            // Field 3
    int64_t pos = 0;            // Field 4
    uint8_t fcorr = 0;         // From CIGAR Field 6
    uint8_t rcorr = 0;         // From CIGAR Field 6
    int avg_qual = 0;           // Field 11
    int mate_diff = 0;

    int mate_pos = 0;
    int mate_chr_num = 0;

    get_record(line, line_size, &flags, &chr_num, &pos, &fcorr, &rcorr,&avg_qual, &mate_diff, &mate_chr_num, &mate_pos);


    /*if(sort_verbose >= 3){
        fprintf(stderr,"%s",line);
        fprintf(stderr,"Fileptr ; %ld\n",fileptr);
        fprintf(stderr,"-----------------------------------------\n");
    }*/

    if(sort_verbose >= 10){
        fprintf(stderr,"[Generated] Flags : %d, chr_num : %d, pos : %ld, fcorr : %d, rcorr : %d, avg_qual : %d, mate_diff : %d\n",(flags & 0xFFF),chr_num,pos,fcorr,rcorr,avg_qual,mate_diff); 
    }

    // TODO create the table

    bseq1s_t * seq;

    if(pos == -1 || chr_num == -1){
        seq = (bseq1s_t *)malloc(sizeof(bseq1s_t));
        memcpy(seq->fileptr,&fileptr,5*sizeof(uint8_t));
        seq->sam_size = (uint16_t)line_size;
        seq->chr_num = chr_num;
        seq->abs_pos = pos;
        return seq;
    }

    seq = (bseq1s_t *)malloc(sizeof(bseq1s_t));
    memcpy(seq->fileptr,&fileptr,5*sizeof(uint8_t));
    seq->sam_size = (uint16_t)line_size;
    seq->chr_num = chr_num;
    seq->is_rev = (flags & 0x10) != 0 ? 1 : 0;
    seq->flags = flags;
    seq->flags |= 0x2000;
    
    if(seq->is_rev){
        // Use rcorr, add to pos
        seq->abs_pos = pos + rcorr;
        seq->correction = rcorr;
    }
    else{
        // Use fcorr. subtract from pos
        seq->abs_pos = pos - fcorr;
        seq->correction = fcorr;
    }
    seq->mate_diff = mate_diff;
    seq->avg_qual = avg_qual;
    seq->last_entry = 0;


    if(sort_verbose >= 10){
        fprintf(stderr,"[Inserted] Flags : %d, chr_num : %d, pos : %ld, correction : %d, avg_qual : %d, mate_diff : %d\n",(seq->flags & 0xFFF),seq->chr_num,seq->abs_pos,seq->correction,seq->avg_qual,seq->mate_diff); 
    }
   
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
    if(sld->fot_length > 0){
        qsort(sld->fot,sld->fot_length,sizeof(ot_entry *),md_comparator);
    }

    if(sld->rot_length > 0){
        qsort(sld->rot,sld->rot_length,sizeof(ot_entry *),md_comparator);
    }

    return;
}






int count_duplicates_from_table(sort_struct_t * sld){
    int64_t i = 0;
    int count = 0;

    for(i=0;i<sld->fot_length;i++){
        if((sld->fot[i]->ote.flags & 0x400) != 0){
            count++;
        }
    }
    for(i=0;i<sld->rot_length;i++){
        if((sld->rot[i]->ote.flags & 0x400) != 0){
            count++;
        }
    }
    return count;

}


void ini_dup_list(dup_list * in, int size){
    in->n = 0;
    in->size = size;
    in->inc = size;
    in->dups = (int64_t *) malloc(size * sizeof(int64_t));
    memset(in->dups,0,size * sizeof(int64_t));
}

void reset_dup_list(dup_list * in){
    memset(in->dups,0,(in->size * sizeof(int64_t)));
    in->n = 0;
}

void add_dup_list(dup_list * in, int64_t num){
    if(in->n == in->size){
        in->dups = realloc(in->dups,(in->size + in->inc) * sizeof(int64_t));
        in->size = in->size + in->inc;
    }
    in->dups[in->n] = num;
    in->n = in->n + 1;
}

void free_dup_list(dup_list * in){
    free(in->dups);
}



void mark_duplicates_in_ot(ot_entry ** ot, int64_t len, int rev, half_mt_entry ** mt){
    int64_t i,j = 0;
    int k = 0;
    ot_entry * e1 = 0;
    ot_entry * e2 = 0;
    int dup = 0;

    dup_list l;

    ini_dup_list(&l,100);
    int mode = 0;
    int64_t prev_ref = -1;



    for(i=0;i<len-1;i++){
        e1 = ot[i];
        half_mt_entry * e1_ote = &e1->ote;
        half_mt_entry * e1_ote_m = (half_mt_entry *)e1->ote.mate;

        if(prev_ref != e1->ref_pos){
            // TODO: Slow, try to change once functionality is correct
            if(mt[e1->ref_pos] != NULL){
                if((mt[e1->ref_pos]->flags & 0x2) == 0){
                    e1_ote = mt[e1->ref_pos];
                    e1_ote_m = (half_mt_entry *)mt[e1->ref_pos]->mate;
                    i = i - 1;
                }
            }
            prev_ref = e1->ref_pos;
        }




        if(((e1_ote->flags & 0x400) != 0) || ((e1_ote->flags & 0x100) != 0) || ((e1_ote->flags & 0x4) != 0)){
            // Ignore entry if it already marked as a duplicate, is a secondary alignment and is an unmapped entry
            continue;
        }
        if((e1_ote->flags & 0x2) == 0){
            mode = 1;
        }
        else{
            mode = 0;
        }
        for(j=i+1;j<len;j++){
            e2 = ot[j];
            half_mt_entry * e2_ote = &e2->ote;
            half_mt_entry * e2_ote_m = (half_mt_entry *)e2->ote.mate;
            if((e1->ref_pos != e2->ref_pos)){
                // Crossed over into another ref_pos territory without seeing entries in between that were not duplicates
                break;
            }
            else{
                
                // check if the sorted table has unpaired entries after a paired entry
                if(mode == 0){
                    assert((e2_ote->flags & 0x2) != 0);
                }

                // Check for duplicates here

                if(e1_ote_m == e2_ote){
                    assert(e2_ote_m == e1_ote);
                    dup = 0;
                }
                else{
                    // Test for Mark duplicate only if the two reads are not mates of each other
                    dup = is_duplicate_ens(e1_ote, e2_ote, rev, e1->ref_pos);
                }

                if(dup == 1){
                    // Restart run from i + 1;
                    
                    // e1 was duplicate, if mode was 1 then we were comparing unpaired read with all, ignore all duplicates stored in the list
                    mode = 0;
                    e1_ote->flags |= 0x400;
                    e1_ote->flags |= 0x8000;
                    break;
                }
                else if(dup == -1){
                    // e2 was duplicate
                    if(mode == 1){
                        // if mode was 1 then we were comparing with an unpaired read,
                        // dont mark e2 as duplicate but store it in the dup_list
                        add_dup_list(&l, j);
                    }
                    else{
                        // mode was 0, safe to add flags for e2
                        e2_ote->flags |= 0x400;
                        e2_ote->flags |= 0x8000;
                    }
                }
            }
            
        }

        // check mode
        if(mode == 1){
            for(k = 0;k<l.n;k++){
                int64_t l1 = l.dups[k];
                half_mt_entry * tmp = &ot[l1]->ote;
                tmp->flags |= 0x400;
                tmp->flags |= 0x8000;
            }
        }

        mode = 0;
        reset_dup_list(&l);

    }

    free_dup_list(&l);
}

void mark_duplicates(sort_struct_t * sld){
    if(sld->fot_length > 0){
        mark_duplicates_in_ot(sld->fot, sld->fot_length, 0, sld->fmt);
    }
    if(sld->rot_length > 0){
        mark_duplicates_in_ot(sld->rot, sld->rot_length, 1, sld->rmt);
    }
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
            if(seq->mate != NULL){
                free(seq->mate);
            }
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
    
    // Sorting the overflow table

    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&proc_start);
    mark_duplicates(sld);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&proc_end);
    s->md_time = (double)(proc_end.tv_sec - proc_start.tv_sec) + ((double)(proc_end.tv_nsec - proc_start.tv_nsec)/(double)(1000000000));
    // Sorting end
    s->num_duplicates = count_duplicates_from_table(sld);

    s->sld = sld;

    free(queue_out);
    pthread_exit(0);
}





// If return value is 0, both are equivalent
// If return value is 1, len1 should go before len0
// If return value is -1, len0 should go before len1


int sort_comparator(const void **p, const void **q) 
{
    sort_list_t * l = *((sort_list_t **)p);
    sort_list_t * r = *((sort_list_t **)q); 
    
    int dec = 0;




    if((l->sort_ref_pos) < (r->sort_ref_pos)){
        return -1;
    }
    else if((r->sort_ref_pos) < (l->sort_ref_pos)){
        return 1;
    }
    else{
        // Both references are equal, compare read names in lexicographic order

        // Check for orientation, forward appears before reverse
        if(((l->ote->flags & 0x10) == 0) && ((r->ote->flags & 0x10) != 0)){
            return -1;
        }
        else if(((r->ote->flags & 0x10) == 0) && ((l->ote->flags & 0x10) != 0)){
            return 1;
        }
        else{

            // Orientations are the same, check for lexicographic ordering

            dec = lexicographic_comparator(l->sam,l->name_len, r->sam,r->name_len);
            if(dec != 0){
                return dec;
            }
            else{
                // First in pair goes first
                if(((l->ote->flags & 0x40) == 0) && ((r->ote->flags & 0x40) != 0)){
                    // r is first in pair
                    return 1;
                }
                if(((l->ote->flags & 0x40) != 0) && ((r->ote->flags & 0x40) == 0)){
                    return -1;
                }
                if(((l->ote->flags & 0x100) == 0) && ((r->ote->flags & 0x100) != 0)){
                    // r is a secondary alignment
                    // l goes before r
                    return -1;
                }
                else if(((r->ote->flags & 0x100) == 0) && ((l->ote->flags & 0x100) != 0)){
                    // l is a secondary alignment
                    // r goes before l
                    return 1;
                }
                else if(((l->ote->flags & 0x4) == 0) && ((r->ote->flags & 0x4) != 0)){
                    // r is an unmapped alignment
                    // l goes before r
                    return -1;
                }
                else if(((r->ote->flags & 0x4) == 0) && ((l->ote->flags & 0x4) != 0)){
                    // l is an unmapped alignment
                    // r goes before l
                    return 1;
                }
                else{
                    return 0;
                }
            }
        }


    }
}

int umt_comparator(const void **p, const void **q) 
{
    unmapped_entry * l = *((unmapped_entry **)p);
    unmapped_entry * r = *((unmapped_entry **)q); 
    
    int dec = 0;
    dec = lexicographic_comparator(l->sam,l->name_len, r->sam,r->name_len);
    return dec;

}

void add_entry_to_sort_list(half_mt_entry * ote, int64_t sort_ref_pos, int64_t md_ref_pos, sort_list * l_in){
    if(l_in->n != 0 && (l_in->n % 100) == 0){
        l_in->list = (sort_list_t**)realloc(l_in->list,(l_in->n + 100)*sizeof(sort_list_t*));
    }

    l_in->list[l_in->n] = (sort_list_t*) malloc(sizeof(sort_list_t));
    l_in->list[l_in->n]->sort_ref_pos = sort_ref_pos;
    l_in->list[l_in->n]->md_ref_pos = md_ref_pos;
    l_in->list[l_in->n]->ote = ote;
    l_in->n = l_in->n + 1;
}

void print_sort_list(sort_list * l){
    int i = 0;
    int64_t fp = 0;
    for(i = 0;i<l->n;i++){
        memcpy(&fp,l->list[i]->ote->fileptr,5*sizeof(uint8_t));
        fprintf(stderr,"fp : %ld\n",fp);
    }
    if(l->n != 0){
        fprintf(stderr,"===================================\n");
    }
}


void get_sam_umt(FILE * in_sam, unmapped_entry * en){
   int64_t fp = 0;
   memcpy(&fp,en->fileptr,5*sizeof(uint8_t));
   fseek(in_sam,fp,SEEK_SET);
   en->sam = (char *) malloc((en->sam_size + 1) * sizeof(char)); 

   memset(en->sam,0,en->sam_size+1);
   fread(en->sam,sizeof(char),en->sam_size,in_sam);
   //sam_str[en->sam_size] = '\0';

   // Print name

   int i = 0;
   for(i=0;i<en->sam_size;i++){
       if(en->sam[i] == '\t'){
           en->name_len = i;
           break;
       }
   }

}






void get_sam(FILE * in_sam, half_mt_entry * en, char ** out_sam_str, int *name_len){
    int64_t fp = 0;
    memcpy(&fp,en->fileptr,5*sizeof(uint8_t));
    fseek(in_sam,fp,SEEK_SET);
    char * sam_str = (char *) malloc((en->sam_size + 1) * sizeof(char)); 

    memset(sam_str,0,en->sam_size+1);
    fread(sam_str,sizeof(char),en->sam_size,in_sam);
    sam_str[en->sam_size] = '\0';

    // Print name

    int i = 0;
    for(i=0;i<en->sam_size;i++){
        if(sam_str[i] == '\t'){
            *name_len = i;
            break;
        }
    }

    *out_sam_str = sam_str;
    en->flags |= 0x4000;
}


void update_sam_flags(char ** sam, half_mt_entry * en){
    char * old_sam = *sam;
    char * new_sam = (char *) malloc((en->sam_size + 4) * sizeof(char));

    //char * new_sam = *sam;
    memset(new_sam,0,en->sam_size+4);

    int i = 0;
    int flag_detected = 0;
    int out_sam_index = 0;

    for(i=0;i<en->sam_size;i++){
        if(flag_detected == 0){
            new_sam[i] = old_sam[i];
            //sprintf(out_sam_str,"%c",sam_str[i]);
        }
        if(old_sam[i] == '\t'){
            if(flag_detected == 0){
                flag_detected = 1;
                out_sam_index = i + 1;
            }
            else{
                break;
            }
        }
    }
    sprintf(new_sam + out_sam_index,"%d%s",(en->flags & 0xFFF),old_sam+i);

    *sam = new_sam;
    free(old_sam);

    return;
}

void update_sam_records_for_list(sort_list * l){
    int i = 0;
    for(i=0;i<l->n;i++){
        update_sam_flags(&l->list[i]->sam, l->list[i]->ote);
    }
}



void get_sam_records_for_list(FILE * in_sam,sort_list * l){
    int i = 0;
    for(i=0;i<l->n;i++){
        get_sam(in_sam,l->list[i]->ote,&l->list[i]->sam, &l->list[i]->name_len);
    }
    return;
}


void get_valid_entries_list(FILE * in_sam,int64_t ref_pos,sort_struct_t * sld,int64_t head_fot, int64_t head_rot,sort_list * l, half_mt_entry * in_en){
    //sort_list * l = (sort_list *) malloc(sizeof(sort_list));
    l->list = (sort_list_t **) malloc(100 * sizeof(sort_list_t *));
    l->n = 0;
    l->ot_entries_start = 0;
   

    int max_correction = 0;
    int mandatory_fcorr = 101;
    int mandatory_rcorr = 256;

    if(in_en == NULL){
        max_correction = 0;
    }
    else{
        max_correction = in_en->correction;
        if((in_en->flags & 0x4000) != 0){
            return;
        }
        add_entry_to_sort_list(in_en, ref_pos + in_en->correction, ref_pos,l);
        in_en->flags |= 0x4000;
    }

    int64_t i = 0;
    int64_t i1 = 0;
    // Get entries from forward main table
    // First get all entries that might be duplicates of this entry





    for(i = ref_pos;i<=ref_pos + mandatory_fcorr;i++){


        // Check if there is any entry in the ot table that can go before the mt entry

        if(head_fot < sld->fot_length){
            while(sld->fot[head_fot]->ref_pos <= i){
                i1 = sld->fot[head_fot]->ref_pos + sld->fot[head_fot]->ote.correction;
                if(((sld->fot[head_fot]->ote.flags & 0x4000) == 0) && (i1 <= (ref_pos + max_correction))){
                    add_entry_to_sort_list(&sld->fot[head_fot]->ote,i1,sld->fot[head_fot]->ref_pos,l);
                }
                head_fot = head_fot + 1;
                if(head_fot >= sld->fot_length){
                    break;    
                }
            }
        }

        if((sld->fmt[i] != 0)){
            if((sld->fmt[i]->flags & 0x4000) == 0){
            // Entry hasnt been seen yet
                if((sld->fmt[i]->correction + i) <= (ref_pos + max_correction)){
                    // Add entry to list
                    add_entry_to_sort_list(sld->fmt[i],i+sld->fmt[i]->correction,i,l);
                }
            }
        }
    }

    // Get entries from reverse main table
    for(i = ref_pos;i<=ref_pos + mandatory_rcorr;i++){

        if(head_rot < sld->rot_length){
            while(sld->rot[head_rot]->ref_pos <= i){
                i1 = sld->rot[head_rot]->ref_pos - sld->rot[head_rot]->ote.correction;
                if(((sld->rot[head_rot]->ote.flags & 0x4000) == 0) && (i1 <= (ref_pos + max_correction))){
                    add_entry_to_sort_list(&sld->rot[head_rot]->ote,i1,sld->rot[head_rot]->ref_pos,l);
                }
                head_rot = head_rot + 1;
                if(head_rot >= sld->rot_length){
                    break;    
                }
            }
        }

        if((sld->rmt[i] != 0)){
            if(((sld->rmt[i]->flags & 0x4000) == 0)){
            // Entry hasnt been seen yet
                if((i - sld->rmt[i]->correction) <= (ref_pos + max_correction)){
                    // Add entry to list
                    add_entry_to_sort_list(sld->rmt[i],i-sld->rmt[i]->correction,i,l);
                }
            }
        }
    }

    
    if(l->n != 0){
        // Get all sam records before sorting to compare names as well
        get_sam_records_for_list(in_sam,l);
        qsort(l->list,l->n,sizeof(sort_list_t *),sort_comparator);
        update_sam_records_for_list(l);

    }

    //*l_in = l;
    return;

}





/*void print_sam(FILE * in_sam, FILE * out_sam, sort_list_t * l_en){
   int64_t fp = 0;
   memcpy(&fp,en->fileptr,5*sizeof(uint8_t));
   fseek(in_sam,fp,SEEK_SET);
   char * sam_str = (char *) malloc((en->sam_size + 1) * sizeof(char));
   memset(sam_str,0,en->sam_size+1);
   fread(sam_str,sizeof(char),en->sam_size,in_sam);
   sam_str[en->sam_size] = '\0';

   // Print name

   int i = 0;
   int flag_detected = 0;
   for(i=0;i<en->sam_size;i++){
       if(flag_detected == 0){
           fprintf(out_sam,"%c",sam_str[i]);
       }
       if(sam_str[i] == '\t'){
           if(flag_detected == 0){
               flag_detected = 1;
           }
           else{
               break;
           }
       }
   }

   fprintf(out_sam,"%d",en->flags & 0xFFF);
   fprintf(out_sam,"%s",sam_str + i);
   free(sam_str);

   
   en->flags |= 0x4000;

}*/

void free_sort_list(sort_list * l){
    int i = 0;
    if(l->list){
        for(i=0;i<l->n;i++){
            if(l->list[i]){
                l->list[i]->ote = NULL;
                if(l->list[i]->sam){
                    free(l->list[i]->sam);
                }
                free(l->list[i]);
            }
        }
        free(l->list);
    }
}


void update_fot_head(int64_t * head_fot, sort_struct_t * sld){
    
    if(sld->fot_length == 0){
        return;
    }

    if(*head_fot >= sld->fot_length){
        *head_fot = sld->fot_length;
        return;
    }

    //Increment head till it points to the first not seen entry

    while((sld->fot[*head_fot]->ote.flags & 0x4000) != 0){
        *head_fot = *head_fot + 1;
        if(*head_fot >= sld->fot_length){
            *head_fot = sld->fot_length;
            return;
        }
    }

}


void update_rot_head(int64_t * head_rot, sort_struct_t * sld, int64_t ref_pos){
    
    if(sld->rot_length == 0){
        return;
    }

    if(*head_rot >= sld->rot_length){
        *head_rot = sld->rot_length;
        return;
    }


    if((sld->rot[*head_rot]->ref_pos) < ref_pos){
        // Do first try only if required ref pos is greater than current rot head ref pos
        while((sld->rot[*head_rot]->ref_pos) < ref_pos){
            *head_rot = *head_rot + 1;
            if(*head_rot >= sld->rot_length){
                *head_rot = sld->rot_length;
                return;
            }
        }
    }

    //Increment head till it points to the first not seen entry

    while((sld->rot[*head_rot]->ote.flags & 0x4000) != 0){
        *head_rot = *head_rot + 1;
        if(*head_rot >= sld->rot_length){
            *head_rot = sld->rot_length;
            return;
        }
    }

}



void generate_sorted_sam(FILE * in_sam, FILE * out_sam, sort_slave_t * sl){
    // Resort the OT table
    sort_struct_t * sld = sl->sld;
    int64_t mt_length = sl->length;

    int64_t i = 0;
    int64_t head_fot = 0;
    int64_t head_rot = 0;
    int64_t head_fot_ref_pos = 0;

    int j = 0;

    sort_list * l = (sort_list *) malloc(sizeof(sort_list));


    // Only go through the FMT for now. Hope it can catch all reverse main table entries as well
    for(i=0;i<mt_length;i++){
        
        //Check if any entry in the mt at i is valid
       

        //update_rot_head(&head_rot,sld,i);


        if(head_fot < sld->fot_length){
            head_fot_ref_pos = sld->fot[head_fot]->ref_pos;
            while(head_fot_ref_pos < i){
                update_rot_head(&head_rot,sld,head_fot_ref_pos);
                get_valid_entries_list(in_sam,head_fot_ref_pos,sld,head_fot,head_rot,l, &sld->fot[head_fot]->ote);

                for(j=0;j<l->n;j++){
                    // Print sam record
                    fprintf(out_sam,"%s",l->list[j]->sam);
                    //print_sam(in_sam,out_sam,l->list[j]->ote);
                }
                free_sort_list(l);
                l->n = 0;

                head_fot++;
                if(head_fot >= sld->fot_length){
                    head_fot = sld->fot_length;
                    break;
                }

                head_fot_ref_pos = sld->fot[head_fot]->ref_pos;
            }
        }

        update_fot_head(&head_fot,sld);


        if(sld->fmt[i] != 0){
            update_rot_head(&head_rot,sld,i);
            //fmt entry is valid
            get_valid_entries_list(in_sam,i,sld,head_fot,head_rot,l, sld->fmt[i]);
            for(j=0;j<l->n;j++){
                // Print sam record
                fprintf(out_sam,"%s",l->list[j]->sam);
                //print_sam(in_sam,out_sam,l->list[j]->ote);
            }
            free_sort_list(l);
            l->n = 0;
        }
        else{
            if(sld->rmt[i] != 0){
                update_rot_head(&head_rot,sld,i);
                get_valid_entries_list(in_sam,i,sld,head_fot,head_rot,l, NULL);
                for(j=0;j<l->n;j++){
                    // Print sam record
                    fprintf(out_sam,"%s",l->list[j]->sam);
                    //print_sam(in_sam,out_sam,l->list[j]->ote);
                }
                free_sort_list(l);
                l->n = 0;
            }
        }

    }
    free(l);
    
}





void add_entry_to_umt(bseq1s_t * seq, unmapped_entry_v * umt){

    if(umt->n != 0 && (umt->n % INIT_OVERFLOW_TABLE_SIZE == 0)){
        umt->list = realloc(umt->list,(umt->n + INIT_OVERFLOW_TABLE_SIZE) * sizeof(unmapped_entry));
    }

    int index = umt->n;
    umt->list[index] = malloc(sizeof(unmapped_entry));
    memcpy(umt->list[index]->fileptr,seq->fileptr,5*sizeof(uint8_t));
    umt->list[index]->sam_size = seq->sam_size;
    umt->list[index]->sam = NULL;
    umt->list[index]->name_len = 0;
    umt->n = umt->n + 1;
    return;
}

void count_valid_entries(sort_struct_t * sld){
    int64_t i = 0;
    int64_t count_valid_entries = 0;
    fprintf(stderr,"mt_length allocated : %ld\n",sld->mt_length);
    for(i=0;i<sld->mt_length;i++){
        if(sld->fmt[i] != 0){
            count_valid_entries++;
        }
    }
    fprintf(stderr,"Total valid entries in fmt: %ld\n",count_valid_entries);
    count_valid_entries = 0;
    for(i=0;i<sld->mt_length;i++){
        if(sld->rmt[i] != 0){
            count_valid_entries++;
        }
    }
    fprintf(stderr,"Total valid entries in rmt: %ld\n",count_valid_entries);
    fprintf(stderr,"Total entries in ots : %ld\n",sld->fot_length + sld->rot_length);

}



void sort_MT(char * in_sam_filename,char * out_sam_filename,int num_threads, int pair_ended){
    int64_t * chr_len;
    char * line = NULL;
    size_t line_size = 0;
    size_t ret = 0;
    uint32_t read_id = 0;

    int num_chromosomes = 0;
    int i = 0;

    FILE * in_sam = fopen(in_sam_filename,"r");
    if(in_sam == NULL){
        fprintf(stderr,"File %s not found\n",in_sam_filename);
        return ;
    }

    FILE * out_sam = NULL;



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
    bseq1s_t * m;

    unmapped_entry_v umt;
    umt.n = 0;
    umt.list = (unmapped_entry **) malloc(INIT_OVERFLOW_TABLE_SIZE * sizeof(unmapped_entry*));

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

    uint64_t fileptr = 0;

    queue * sending_queue = NULL;
    if(sort_verbose >= 3){
        clock_gettime(CLOCK_REALTIME,&read_start);
    }
    while(!feof(in_sam)){
        fileptr = (uint64_t)ftell(in_sam);
        ret = getline(&line,&line_size,in_sam);
        if(ret != -1 && line[0] != '@'){
            s = get_sam_record(line, strlen(line), fileptr);
            read_id++;
            if(s->chr_num != -1 && s->abs_pos != -1){
            //if(s != NULL){

                // Entry has valid sam record
                // Search for mat if the entry is properly paired and not a secondary alignment
                
                if(((s->flags & 0x1) != 0) && ((s->flags & 0x100) == 0)){
                    // s is not a secondary alignment and is a pair ended read
                    while(1){
                        fileptr = (uint64_t)ftell(in_sam);
                        ret = getline(&line,&line_size,in_sam);
                        if(ret != -1 && line[0] != '@'){
                            m = get_sam_record(line,strlen(line),fileptr);
                            read_id++;
                        }

                        // The mate should go to the same chromosome and also be properly paired and not a secondary alignment;

                        if(((m->flags & 0x100) == 0)){

                            s->mate_chr_num = m->chr_num;
                            s->mate_pos = m->abs_pos;
                            m->mate_chr_num = s->chr_num;
                            m->mate_pos = s->abs_pos;
                            if(m->chr_num != -1 && m->abs_pos != -1){
                                if(m->chr_num != s->chr_num){
                                    m->abs_pos += chr_len[m->chr_num - 1];
                                    m->mate = NULL; 
                                    tid = chr_thread_id[m->chr_num - 1];
                                    sending_queue = qs[tid];
                                    addElement(sending_queue, (void *)m);
                                    
                                    s->mate = NULL;
                                    
                                }
                                else{
                                    s->mate = (void *)m;
                                    m->mate = (void *)s;
                                }
                            }
                            else{
                                add_entry_to_umt(m,&umt);
                                free(m);
                            }

                            break; 
                        }
                        else{
                            // mate is secodary alignment send it independantly
                            if(m->chr_num != -1 && m->abs_pos != -1){
                                m->abs_pos += chr_len[m->chr_num - 1];
                                m->mate = NULL; 
                                tid = chr_thread_id[m->chr_num - 1];
                                sending_queue = qs[tid];
                                addElement(sending_queue, (void *)m);
                            }
                            else{
                                add_entry_to_umt(m,&umt);
                                free(m);
                            }
                        }
                    }
                }
                else{
                    s->mate = NULL;
                }

                s->abs_pos += chr_len[s->chr_num - 1];
                if(s->mate != NULL){
                    //m = (bseq1s_t * ) s->mate;
                    m->abs_pos += chr_len[s->chr_num - 1];
                }
                tid = chr_thread_id[s->chr_num - 1];
                sending_queue = qs[tid];
                addElement(sending_queue, (void *)s);
            }
            else{
                // Add entry to unmapped table
                add_entry_to_umt(s,&umt);
                free(s);
            } 
            if(sort_verbose >= 3 && read_id % 10000000 == 0){
                clock_gettime(CLOCK_REALTIME,&read_end);
                read_time = (double)(read_end.tv_sec - read_start.tv_sec) + ((double)(read_end.tv_nsec - read_start.tv_nsec)/(double)(1000000000));
                fprintf(stderr,"Read %d reads in %f secs\n",read_id,read_time);
            }
        }
    }

    fclose(in_sam);

    for(i=0;i<num_threads;i++){
        sending_queue = qs[i];
        s = (bseq1s_t *) malloc( sizeof(bseq1s_t));
        s->fileptr[0] = 0;
        s->fileptr[1] = 0;
        s->fileptr[2] = 0;
        s->fileptr[3] = 0;
        s->fileptr[4] = 0;
        s->sam_size = 0;
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
        fprintf(stderr,"Read %d reads in %f secs\n",read_id,read_time);
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
    for(i=0;i<num_threads;i++){
        total_duplicates += slaves[i].num_duplicates;
        if(sort_verbose >= 5){
            fprintf(stderr,"For tid : %d\n",i+1);
            fprintf(stderr,"\t%ld Reads processed in %f secs\n",slaves[i].num_reads,slaves[i].write_time);
            fprintf(stderr,"\t%d Duplicates detected\n",slaves[i].num_duplicates);

            fprintf(stderr,"\tSorting time : %f\n",slaves[i].sort_time);
            fprintf(stderr,"\tMD time : %f\n",slaves[i].md_time);
            fprintf(stderr,"\tMark duplicates time : %f\n",slaves[i].md_time);
        }
    }

    fprintf(stderr,"Total duplicates detected : %d\n",total_duplicates);


    // Print output sam file

    if(out_sam_filename == NULL){
        out_sam = stdout;
    }
    else{
        out_sam = fopen(out_sam_filename,"w");
    }
    in_sam = fopen(in_sam_filename,"r");

    for(i=0;i<num_threads;i++){
        generate_sorted_sam(in_sam, out_sam,&slaves[i]);
        //fprintf(stderr,"Statistics for entries for tid : %d\n",i);
        //count_valid_entries(slaves[i].sld);
        struct_delete(slaves[i].sld);
        free(slaves[i].sld);
        //fprintf(stderr,"=================================\n");
    }


    // Get all strings for unmapped entries

    /*for(i=0;i<umt.n;i++){
        get_sam_umt(in_sam, umt.list[i]);
    }

    if(umt.n > 0){
        qsort(umt.list,umt.n,sizeof(unmapped_entry *),umt_comparator);
    }*/

    for(i=0;i<umt.n;i++){
        //get_sam_umt(in_sam, umt.list[i]);
        //fprintf(out_sam,"%s",umt.list[i]->sam);
        free(umt.list[i]->sam);
        free(umt.list[i]);
    }

    // Free umt

    free(umt.list);


    free(slaves);
    free(pts);
    free(qs);

    free(lens);
    free(chr_thread_id);
    free(chr_len);
    fclose(in_sam);
    fclose(out_sam);
    return;
}

int main_memsort(int argc, char *argv[]){
    
    int rargc = 2;
    int i = 0;
    int output_file = 0; 
    int num_threads = 1;
    int pair_ended = 0;
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
        if(strcmp(argv[i],"-pe") == 0){
            pair_ended = 1;
            continue;
        }
    }


    if(output_file == 1){
        sort_MT(temp_unsorted_file_name,sorted_file_name,num_threads, pair_ended);
    }
    else{
        sort_MT(temp_unsorted_file_name,NULL,num_threads,pair_ended);
    }

    //fclose(temp_unsorted_file);
    return 0;

}
