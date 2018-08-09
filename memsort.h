#include "bwa.h"


typedef struct {
    uint8_t file_ptr[5];
    int16_t sam_size;
    uint8_t correction; 
    uint8_t avg_qual;

    // Bit 0 : Is duplicate or not
    // Bit 1 : Is properly paired or not
    uint8_t flags;
    uint16_t mate_diff;
} sort_info;

typedef struct {
    uint8_t avg_qual;
    uint8_t correction;
    uint16_t mate_diff;
} md_info;


// Each entry should have sorter info and MD info
typedef struct {
    sort_info si;
    md_info mdi;
}half_mt_entry;

typedef struct {
    half_mt_entry fe;
    half_mt_entry re;
} mt_entry;


typedef struct {
    md_info mdi;
    int64_t ref_pos;    // Corrected with CIGAR
    uint8_t is_rev;
} md_ot_entry;


typedef struct {
    sort_info si;
    int64_t ref_pos;    // Not corrected with CIGAR            
    uint8_t is_rev;
} sort_ot_entry;

void sorting_init(int64_t l_pac);
void sorting_close();
void process_sam_record(bseq1_t *seq);
void get_sorted_sam();
