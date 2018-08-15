#include "bwa.h"
#include "pc_queue.h"

typedef struct {
    uint32_t read_id;
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

typedef struct {
    mt_entry * mt;
    int64_t mt_length;

    sort_ot_entry ** sort_ot;
    int sort_ot_size;
    int sort_ot_length;

    md_ot_entry ** md_ot;
    int md_ot_size;
    int md_ot_length;
} sort_struct_t;

typedef struct {
    int thread_id;
    int64_t length;
    double time;
    int64_t num_reads;
    queue * q;
} sort_slave_t;

void sorting_init(int64_t l_pac);
void sorting_close();
void process_sam_record(bseq1_t *seq);
void get_sorted_sam();
