#include "bwa.h"
#include "pc_queue.h"

typedef struct {
    uint32_t read_id;
    uint8_t correction; 
    uint8_t avg_qual;

    // Bit 10 : Is duplicate or not(0x400)
    // Bit 1 : Is properly paired or not(0x2)
    // Bit 0 : Read paired (0x1)
    // Bit 4 : Is reverse (0x10)
    // Bit 15: Seen for mark duplicate processing(0x8000)
    uint16_t flags;
    
    uint16_t mate_diff;
} half_mt_entry;


typedef struct {
    half_mt_entry fe;
    half_mt_entry re;
} mt_entry;

typedef struct {
    half_mt_entry ote;
    int64_t ref_pos;    // Not corrected with CIGAR            
} ot_entry;

typedef struct {
    mt_entry * mt;
    int64_t mt_length;

    ot_entry ** sort_ot;
    int64_t sort_ot_size;
    int64_t sort_ot_length;

} sort_struct_t;

typedef struct {
    int thread_id;
    int64_t length;
    int num_duplicates;

    double write_time;
    double sort_time;
    double md_time;
    int64_t num_reads;

    queue * q;
} sort_slave_t;

void sorting_init(int64_t l_pac);
void sorting_close();
void process_sam_record(bseq1_t *seq);
void get_sorted_sam();
