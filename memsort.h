#include "bwa.h"
#include "pc_queue.h"

typedef struct {
    uint8_t fileptr[5];
    uint16_t sam_size;
    uint8_t correction; 
    uint16_t avg_qual;

    // Bit 10 : Is duplicate or not(0x400)
    // Bit 1 : Is properly paired or not(0x2)
    // Bit 0 : Read paired (0x1)
    // Bit 4 : Is reverse (0x10)
    // Bit 15: Seen for mark duplicate processing(0x8000)
    // Bit 14: Seen for generating sam(0x4000)
    // Bit 13: valid(0x2000)
    uint16_t flags;
    uint16_t mate_diff;

    // mate entry's pointer
    void * mate;
} half_mt_entry;


/*typedef struct {
    half_mt_entry fe;
    half_mt_entry re;
} mt_entry;*/

typedef struct {
    half_mt_entry ote;
    int64_t ref_pos;    // Not corrected with CIGAR            
} ot_entry;

typedef struct {
    uint8_t fileptr[5];
    uint16_t sam_size;
    char * sam;
    int name_len;
} unmapped_entry;

typedef struct {
    int n;
    unmapped_entry ** list;
} unmapped_entry_v;

typedef struct {
    half_mt_entry ** fmt;
    half_mt_entry ** rmt;
    int64_t mt_length;

    ot_entry ** fot;
    int64_t fot_size;
    int64_t fot_length;

    ot_entry ** rot;
    int64_t rot_size;
    int64_t rot_length;

} sort_struct_t;

typedef struct {
    int thread_id;
    int64_t length;
    int num_duplicates;
    double write_time;
    double sort_time;
    double md_time;
    int64_t num_reads;

    sort_struct_t * sld;
    queue * q;
} sort_slave_t;

typedef struct {
    int64_t sort_ref_pos;
    int64_t md_ref_pos;
    char * sam;
    int name_len;
    half_mt_entry * ote;
} sort_list_t;


typedef struct {
    int n;
    int ot_entries_start;
    sort_list_t ** list;
} sort_list;


void sorting_init(int64_t l_pac);
void sorting_close();
void process_sam_record(bseq1_t *seq);
void get_sorted_sam();
