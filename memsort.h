#include "bwa.h"

typedef struct {
    uint8_t file_ptr[5];
    uint8_t avg_qual;       // If avg quality score is zero, entry is invalid
    int16_t sam_size;
}half_mt_entry;

typedef struct {
    half_mt_entry fe;
    half_mt_entry re;
} mt_entry;

/*typedef struct {
    uint8_t file_ptr[5];
    uint8_t avg_qual;      
    int16_t sam_size;
} of_entry;*/

typedef struct {
    int64_t ref_pos;
    uint8_t file_ptr[5];
    int16_t sam_size;
} dpt_entry;

typedef struct {
    uint8_t file_ptr[5];
    int16_t sam_size;
} umt_entry;
// TODO: Printing order is different from order required by mark duplicates

void sorting_init(int64_t l_pac);
void sorting_close();
void process_sam_record(bseq1_t *seq);
void get_sorted_sam();
