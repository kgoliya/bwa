#include <pthread.h>
typedef struct {
    void **buf;
    long head, tail;
    int full, empty;
    int queue_size;
    pthread_mutex_t *mut;
    pthread_cond_t *notFull, *notEmpty;
} queue;


queue *queueInit (int queue_size);
void queueDelete (queue *q);
void queueAdd (queue *q, void* in);
void queueDel (queue *q, void ** out);
void getElement(queue * q, void ** out);
void addElement(queue * q, void * in);
