#include <pthread.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <poll.h>
#include "pc_queue.h"

queue *queueInit (int queue_size)
{
    queue *q;

    q = (queue *)malloc (sizeof (queue));
    if (q == NULL) return (NULL);

    q->empty = 1;
    q->full = 0;
    q->head = 0;
    q->tail = 0;
    q->queue_size = queue_size;
    q->mut = (pthread_mutex_t *) malloc (sizeof (pthread_mutex_t));
    pthread_mutex_init (q->mut, NULL);
    q->notFull = (pthread_cond_t *) malloc (sizeof (pthread_cond_t));
    pthread_cond_init (q->notFull, NULL);
    q->notEmpty = (pthread_cond_t *) malloc (sizeof (pthread_cond_t));
    pthread_cond_init (q->notEmpty, NULL);

    q->buf = (void**) malloc(queue_size * sizeof(void*));
    
    return (q);
}
void queueDelete (queue *q)
{
    pthread_mutex_destroy (q->mut);
    free (q->mut);      
    pthread_cond_destroy (q->notFull);
    free (q->notFull);
    pthread_cond_destroy (q->notEmpty);
    free (q->notEmpty);
    free (q->buf);
    free (q);
}
void queueAdd (queue *q, void* in)
{
    q->buf[q->tail] = in;
    q->tail++;
    if (q->tail == q->queue_size)
        q->tail = 0;
    if (q->tail == q->head)
        q->full = 1;
    q->empty = 0;

    return;
}
void queueDel (queue *q, void ** out)
{
    *out = q->buf[q->head];

    q->head++;
    if (q->head == q->queue_size)
        q->head = 0;
    if (q->head == q->tail)
        q->empty = 1;
    q->full = 0;
    return;
}

void getElement(queue * q, void ** out){
    pthread_mutex_lock (q->mut);
    while (q->empty) {
        pthread_cond_wait (q->notEmpty, q->mut);
    }
    
    queueDel (q, out);
    pthread_mutex_unlock (q->mut);
    pthread_cond_signal (q->notFull);
    return;
}

void addElement(queue * q, void * in){

    pthread_mutex_lock (q->mut);
    while (q->full) {
        pthread_cond_wait (q->notFull, q->mut);
    }
    queueAdd (q, in);
    pthread_mutex_unlock (q->mut);
    pthread_cond_signal (q->notEmpty);
    return;
}
