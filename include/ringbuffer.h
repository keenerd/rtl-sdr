/*
 * Lockfree high-throughput inter-thread ringbuffer.
 * Written by Elias Önal <EliasOenal@gmail.com>, released as public domain.
 */

#ifndef RINGBUFFER_H
#define RINGBUFFER_H

#include <stdint.h>
#include <signal.h>

typedef struct {
    volatile sig_atomic_t total_size;
    volatile sig_atomic_t element_size;
    volatile sig_atomic_t head;
    volatile sig_atomic_t tail;
    volatile unsigned char elems[1];
} ringbuffer;

void ringbuffer_init(ringbuffer *rb, int total_size, int element_size);
int ringbuffer_is_empty(volatile const ringbuffer *rb);
int ringbuffer_is_full(volatile const ringbuffer *rb);
void ringbuffer_write(volatile ringbuffer *rb, const unsigned char *element);
void ringbuffer_read(volatile ringbuffer *rb, unsigned char *element);

#endif
