// See header for information
#include <ringbuffer.h>
#include <string.h>

void ringbuffer_init(ringbuffer *rb, int total_size, int element_size)
{
    rb->total_size = total_size;
    rb->element_size = element_size;
    rb->head = 0;
    rb->tail = 0;
}

int ringbuffer_is_full(volatile const ringbuffer *rb)
{
    return (rb->tail + rb->element_size) % rb->total_size == rb->head;
}

int ringbuffer_is_empty(volatile const ringbuffer *rb)
{
    return rb->tail == rb->head;
}

void ringbuffer_write(volatile ringbuffer *rb, const unsigned char *element)
{
    memcpy((void*)&rb->elems[rb->tail], (void*)element, rb->element_size);
    __sync_synchronize();
    rb->tail = (rb->tail + rb->element_size) % rb->total_size;
    if(ringbuffer_is_empty(rb)) rb->head = (rb->head + 1) % rb->total_size;
}

void ringbuffer_read(volatile ringbuffer *rb, unsigned char *element)
{
    memcpy((void*)element, (void*)&rb->elems[rb->head], rb->element_size);
    __sync_synchronize();
    rb->head = (rb->head + rb->element_size) % rb->total_size;
}
