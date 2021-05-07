#include "task_list.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

TaskNode* tn_init(size_t task) {
    TaskNode* this = malloc(sizeof(*this));
    if (this == NULL) {
        perror("malloc");
        return NULL;
    }

    memcpy(&this->task, &task, sizeof(task));
    this->next = NULL;
    return this;
}

bool tl_init(TaskList* this, size_t starving_threshold) {
    if (this == NULL) return false;
    memset(this, 0, sizeof(*this));

    int err;
    if ((err = pthread_mutex_init(&this->lock, NULL)) != 0) {
        fprintf(stderr, "pthread_mutex_init: %s\n", strerror(err));
        return false;
    }

    this->starving_threshold = starving_threshold;
    this->starving = false;
    
    return true;
}

bool tl_add(TaskList* this, size_t task) {
    if (this == NULL) return false;

    TaskNode* node = tn_init(task);
    if (node == NULL) return false;

    pthread_mutex_lock(&this->lock);

    if (this->size == 0) {
        this->head = node;
        this->tail = node;
        this->size++;
        pthread_mutex_unlock(&this->lock);
        return true;
    }

    this->tail->next = node;
    this->tail = this->tail->next;
    this->tail->next = NULL;
    this->size++;

    pthread_mutex_unlock(&this->lock);

    return true;
}

static bool remove_unsafe(TaskList* this) {
    if (this->head == NULL) return false;

    TaskNode* new_head = this->head->next;
    free(this->head);
    this->head = new_head;

    this->size--;

    if (this->size == 0) {
        this->tail = NULL;
    }

    return true;
}

bool tl_pop(TaskList* this, size_t* task) {
    if (this == NULL) return false;

    pthread_mutex_lock(&this->lock);

    if (this->head == NULL) {
        pthread_mutex_unlock(&this->lock);
        return false;
    }

    memcpy(task, &this->head->task, sizeof(*task));
    remove_unsafe(this);

    if (this->size <= this->starving_threshold) this->starving = true;

    pthread_mutex_unlock(&this->lock);

    return true;
}

void tl_free(TaskList* this) {
    pthread_mutex_lock(&this->lock);

    while (remove_unsafe(this)) 
        ;

    pthread_mutex_unlock(&this->lock);

    pthread_mutex_destroy(&this->lock);
}

size_t tl_size(TaskList* this) {
    pthread_mutex_lock(&this->lock);
    const size_t size = this->size;
    pthread_mutex_unlock(&this->lock);

    return size;
}

bool tl_starving(TaskList* this) {    
    pthread_mutex_lock(&this->lock);
    const bool starving = this->starving;
    pthread_mutex_unlock(&this->lock);

    return starving;
}

void tl_reset_starving(TaskList* this) {
    pthread_mutex_lock(&this->lock);
    this->starving = false;
    pthread_mutex_unlock(&this->lock);
}
