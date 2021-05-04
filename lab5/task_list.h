#ifndef LAB5_TASK_LIST_H
#define LAB5_TASK_LIST_H

#include <stddef.h>
#include <pthread.h>
#include <stdbool.h>

typedef struct TaskNode {
    size_t task;
    struct TaskNode* next;
} TaskNode;

typedef struct {
    TaskNode* head;
    TaskNode* tail;

    size_t size;

    pthread_mutex_t lock;
} TaskList;

bool tl_init(TaskList* this);
void tl_free(TaskList* this);

bool tl_add(TaskList* this, size_t task);
bool tl_pop(TaskList* this, size_t* task);

size_t tl_size(TaskList* this);

#endif // !LAB5_TASK_LIST_H