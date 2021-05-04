#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdarg.h>
#include <limits.h>

#include "task_list.h"

double load_cpu(size_t i) {
    return sin(cos(2.123 * i / (i + 0.0123432)));
}

void make_iterations(size_t count) {
    for (size_t i = 0; i < count; ++i) {
        printf("%f ", load_cpu(i));
    }
}

bool generate_iteration(TaskList* list, size_t min, size_t max) {
    if (max <= min) {
        fprintf(stderr, "max (%d) must be > min (%d)\n", max, min);
        return false;
    }

    const size_t task = min + (rand() % (max - min));
    if (!tl_add(list, task)) return false;

    return true;
}

bool generate_iterations(TaskList* list, size_t min, size_t max, size_t task_count) {
    for (size_t i = 0; i < task_count; ++i) {
        if (!generate_iteration(list, min, max)) return false;
    }

    return true;
}

int get_rank(MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

int get_size(MPI_Comm comm) {
    int size;
    MPI_Comm_size(comm, &size);
    return size;
}

void debug_on_root(const char* format, ...) {
    if (get_rank(MPI_COMM_WORLD) != 0) return;

    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
}

#define REQUEST_TAG 1
#define JOB_SEND_TAG 2

#define REQUEST_JOB 0
#define REQUEST_NO_JOBS 1
#define REQUEST_CANCEL 2

#define ERR_TASK LONG_MAX

typedef struct {
    TaskList list;
    bool* has_jobs;
    size_t proc_count;
    pthread_mutex_t lock;
} TaskData;

bool td_init(TaskData* this) {
    this->proc_count = get_size(MPI_COMM_WORLD);
    this->has_jobs = calloc(this->proc_count, sizeof(*this->has_jobs));

    if (this->has_jobs == NULL) {
        perror("calloc");
        return false;
    }

    if (!tl_init(&this->list)) {
        free(this->has_jobs);
        return false;
    }

    int err;
    if ((err = pthread_mutex_init(&this->lock, NULL)) != 0) {
        fprintf(stderr, "pthread_mutex_init: %s\n", strerror(err));
        free(this->has_jobs);
        tl_free(&this->list);
        return false;
    }

    return true;
}

void td_free(TaskData* this) {
    free(this->has_jobs);
    tl_free(&this->list);
    pthread_mutex_destroy(&this->lock);
}

bool td_refresh(TaskData* this) {
    MPI_Barrier(MPI_COMM_WORLD);

    pthread_mutex_lock(&this->lock);

    //if (!generate_iterations(&this->list, 100, 5000000, 500)) return false;
    for (size_t i = 0; i < 50; ++i) {
        if (get_rank(MPI_COMM_WORLD) == 0) {
            tl_add(&this->list, 500000);
        } else {
            tl_add(&this->list, 50000);
        }
    }

    for (size_t i = 0; i < this->proc_count; ++i) {
        this->has_jobs[i] = true;
    }

    pthread_mutex_unlock(&this->lock);

    return true;
}

bool td_is_empty(TaskData* this) {
    pthread_mutex_lock(&this->lock);

    for (size_t i = 0; i < this->proc_count; ++i) {
        if (this->has_jobs[i]) {
            pthread_mutex_unlock(&this->lock);
            return false;
        }
    }

    pthread_mutex_unlock(&this->lock);
    return true;
}

bool td_has_jobs(TaskData* this, size_t i) {
    pthread_mutex_lock(&this->lock);
    const bool has_jobs = this->has_jobs[i];
    pthread_mutex_unlock(&this->lock);

    return has_jobs;
}

void td_set_has_jobs(TaskData* this, size_t i, bool has_jobs) {
    pthread_mutex_lock(&this->lock);
    this->has_jobs[i] = has_jobs;
    pthread_mutex_unlock(&this->lock);
}

bool request_jobs(TaskData* this, size_t* task) {
    if (td_is_empty(this)) return false;

    const size_t rank = get_rank(MPI_COMM_WORLD);
    for (size_t i = 0; i < this->proc_count; ++i) {
        if (i == rank) continue;
        if (!td_has_jobs(this, i)) continue;

        size_t job_flag = REQUEST_JOB;
        MPI_Send(&job_flag, 1, MPI_LONG, i, REQUEST_TAG, MPI_COMM_WORLD);
        MPI_Recv(task, 1, MPI_LONG, i, JOB_SEND_TAG, MPI_COMM_WORLD, NULL);
        if (*task == ERR_TASK) {
        } else {
            return true;
        }
    }

    return false;
}

void notify_no_jobs(TaskData* this) {
    const size_t rank = get_rank(MPI_COMM_WORLD);
    td_set_has_jobs(this, rank, false);
    for (size_t i = 0; i < this->proc_count; ++i) {
        if (i == rank) continue;

        size_t no_job_flag = REQUEST_NO_JOBS;
        MPI_Send(&no_job_flag, 1, MPI_LONG, i, REQUEST_TAG, MPI_COMM_WORLD);
    }
}

void* accept_loop(TaskData* this) {
    for (;;) {
        size_t request;
        MPI_Status status;
        MPI_Recv(&request, 1, MPI_LONG, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        fprintf(stderr, "request: (%d, %d)\n", status.MPI_SOURCE, request);
        switch (request) {
            case REQUEST_JOB: {
                size_t task;
                if (!tl_pop(&this->list, &task)) {
                    task = ERR_TASK;
                }
                MPI_Send(&task, 1, MPI_LONG, status.MPI_SOURCE, JOB_SEND_TAG, MPI_COMM_WORLD);
                break;
            }
        
            case REQUEST_NO_JOBS: {
                td_set_has_jobs(this, status.MPI_SOURCE, false);
                break;
            } 

            case REQUEST_CANCEL: {
                return NULL;
            }

            default: break;
        }
    }

    return NULL;
}

int main_loop() {
    TaskData data;
    if (!td_init(&data)) return EXIT_FAILURE;

    pthread_t accept_thread;
    pthread_create(&accept_thread, NULL, (void* (*)(void*)) accept_loop, &data);

    for (size_t i = 0; i < 3; ++i) {
        if (!td_refresh(&data)) return EXIT_FAILURE;

        size_t task;
        size_t j = 0;

        while (tl_pop(&data.list, &task)) {
            make_iterations(task);
        }

        notify_no_jobs(&data);

        while (request_jobs(&data, &task)) {
            make_iterations(task);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    size_t cancel_flag = REQUEST_CANCEL;
    MPI_Send(&cancel_flag, 1, MPI_LONG, get_rank(MPI_COMM_WORLD), REQUEST_TAG, MPI_COMM_WORLD);
    pthread_join(accept_thread, NULL);

    td_free(&data);
    return EXIT_SUCCESS;
}

void finalize_wrapper() {
    MPI_Finalize();
}

int main(int argc, char* argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    atexit(finalize_wrapper);

    if (provided != MPI_THREAD_MULTIPLE) {
        debug_on_root("MPI_THREAD_MULTIPLE is not supported!\n");
        return EXIT_FAILURE;
    }

    return main_loop();
}
