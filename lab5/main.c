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

#define REQUEST_COUNT 1
#define STARVING_THRESHOLD 1
#define PING_THRESHOLD (STARVING_THRESHOLD)
#define L 2000
#define N 1000
#define ITERATIONS 1

double load_cpu(size_t i) {
    return sin(cos(2.123 * i / (i + 0.0123432)));
}

void make_iterations(size_t count) {
    double result = 0.0;
    for (size_t i = 0; i < count; ++i) {
        result = load_cpu(i);
    }
    printf("%f ", result);
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

void log_on_root(const char* format, ...) {
    if (get_rank(MPI_COMM_WORLD) != 0) return;

    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
}

#define REQUEST_TAG 1
#define JOB_SEND_TAG 2
#define COUNT_SEND_TAG 3

#define REQUEST_JOB 0
#define REQUEST_CANCEL 1
#define COUNT_REQUEST 2

#define ERR_TASK LONG_MAX

typedef struct {
    TaskList list;
    size_t proc_count;
    pthread_cond_t start_requests;
    pthread_mutex_t filling_requests[3];
    bool tasks_available;
} TaskData;

bool td_init(TaskData* this) {
    this->proc_count = get_size(MPI_COMM_WORLD);

    if (!tl_init(&this->list, STARVING_THRESHOLD)) return false;

    int err;
    if ((err = pthread_cond_init(&this->start_requests, NULL)) != 0) {
        fprintf(stderr, "pthread_cond_init: %s\n", strerror(err));
        return false;
    }

    for (size_t i = 0; i < 3; ++i) {
        if ((err = pthread_mutex_init(&this->filling_requests[i], NULL)) != 0) {
            pthread_cond_destroy(&this->start_requests);

            for (size_t j = 0; j < i; ++j) {
                pthread_mutex_destroy(&this->filling_requests[j]);
            }

            fprintf(stderr, "pthread_mutex_init: %s\n", strerror(err));
            return false;
        }
    }

    this->tasks_available = true;

    return true;
}

void td_free(TaskData* this) {
    tl_free(&this->list);
    pthread_cond_destroy(&this->start_requests);
    for (size_t i = 0; i < 3; ++i) {
        pthread_mutex_destroy(&this->filling_requests[i]);
    }
}

bool td_refresh(TaskData* this, size_t iteration) {
    const int rank = get_rank(MPI_COMM_WORLD);
    const int size = get_size(MPI_COMM_WORLD);

    size_t factor = 0;
    for (size_t i = 0; i < size; ++i) {
        factor += i;
    }

    for (size_t i = 0; i < N; ++i) {
        tl_add(&this->list, 
            //(size_t) L * abs((N / 2) - i) * size * (abs(rank - (iteration % size)) / (double) factor));
            (size_t) L * abs((N / 2) - i) * (abs(rank - (iteration % size)) / (double) factor));
    }

    this->tasks_available = true;
    tl_reset_starving(&this->list);

    return true;
}

void* accept_loop(TaskData* this) {
    for (;;) {
        size_t request;
        MPI_Status status;
        MPI_Recv(&request, 1, MPI_LONG, MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        switch (request) {
            case REQUEST_JOB: {
                size_t task;
                if (tl_starving(&this->list)) {
                    task = ERR_TASK;
                } else {
                    if (!tl_pop(&this->list, &task)) {
                        task = ERR_TASK;
                    }
                }

                MPI_Send(&task, 1, MPI_LONG, status.MPI_SOURCE, JOB_SEND_TAG, MPI_COMM_WORLD);
                break;
            }

            case COUNT_REQUEST: {
                size_t count = tl_size(&this->list);
                if (tl_starving(&this->list)) {
                    count = 0;
                }
                MPI_Send(&count, 1, MPI_LONG, status.MPI_SOURCE, COUNT_SEND_TAG, MPI_COMM_WORLD);
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

typedef struct {
    int rank;
    size_t weight;
} WeightedRank;

int compare_ranks(WeightedRank* lhs, WeightedRank* rhs) {
    if (rhs->weight > lhs->weight) return 1;
    if (rhs->weight < lhs->weight) return -1;
    return 0;
}

size_t try_make_requests(size_t* tasks, size_t count) {
    MPI_Status unused;

    const size_t rank = get_rank(MPI_COMM_WORLD);
    const size_t proc_count = get_size(MPI_COMM_WORLD);

    WeightedRank* ordered_ranks = malloc((proc_count - 1) * sizeof(*ordered_ranks));

    size_t j = 0;
    for (size_t i = 0; i < proc_count; ++i) {
        if (i == rank) continue;

        size_t job_flag = COUNT_REQUEST;
        MPI_Send(&job_flag, 1, MPI_LONG, i, REQUEST_TAG, MPI_COMM_WORLD);
        MPI_Recv(&ordered_ranks[j].weight, 1, MPI_LONG, i, COUNT_SEND_TAG, MPI_COMM_WORLD, &unused);
        ordered_ranks[j].rank = i;
        j++;
    }

    qsort(ordered_ranks, proc_count - 1, sizeof(*ordered_ranks), (__compar_fn_t) compare_ranks);

    size_t found_count = 0;
    for (size_t i = 0; i < proc_count - 1 && found_count < count;) {
        size_t job_flag = REQUEST_JOB;
        MPI_Send(&job_flag, 1, MPI_LONG, ordered_ranks[i].rank, REQUEST_TAG, MPI_COMM_WORLD);
        MPI_Recv(&tasks[found_count], 1, MPI_LONG, ordered_ranks[i].rank, JOB_SEND_TAG, MPI_COMM_WORLD, &unused);

        if (tasks[found_count] != ERR_TASK) {
            found_count++;
        } else {
            i++;
        }
    }

    free(ordered_ranks);
    return found_count;
}

size_t mod_inc(size_t num, size_t mod) {
    return (num + 1) % mod;
} 

size_t mod_dec(size_t num, size_t mod) {
    return (num - 1) % mod;
}

bool fetch_requests(TaskData* this) {
    static size_t tasks[REQUEST_COUNT];
    size_t count;
    if ((count = try_make_requests(tasks, REQUEST_COUNT)) == 0) {
        return false;
    }

    for (size_t i = 0; i < count; ++i) {
        tl_add(&this->list, tasks[i]);
    }

    return true;
}

void* request_loop(TaskData* this) {
    double actual_time = 0.0;
    
    size_t my_lock = 1;
    pthread_mutex_lock(&this->filling_requests[0]);
    pthread_mutex_lock(&this->filling_requests[1]);
    pthread_mutex_unlock(&this->filling_requests[0]);

    pthread_cond_signal(&this->start_requests);
    pthread_mutex_lock(&this->filling_requests[2]);

    for (;;) {
        const double start = MPI_Wtime();

        this->tasks_available = false;
        bool tried_requests = false;

        while (tl_size(&this->list) <= STARVING_THRESHOLD) {
            tried_requests = true;
            if (!fetch_requests(this)) {
                break;
            }
            this->tasks_available = true;
        }

        if (!tried_requests) {
            this->tasks_available = true;
        }
        
        const double end = MPI_Wtime();
        actual_time += (end - start);

        pthread_mutex_unlock(&this->filling_requests[my_lock]);
        my_lock = mod_inc(my_lock, 3);
        pthread_mutex_lock(&this->filling_requests[mod_inc(my_lock, 3)]);

        if (!this->tasks_available) {
            pthread_mutex_unlock(&this->filling_requests[my_lock]);
            pthread_mutex_unlock(&this->filling_requests[mod_inc(my_lock, 3)]);
            fprintf(stderr, "%d: request thread time: %f\n", get_rank(MPI_COMM_WORLD), actual_time);
            return NULL;
        }
    }

    return NULL;
}

size_t task_loop(TaskData* this) {
    size_t tasks_done = 0;
    size_t my_lock = 0;

    for (;;) {
        size_t task;
        if (tl_pop(&this->list, &task)) {
            tasks_done++;
            make_iterations(task);
        }

        if (tl_size(&this->list) <= PING_THRESHOLD) {
            if (tl_size(&this->list) == 0) {
                pthread_mutex_lock(&this->filling_requests[mod_inc(my_lock, 3)]);

                if (!this->tasks_available) {
                    pthread_mutex_unlock(&this->filling_requests[mod_inc(my_lock, 3)]);
                    pthread_mutex_unlock(&this->filling_requests[my_lock]);
                    break;
                }
            } else {
                if (pthread_mutex_trylock(&this->filling_requests[mod_inc(my_lock, 3)]) != 0) {
                    continue; 
                }
            }

            pthread_mutex_unlock(&this->filling_requests[my_lock]);
            my_lock = mod_inc(my_lock, 3);
        }            
    }

    return tasks_done;
}

void print_disbalance(double time, size_t iteration) {
    const size_t size = get_size(MPI_COMM_WORLD);
    double* proc_times = NULL;
    if (get_rank(MPI_COMM_WORLD) == 0) {
        proc_times = malloc(size * sizeof(*proc_times));
    }

    MPI_Gather(&time, 1, MPI_DOUBLE, proc_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (get_rank(MPI_COMM_WORLD) == 0) {
        double max_delta = 0.0;
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                const double delta = fabs(proc_times[i] - proc_times[j]);
                if (delta > max_delta) {
                    max_delta = delta;
                }
            }
        }

        double max_time = 0.0;
        for (size_t i = 0; i < size; ++i) {
            if (proc_times[i] > max_time) {
                max_time = proc_times[i];
            }
        }

        fprintf(stderr, "Iteration %d: max disbalance: %f; disbalance percentage: %.2f%%\n", 
            iteration, max_delta,
            max_delta / max_time * 100.0);
    }

    free(proc_times);
}

int main_loop() {
    TaskData data;
    if (!td_init(&data)) return EXIT_FAILURE;

    pthread_t accept_thread;
    pthread_create(&accept_thread, NULL, (void* (*)(void*)) accept_loop, &data);

    for (size_t i = 0; i < ITERATIONS; ++i) {
        if (!td_refresh(&data, i)) {
            return EXIT_FAILURE;
        }
        MPI_Barrier(MPI_COMM_WORLD);


        pthread_mutex_lock(&data.filling_requests[0]);
        pthread_mutex_lock(&data.filling_requests[2]);
        pthread_t request_thread;
        pthread_create(&request_thread, NULL, (void* (*)(void*)) request_loop, &data);

        pthread_cond_wait(&data.start_requests, &data.filling_requests[0]);
        pthread_mutex_unlock(&data.filling_requests[2]);

        const double start = MPI_Wtime();
        const size_t tasks_done = task_loop(&data);
        const double end = MPI_Wtime();

        pthread_join(request_thread, NULL);
        fprintf(stderr, "Process %d: Iteration %d: Tasks: %d: Time: %f\n", 
            get_rank(MPI_COMM_WORLD), i, tasks_done, end - start);
        MPI_Barrier(MPI_COMM_WORLD);
        print_disbalance(end - start, i);
    }

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
        log_on_root("MPI_THREAD_MULTIPLE is not supported!\n");
        return EXIT_FAILURE;
    }

    return main_loop();
}
