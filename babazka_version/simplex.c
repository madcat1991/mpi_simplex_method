#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <pcontrol.h>

#define MASTER 0

int voidprintf(char *fmt, ...) {
}

#define log voidprintf

/*
 * Схема хранения симплекс-таблицы
 * ================================
 *  Для задачи
 * 
 *  maximize
 *  +2531920636.634x0 +4003819341.584x1 +1086267303.960x2

 *  subject to
 *  +40.058x0 + 6.151x1 +33.627x2 <= +2773.267
 *  + 6.151x0 +33.627x1 +29.675x2 <= +2049.243
 *  +33.627x0 +29.675x1 + 4.044x2 <= +2027.762
 *  +29.675x0 + 4.044x1 +32.659x2 <= +1963.574
 * в которой n=3 переменных и m=4 ограничений,
 * 
 * симплекс-таблица (одна большая матрица) есть следующий список из m+1 строк (каждая длины m+n+1):
 * [
 *  [0.0,                 -2531920636.6336632,  -4003819341.5841589,  -1086267303.9603961,  0.0,  0.0,  0.0,  0.0]
 *  [2773.2670739999999,  40.057912000000002,   6.1513150000000003,   33.627251000000001,   1.0,  0.0,  0.0,  0.0]
 *  [2049.242808,         6.1513150000000003,   33.627251000000001,   29.675108999999999,   0.0,  1.0,  0.0,  0.0]
 *  [2027.7617310000001,  33.627251000000001,   29.675108999999999,   4.0444529999999999,   0.0,  0.0,  1.0,  0.0]
 *  [1963.5742769999999,  29.675108999999999,   4.0444529999999999,   32.658912999999998,   0.0,  0.0,  0.0,  1.0]
 *  ]
 * 
 * Здесь мы храним таблицу по столбцам (т.е. tab[column_index][row_index]),
 * и при распараллеливании каждый узел получает свое подмножество столбцов для обработки. (+ нулевой столбец)
 */

/* симплекс-таблица */
struct SimplexTableau {
    int n;
    int n_cols;
    int m;
    int n_rows;
    /* номера переменных, отданных текущем узлу
     * (номер столбца = номер переменной + 1)
     */
    int x_lo;
    int x_hi;
    /* хранится по столбцам, причем каждый узел хранит
     * только столбец 0 + выбранный на текущей итерации
     * столбец + свой диапазон */
    double **tab;
    /* tab[column_index][row_index] */

    /* контекст MPI */
    int my_rank;
    int size;
};

/* предложение выбрать вот эту строку и столбец */
struct SimplexRowColumnVote {
    int my_rank;
    int p_xi;  /* i = p_xi + 1 */
    int p_row; /* j = row */
    double tab_xi_0; /* значение tab[p_xi + 1][0] */
};



#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


/* выбор диапазона, достающегося одному узлу для обработки  */
void select_chunk(int total_items, int i, int size, int *lo, int *hi) {
    int items_per_child = (int)(ceil((float)(total_items) / (float)(size)));
    *lo = items_per_child * i;
    *hi = MIN(total_items - 1, items_per_child * (i + 1) - 1);
}

/* инициализация, выделение памяти под симплекс-таблицу */
void spx_init_tableau(struct SimplexTableau *spx, int n, int m, int x_lo, int x_hi) {
    int i, j;
    spx->n = n;
    spx->m = m;
    spx->n_cols = n + m + 1;
    spx->n_rows = m + 1;

    spx->x_hi = x_hi;
    spx->x_lo = x_lo;

    spx->tab = malloc(sizeof(double *) * spx->n_cols);
    for (i = 0; i < spx->n_cols; i++) spx->tab[i] = NULL;
    spx->tab[0] = malloc(sizeof(double) * spx->n_rows);
    for (j = 0; j < spx->n_rows; j++) spx->tab[0][j] = 0.0;
    for (i = x_lo; i <= x_hi; i++) {
        spx->tab[i+1] = malloc(sizeof(double) * spx->n_rows);
        for (j = 0; j < spx->n_rows; j++) spx->tab[i+1][j] = 0.0;
    }
}

void spx_free_column(struct SimplexTableau *spx, int xi) {
    free(spx->tab[xi+1]);
    spx->tab[xi+1] = NULL;
}

/* считывание задачи линейного программирования из файла */
void spx_fill_from_file(struct SimplexTableau *spx, FILE *f, int my_rank, int size) {
    int m, n, i, j;
    int x_hi, x_lo;
    double f_d;

    spx->my_rank = my_rank;
    spx->size = size;

    fscanf(f, "%d", &n);
    fscanf(f, "%d", &m);

    select_chunk(m + n, my_rank, size, &x_lo, &x_hi);
    printf("rank %d: my chunk is %d..%d\n", my_rank, x_lo, x_hi);
    spx_init_tableau(spx, n, m, x_lo, x_hi);

    /* vector c */
    for (i = 0; i < n; i++) {
        fscanf(f, "%lf", &f_d);
        if (x_lo <= i && i <= x_hi) {
            spx->tab[i+1][0] = -f_d;
        }
    }
    for (j = 1; j <= m; j++) {
        /* row of matrix A */
        for (i = 0; i < n; i++) {
            fscanf(f, "%lf", &f_d);
            if (x_lo <= i && i <= x_hi) {
                spx->tab[i+1][j] = f_d;
            }
        }
        /* element of vector b */
        fscanf(f, "%lf", &f_d);
        spx->tab[0][j] = f_d;
        /* slack variables */
        i = n + j - 1;
        if (x_lo <= i && i <= x_hi) {
            spx->tab[i+1][j] = 1.0;
        }
    }
}

/* вывод симплекс-таблицы в stdout */
void spx_print(struct SimplexTableau *spx) {
    int i, j;
    double v;
    printf("[\n");
    for (j = 0; j < spx->n_rows; j++) {
        printf("[");
        for (i = 0; i < spx->n_cols; i++) {
            if (spx->tab[i]) {
                v = spx->tab[i][j];
            } else {
                v = 0.0;
            }
            if (v >= 0.0) {
                printf(" +%12.3lf", v);
            } else {
                printf(" -%12.3lf", -v);
            }
        }
        printf("]\n");
    }
    printf("]\n");
}

/* вывод решения (список иксов) */
void spx_print_xs(struct SimplexTableau *spx) {
    int xi;
    double x;
    int i, j;

    printf("[");
    for (xi = 0; xi < spx->n; xi++) {
        i = xi + 1;
        x = 0.0;
        if (!spx->tab[i]) {
            printf(" NO COL");
            continue;
        }

        if (spx->tab[i][0] == 0.0) {
            for (j = 1; j < spx->n_rows; j++) {
                x += spx->tab[0][j] * spx->tab[i][j];
            }
        } 
        printf(" %lf", x);
    }
    printf("]\n");
}

/* выбор строки и столбца для новой итерации */
struct SimplexRowColumnVote spx_select_col_and_row(struct SimplexTableau *spx) {
    struct SimplexRowColumnVote vote;
    int max_xi, i, j;
    double xij, ratio, min_ratio;
    vote.my_rank = spx->my_rank;
    vote.p_xi = -1;
    vote.p_row = -1;
    vote.tab_xi_0 = 0.0;

    max_xi = MIN(spx->n - 1, spx->x_hi);
    for (i = spx->x_lo; i <= max_xi; i++) {
        if (spx->tab[i+1][0] < vote.tab_xi_0) {
            vote.p_xi = i;
            vote.tab_xi_0 = spx->tab[i+1][0];
        }
    }

    if (vote.p_xi == -1) return vote;
    for (j = 1; j < spx->m + 1; j++) {
        xij = spx->tab[vote.p_xi+1][j];
        if (xij > 0) {
            ratio = spx->tab[0][j] / xij;
            if (vote.p_row == -1 || ratio < min_ratio) {
                vote.p_row = j;
                min_ratio = ratio;
            }
        }
    }
    return vote;
}


/* выбор наилучшего столбца по массиву голосов.
 * Если текущее решение оптимально, то поле my_rank результата равно -1.
*/
struct SimplexRowColumnVote spx_select_best_vote(struct SimplexRowColumnVote *votes, int n_votes) {
    struct SimplexRowColumnVote best;
    int i;
    best.my_rank = -1;
    best.tab_xi_0 = 1.0; /* выбираются столбцы с наименьшим отрицательным */
    for (i = 0; i < n_votes; i++) {
        if (votes[i].p_xi != -1 && votes[i].p_row != -1) {
            if (votes[i].tab_xi_0 < best.tab_xi_0) {
                best = votes[i];
            }
        }
    }
    return best;
}

/* внесение переменной p_xi в базис и исключение переменной p_row-1 из базиса */
void spx_pivot(struct SimplexTableau *spx, int p_xi, int p_row) {
    double pivot_value;
    double ratio;
    int xi, row;

    pivot_value = spx->tab[p_xi+1][p_row];
    for (xi = spx->x_lo; xi <= spx->x_hi; xi++) {
        spx->tab[xi+1][p_row] /= pivot_value;
    }
    spx->tab[0][p_row] /= pivot_value;

    #pragma omp parallel for private(row, ratio, xi)
    for (row = 0; row < spx->m + 1; row++) {
        if (row == p_row) continue;
        ratio = spx->tab[p_xi+1][row];
        for (xi = spx->x_lo; xi <= spx->x_hi; xi++) {
            spx->tab[xi+1][row] -= ratio * spx->tab[xi+1][p_row];
        }
        spx->tab[0][row] -= ratio * spx->tab[0][p_row];
    }

    if (!(spx->x_lo <= p_xi && p_xi <= spx->x_hi)) {
        spx_free_column(spx, p_xi);
    }
}

/* получение голоса от узла */
struct SimplexRowColumnVote par_receive_vote_unicast(int from_node) {
    struct SimplexRowColumnVote vote;
    MPI_Status status;
    int votesize = sizeof(struct SimplexRowColumnVote);
    int r;
    log("{0} recv vote from %d into %08x (%d bytes) \n", from_node, &vote, votesize);
    r = MPI_Recv(&vote, votesize, MPI_CHAR, from_node, 0, MPI_COMM_WORLD, &status);
    log("{0} recv returned %d\n", r);
    return vote;
}

/* отправка голоса узлу */
void par_send_vote_unicast(struct SimplexRowColumnVote vote, int to_node) {
    int votesize = sizeof(struct SimplexRowColumnVote);
    int r = MPI_Send(&vote, votesize, MPI_CHAR, to_node, 0, MPI_COMM_WORLD);
    log("[] send returned %d\n", r);
}

/* отправка или получение широковещательной рассылки с выбранным голосом */
void par_vote_rendezvous_broadcast(struct SimplexRowColumnVote *vote, int sender) {
    MPI_Bcast(vote, sizeof(struct SimplexRowColumnVote), MPI_CHAR, sender, MPI_COMM_WORLD);
}

/* отправка или получение широковещательной рассылки с одним столбцом */
void par_column_rendezvous_broadcast(struct SimplexTableau *spx, int col_i, int sender) {
    int row_size_in_bytes = spx->n_rows * sizeof(double);
    if (sender == spx->my_rank) {
        /* pass */
    } else {
        if (!(col_i == 0 || (spx->x_lo <= col_i - 1 && col_i - 1 <= spx->x_hi))) {
            spx->tab[col_i] = malloc(row_size_in_bytes);
        }
    }
    MPI_Bcast(spx->tab[col_i], row_size_in_bytes, MPI_CHAR, sender, MPI_COMM_WORLD);
}

/* получение недостающих столбцов со всех узлов */
void par_collect_all_columns(struct SimplexTableau *spx, int size) {
    MPI_Status status;
    int node_i;
    int row_size_in_bytes = spx->n_rows * sizeof(double);
    int x_lo, x_hi, col_i;
    for (node_i = 0; node_i < size; node_i++) {
        if (node_i == spx->my_rank) continue;
        select_chunk(spx->m + spx->n, node_i, size, &x_lo, &x_hi);
        log("{0} getting cols %d..%d from node %d\n", x_lo+1, x_hi+1, node_i);
        for (col_i = x_lo + 1; col_i <= x_hi + 1; col_i++) {
            spx->tab[col_i] = malloc(row_size_in_bytes);
            MPI_Recv(spx->tab[col_i], row_size_in_bytes, MPI_CHAR, node_i, 0, MPI_COMM_WORLD, &status);
            log("{0} got col %d\n", col_i);
        }
    }
    for (col_i = 0; col_i < spx->n_cols; col_i++) {
        if (!spx->tab[col_i]) {
            printf("{0} achtung! column %d is missing!\n", col_i);
        }
    }
}

/* рассылка своих столбцов */
void par_send_my_columns(struct SimplexTableau *spx, int receiver) {
    int col_i;
    int row_size_in_bytes = spx->n_rows * sizeof(double);
    for (col_i = spx->x_lo + 1; col_i <= spx->x_hi + 1; col_i++) {
        MPI_Send(spx->tab[col_i], row_size_in_bytes, MPI_CHAR, receiver, 0, MPI_COMM_WORLD);
    }
}

/* что делает главный узел */
void par_master(struct SimplexTableau *spx, int rank, int size) {
    struct SimplexRowColumnVote *votes;
    struct SimplexRowColumnVote best_vote;
    int node_i;
    int step_i;
    double t0, t1, t2;

    votes = malloc(size * sizeof(struct SimplexRowColumnVote));

    log("{0} entered par_master\n");
    t0 = MPI_Wtime();

    for (step_i = 0; ; step_i++) {
        votes[MASTER] = spx_select_col_and_row(spx);
        log("{0} my vote is for (%d %d)\n", votes[MASTER].p_xi, votes[MASTER].p_row);

        MPI_Pcontrol(TRACEEVENT, "entry", 0, 0, "");
        for (node_i = MASTER + 1; node_i < size; node_i++) {
            votes[node_i] = par_receive_vote_unicast(node_i);
        }
        MPI_Pcontrol(TRACEEVENT, "exit", 0, 0, "");

        best_vote = spx_select_best_vote(votes, size);
        log("{0} reached decision (%d %d)\n", best_vote.p_xi, best_vote.p_row);

        /* это мы рассылаем */
        MPI_Pcontrol(TRACEEVENT, "entry", 1, 0, "");
        par_vote_rendezvous_broadcast(&best_vote, MASTER);
        log("{0} transmitted decision\n");

        if (best_vote.my_rank == -1) {
            printf("step %d: optimal\n", step_i);
            break;
        }

        printf("step %d: x %d, row %d  (node %d)\n", step_i, best_vote.p_xi, best_vote.p_row, best_vote.my_rank);
        /* а это мы получаем или рассылаем, в зависимости от выбранного голоса */
        par_column_rendezvous_broadcast(spx, best_vote.p_xi + 1, best_vote.my_rank);
        log("{0} synced column %d\n", best_vote.p_xi + 1);
        MPI_Pcontrol(TRACEEVENT, "exit", 1, 0, "");

        spx_pivot(spx, best_vote.p_xi, best_vote.p_row);
        log("{0} completed pivot\n");
    }
    t1 = MPI_Wtime();
    /* соберем итоговую таблицу по частям */
    MPI_Pcontrol(TRACEEVENT, "entry", 2, 0, "");
    par_collect_all_columns(spx, size);
    MPI_Pcontrol(TRACEEVENT, "exit", 2, 0, "");
    t2 = MPI_Wtime();

    /*spx_print(spx); */
    spx_print_xs(spx);
    printf("F = %lf\n", spx->tab[0][0]);

    printf("calc time:    %20.5lf\n", t1-t0);
    printf("compose time: %20.5lf\n", t2-t1);
    printf("total time:   %20.5lf\n", t2-t0);

    free(votes);
}

/* что делает подчиненный узел */
void par_slave(struct SimplexTableau *spx, int rank, int size) {
    struct SimplexRowColumnVote vote;
    int step_i;
    log("[%d] entered par_slave\n", rank);

    for (step_i = 0; ; step_i++) {
        vote = spx_select_col_and_row(spx);
        log("[%d] my vote is for (%d %d), sending to master \n", rank, vote.p_xi, vote.p_row);
        MPI_Pcontrol(TRACEEVENT, "enter", 0, 0, "");
        par_send_vote_unicast(vote, MASTER);
        MPI_Pcontrol(TRACEEVENT, "exit", 0, 0, "");

        MPI_Pcontrol(TRACEEVENT, "enter", 1, 0, "");
        par_vote_rendezvous_broadcast(&vote, MASTER);
        log("[%d] received decision (%d %d)\n", rank, vote.p_xi, vote.p_row);

        if (vote.my_rank == -1) {
            break;
        }
        /* а это мы получаем или рассылаем, в зависимости от выбранного голоса */
        par_column_rendezvous_broadcast(spx, vote.p_xi + 1, vote.my_rank);
        log("[%d] synced column %d\n", rank, vote.p_xi + 1);
        MPI_Pcontrol(TRACEEVENT, "exit", 1, 0, "");

        spx_pivot(spx, vote.p_xi, vote.p_row);
        log("[%d] completed pivot\n", rank);
    }
    /* соберем итоговую таблицу по частям */
    log("[%d] sending columns\n", rank);
    MPI_Pcontrol(TRACEEVENT, "entry", 2, 0, "");
    par_send_my_columns(spx, MASTER);
    MPI_Pcontrol(TRACEEVENT, "exit", 2, 0, "");
    log("[%d] terminating\n", rank);
}

void main_simplex(int rank, int size, char *filename) {
    struct SimplexTableau spx;
    FILE *f;
    double t_0, t_1;

    t_0 = MPI_Wtime();

    f = fopen(filename, "r");
    if (!f) {
        perror("fopen"); return;
    }
    spx_fill_from_file(&spx, f, rank, size);
    fclose(f);

    printf("My rank is %d of %d\n", rank, size);
    if (rank == MASTER) {
        par_master(&spx, rank, size);
    } else {
        par_slave(&spx, rank, size);
    }

    t_1 = MPI_Wtime();

    if (rank == MASTER) {
        printf("run time:   %20.5lf\n", t_1 - t_0);
    }
}

int main(int argc, char **argv) {
	MPI_Init(NULL, NULL);
    MPI_Pcontrol(TRACEFILES, NULL, "trace", 0);
    MPI_Pcontrol(TRACELEVEL, 1, 1, 1);
    MPI_Pcontrol(TRACENODE, 1000000, 1, 1);

	int rank, size;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc == 1 || !argv[1]) {
        printf("Usage: %s <simplex_task.txt>\n", argv[0]);
    } else {
        main_simplex(rank, size, argv[1]);
    }

    MPI_Finalize();
}
