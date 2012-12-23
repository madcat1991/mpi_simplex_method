#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <float.h>
#include <time.h>
#include <pcontrol.h>
#include <math.h>

#define GEN_RANGE_MAX 100
#define GEN_RANGE_MIN -100

/*
 * Структура для хранения таблицы симплекс метода, представлена следующим образом
 * row_number - количество строк
 * column_number - количество столбцов
 *         0   1  ... n-1  n  ...  n+m-1    
 *  -------------------------------------------------
 *  |    |x1 |x2 |...|xn |xn+1|...|xn+m |b    |z    |
 *  -------------------------------------------------
 *0 |    |-c1|-c2|...|-cn|0   |...|0    |0    |1    |
 *  -------------------------------------------------
 *1 |xn+1|a11|a12|...|a1n|1   |...|0    |b1   |0    |
 *2 |xn+2|a21|a22|...|a2n|0   |...|0    |b2   |0    |
 *. |....|...|...|...|...|....|...|.... |.....|.....|
 *m |xn+m|am1|am2|...|amn|0   |...|1    |bm   |0    |
 *  -------------------------------------------------
 * */
typedef struct{
    int n;
    int m;

    double **table;
    double *b;
    double *z;
} SimplexTable;

typedef struct{
    double **C;
    double *b;

    int row_count;
    int col_count;
} ProcessInitData;

typedef struct{
    int row_no;
    double min_X;
} ProcessWorkResult;

SimplexTable init_random_simplex_table(int n, int m)
{
    //Генерируем рандомную симплех таблицу с количеством переменных n 
    //и количеством услови m
    srand(time(0));

    printf("Initing random simplex method table\n");

    int row_number, column_number;
    //+1 под строку с C
    row_number = m + 1;
    //+2 под столбцы:
    //2-ой с конца, для b, 
    //1-ый с конца, для z 
    column_number = n + m;

    double **table;
    table = malloc(sizeof(double *) * row_number);
 
    double *b, *z;
    b = malloc(sizeof(double) * row_number);
    z = malloc(sizeof(double) * row_number);
    int i, j;
    double sum_sqrs;

    //генерация table
    for(i = 0; i < row_number; i++){
        table[i] = malloc(sizeof(double) * column_number);

        for(j = 0; j < column_number; j++){
            if(j < n) {
                //то что относиться к оригинальным переменным - рандом
                // первая строка инвертирована
                double rand_value = (double)rand() / (GEN_RANGE_MAX + 1) * (GEN_RANGE_MAX - GEN_RANGE_MIN) + GEN_RANGE_MIN;
                sum_sqrs += rand_value * rand_value;

                if(i != 0)
                {
                    table[i][j] = rand_value;
                } else {
                    table[i][j] = -rand_value;
                }
            } else {
                //то что относиться к введенным переменным(базисным)
                //единички по диагонали, кроме первой строки, в которой оставшиеся
                //значения - 0
                if(j - n == i - 1){
                    table[i][j] = 1;
                } else {
                    table[i][j] = 0;
                }
            }
        }

        if(i != 0){
            double rand_value = (double)rand() / (GEN_RANGE_MAX + 1) * (GEN_RANGE_MAX - GEN_RANGE_MIN) + GEN_RANGE_MIN;
            sum_sqrs += rand_value * rand_value;

            b[i] = rand_value;
            z[i] = 0;
        } else {
            z[i] = 1;
        }
    }

    double scale = sqrt(sum_sqrs) / 100.0;
    for(i = 0; i < row_number; i++){
        for(j = 0; j < column_number; j++){
            table[i][j] /= scale;
        }

        b[i] /= scale;
    }
 
    SimplexTable st = {
        .n = n,
        .m = m,
        .table = table,
        .b = b,
        .z = z
    };

    printf("Table has initialized\n");
    return st;
}


SimplexTable init_simplex_table(){
    //Min z = (c^(T)) * x
    //Ax = b
    //x >= 0

    //x = (x0, x1)
    //c = (2, 5)
    //A = (1, 0; 0, 1; 1, 1)
    //b = (40, 30, 50)
    //len(x) = n, n = 2
    //len(b) = m, m = 3

    printf("Initing simplex method table\n");
    double **table;
    int n = 2, m = 3; //n - количество переменных, m - количество условий

    int row_number, column_number;
    //+1 под строку с C
    row_number = m + 1;
    //+2 под столбцы:
    //2-ой с конца, для b, 
    //1-ый с конца, для z 
    column_number = n + m;
    
    table = malloc(sizeof(double *) * row_number);
    int i;
    for(i = 0; i < row_number; i++){
        table[i] = malloc(sizeof(double) * column_number);
    }

    //нулевая строка
    table[0][0] = -2;
    table[0][1] = -5;
    table[0][2] = 0;
    table[0][3] = 0;
    table[0][4] = 0;

    //первая строка
    table[1][0] = 1;
    table[1][1] = 0;
    table[1][2] = 1;
    table[1][3] = 0;
    table[1][4] = 0;

    //вторая строка
    table[2][0] = 0;
    table[2][1] = 1;
    table[2][2] = 0;
    table[2][3] = 1;
    table[2][4] = 0;
    
    //третья строка
    table[3][0] = 1;
    table[3][1] = 1;
    table[3][2] = 0;
    table[3][3] = 0;
    table[3][4] = 1;

    double *b, *z;
    b = malloc(sizeof(double) * row_number);
    b[0] = 0; b[1] = 40; b[2] = 30; b[3] = 50;
    
    z = malloc(sizeof(double) * row_number);
    z[0] = 1; z[1] = 0; z[2] = 0; z[3] = 0;

    SimplexTable st = {
        .n = n,
        .m = m,
        .table = table,
        .b = b,
        .z = z
    };

    printf("Table has initialized\n");
    return st;
}

ProcessInitData get_process_init_data(SimplexTable st, int max_rows_per_proc,
    int row_length, int worker_index)
{
    printf("Initialize process %d\n", worker_index);
    ProcessInitData pid;
    int i;

    int row_count_for_process = max_rows_per_proc;
    /*
     * Так как воркеры нумеруются с 0, а строк может
     * быть нечетное количество
     */
    if ((worker_index + 1) * max_rows_per_proc > st.m){
        row_count_for_process -= 1;
    }

    pid.row_count = row_count_for_process;
    pid.col_count = row_length;
    pid.C = malloc(sizeof(double *) * row_count_for_process);
    pid.b = malloc(sizeof(double) * row_count_for_process);
    
    #pragma omp parallel for private(worker_index)
    for(i = 0; i < pid.row_count; i++)
    {
        pid.C[i] = malloc(sizeof(double) * (pid.col_count));
        //+1 так как надо пропустить строку с c^T
        int row_index = worker_index * max_rows_per_proc + i + 1;
        
        int j;
        for(j = 0; j < row_length; j++){
            pid.C[i][j] = st.table[row_index][j];
        }
        pid.b[i] = st.b[row_index];
    }

    printf("Process %d initialized\n", worker_index);
    return pid;
}

void send_process_init_data(ProcessInitData data, int reciever_index){
    //отправляем данные на инициализацию процессов
    MPI_Send(&(data.row_count), 1, MPI_INT, reciever_index, 0, MPI_COMM_WORLD);
    MPI_Send(&(data.col_count), 1, MPI_INT, reciever_index, 0, MPI_COMM_WORLD);

    int i;
    for(i = 0; i < data.row_count; i++){
        MPI_Send(data.C[i], data.col_count, MPI_DOUBLE, reciever_index, 0, MPI_COMM_WORLD);
    }

    MPI_Send(data.b, data.row_count, MPI_DOUBLE, reciever_index, 0, MPI_COMM_WORLD);
}

ProcessInitData recieve_process_init_data(int sender_index, MPI_Status *status){
    //получение данных для инициализациии процесса
    ProcessInitData data;
    MPI_Recv(&(data.row_count), 1, MPI_INT, sender_index, 0, MPI_COMM_WORLD, status);
    MPI_Recv(&(data.col_count), 1, MPI_INT, sender_index, 0, MPI_COMM_WORLD, status);

    data.C = malloc(sizeof(double *) * data.row_count);
    data.b = malloc(sizeof(double) * data.row_count);

    int i;
    for(i = 0; i < data.row_count; i++){
        data.C[i] = malloc(sizeof(double) * (data.col_count));
        MPI_Recv(data.C[i], data.col_count, MPI_DOUBLE, sender_index, 0, MPI_COMM_WORLD, status);
    }

    MPI_Recv(data.b, data.row_count, MPI_DOUBLE, sender_index, 0, MPI_COMM_WORLD, status);
    return data;
}

void send_process_work_result(ProcessWorkResult data, int reciever_index){
    //отправление данных о работе сделанной процессом
    printf("Sending work result to %d\n", reciever_index); 
    MPI_Send(&(data.min_X), 1, MPI_DOUBLE, reciever_index, 0, MPI_COMM_WORLD);
    MPI_Send(&(data.row_no), 1, MPI_INT, reciever_index, 0, MPI_COMM_WORLD);
}

ProcessWorkResult recieve_process_work_result(int process_index, MPI_Status *status){
    //получение данных о работе сделанной другим процессом
    printf("Recieving process work data from %d\n", process_index);
    ProcessWorkResult data;
    MPI_Recv(&(data.min_X), 1, MPI_DOUBLE, process_index, 0, MPI_COMM_WORLD, status);
    MPI_Recv(&(data.row_no), 1, MPI_INT, process_index, 0, MPI_COMM_WORLD, status);
    return data;
}

ProcessWorkResult process_work(ProcessInitData pid, int column)
{
    printf("Finding min row and min X\n");
    ProcessWorkResult res;
    int j, is_first = 1;
    double X;
    
    for(j = 0; j < pid.row_count; j++){
        if(pid.C[j][column] != 0){
            X = pid.b[j] / pid.C[j][column];
        } else {
            X = DBL_MAX;
        }

        if(is_first){
            is_first = 0;
            res.min_X = X;
            res.row_no = j;            
        } else {
            if(X < res.min_X){
                res.min_X = X;
                res.row_no = j;
            }
        }
    }
    return res;
}

void print_vector(char* name, double *vector, int count){
    printf("Printing vector of '%s'\n", name);
    int i;
    for(i = 0; i < count; i++){
        printf("%4.2f ", vector[i]);
    }
    printf("\n");
}

void print_matrix(char* name, double **matrix, int row_number, int col_number){
    printf("Printing matrix '%s'\n", name);
    int i;
    for(i = 0; i < row_number; i++){
        print_vector(name, matrix[i], col_number);
    }
}

void free_vector(double *vector){
    free(vector);
}

void free_matrix(double **matrix, int row_number){
    int i;
    for(i = 0; i < row_number; i++){
        free_vector(matrix[i]);
    }

    free(matrix);
}

int main(int argc, char **argv)
{
    MPI_Init(NULL, NULL);

    // Find out rank, size
    int instance_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &instance_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    char* tracefile;
    tracefile = getenv("TVTRACE");
    if( tracefile != NULL ){
        printf( "tv tracefile=%s\n", tracefile);
        MPI_Pcontrol(TRACEFILES, NULL, tracefile, 0);
    }
    else{
        MPI_Pcontrol(TRACEFILES, NULL, "trace", 0);
    }

    MPI_Pcontrol(TRACELEVEL, 1,1,1);
    MPI_Pcontrol(TRACENODE, 1000000, 0, 1);     
    
    MPI_Status status;


start:
    printf("Instance rank: %d, world_size: %d. Start\n", instance_rank, world_size);
    
    // Время для замеров
    double start_time, end_time, delta;
    start_time = MPI_Wtime();

    ProcessInitData pid;
    int i, j, k;
    
    //Здесь код симплекс метода
    if(instance_rank == 0){
        SimplexTable st = init_random_simplex_table(3000,6000);

        int row_length = st.n + st.m;
        //чтобы получилось обработать и нечетно количество условий
        int max_rows_per_proc = (st.m + 1) / world_size;

        double *q;
        q = malloc(sizeof(double) * (row_length));
        for(i = 0; i < row_length; i++){
            q[i] = st.table[0][i];
        }

        #pragma omp parallel for private(row_length, max_rows_per_proc)
        for(i = 0; i < world_size; i++){
            ProcessInitData data;
            data = get_process_init_data(st, max_rows_per_proc, row_length, i);
            //рассылаем данные слейвам
            if(i > 0){
                send_process_init_data(data, i);
            } else {
                pid = data;
            }
        }

        //тут уже параллельность
        while(1){
            int column = -1;
            //смотрим только по количеству переменных
            for(i = 0; i < st.n; i++){
                if(q[i] < 0){
                    column = i;
                    break;
                }
            }

            //рассылаем всем column, он же признак завершения работы всего
            MPI_Bcast(&column, 1, MPI_INT, instance_rank, MPI_COMM_WORLD);
            if(column < 0){
                break;
            } 

            //работаем сами
            ProcessWorkResult work_res;
            work_res = process_work(pid, column);

            //получаем результаты, выбираем лучший
            int min_L = 0;
            ProcessWorkResult any_process_res, best_work_res = work_res;

            for(i=1; i < world_size; i++){
                any_process_res = recieve_process_work_result(i, &status);
                if(any_process_res.min_X < best_work_res.min_X){
                    best_work_res = any_process_res;
                    min_L = i;
                }
            }

            //сообщаем слевам, кто лучший
            MPI_Bcast(&min_L, 1, MPI_INT, instance_rank, MPI_COMM_WORLD);
            double *g;
            double h;
            if(min_L != 0){
                //лучший процесс - один из слейвов
                //ждем g и h
                g = malloc(sizeof(MPI_DOUBLE) * (pid.col_count));
                MPI_Bcast(g, pid.col_count, MPI_DOUBLE, min_L, MPI_COMM_WORLD);
                MPI_Bcast(&(h), 1, MPI_DOUBLE, min_L, MPI_COMM_WORLD);
            } else {
                //лучший процесс - мастер
                //отпрака всем g=С[row_no] и h=b[row_no]
                g = pid.C[best_work_res.row_no];
                h = pid.b[best_work_res.row_no];
                MPI_Bcast(g, pid.col_count, MPI_DOUBLE, instance_rank, MPI_COMM_WORLD);
                MPI_Bcast(&(h), 1, MPI_DOUBLE, instance_rank, MPI_COMM_WORLD);
            }

            //а здесь у нас уже есть g и h
            //изменяем C и b

            if (min_L == instance_rank){
                //особенности процесса с минимумом
                
                #pragma omp parallel for
                for(i = 0; i < row_length; i++){
                    pid.C[best_work_res.row_no][i] = g[i]/g[column];
                }
                pid.b[best_work_res.row_no] = h / g[column];
            }

            //особенности мастера
            double a;
            if(g[column] != 0){
                a = q[column] / g[column];
            } else {
                a = 0;
            }

            #pragma omp parallel for
            for(i = 0; i < row_length; i++){
                q[i] = q[i] - (a * g[i]);
            }
            st.z[0] = st.z[0] - a * h;
        }
        
        #pragma omp parallel for
        for(i = 0; i < row_length; i++){
            st.table[0][i] = q[i];
        }

        for(k = 0; k < world_size; k++){
            ProcessInitData data;
            int max_rows_per_proc;
            if(!k){
                data = pid;
                max_rows_per_proc = pid.row_count;
            } else {
                data = recieve_process_init_data(k, &status);
            }
           
            #pragma omp parallel for private(max_rows_per_proc, k)
            for(i = 0; i < data.row_count; i++){
                int row_index = k * max_rows_per_proc + i + 1;
                for(j = 0; j < data.col_count; j++){
                    st.table[row_index][j] = data.C[i][j];
                }
                st.z[row_index] = data.b[i];
            }

            if(k != 0)
            {
                free_matrix(data.C, data.row_count);
                free_vector(data.b);
            }
        }

//        print_matrix("table", st.table, st.m + 1, st.m + st.n);
//        print_vector("table z", st.z, st.m + 1);
//        print_vector("table b", st.b, st.m + 1);

        //чистим данные!!!
        free_matrix(st.table, st.m + 1);
        free_vector(st.z);
        free_vector(st.b);
    } else {
        //получаем информацию для инициализации
        pid = recieve_process_init_data(0, &status);
        while(1){
            int column;
            //получаем column или финалайз
            MPI_Bcast(&(column), 1, MPI_INT, 0, MPI_COMM_WORLD);
            printf("Instance %d, column %d\n", instance_rank, column);
            if(column < 0){
                break;
            }

            ProcessWorkResult work_res;
            work_res = process_work(pid, column);
            //шлем результат мастеру
            send_process_work_result(work_res, 0);

            //получаем min_L
            int min_L;
            MPI_Bcast(&(min_L), 1, MPI_INT, 0, MPI_COMM_WORLD);

            double *g;
            double h;
            if(min_L != instance_rank){
                //лучший процесс не мы
                //ждем g и h
                g = malloc(sizeof(double) * (pid.col_count));
                MPI_Bcast(g, pid.col_count, MPI_DOUBLE, min_L, MPI_COMM_WORLD);
                MPI_Bcast(&(h), 1, MPI_DOUBLE, min_L, MPI_COMM_WORLD);
            } else {
                //лучший процесс - мы
                //отпрака всем g=С[row_no] и h=b[row_no]
                g = pid.C[work_res.row_no];
                h = pid.b[work_res.row_no];
                MPI_Bcast(g, pid.col_count, MPI_DOUBLE, instance_rank, MPI_COMM_WORLD);
                MPI_Bcast(&(h), 1, MPI_DOUBLE, instance_rank, MPI_COMM_WORLD);
            }

            //а здесь у нас уже есть g и h
            //изменяем C и b
            //особенность процесс с минимумом
            if(min_L == instance_rank){
                #pragma omp parallel for
                for(i = 0; i < pid.col_count; i++){
                    pid.C[work_res.row_no][i] = g[i]/g[column];
                }
                pid.b[work_res.row_no] = h / g[column];
            }

            if(min_L != instance_rank){
                free_vector(g);
            }
        }

        //отсылаем C и b
        send_process_init_data(pid, 0);
        //чистим данные
        free_matrix(pid.C, pid.row_count);
        free_vector(pid.b);
    }

finalize:
    end_time = MPI_Wtime();
    delta = end_time - start_time;
    printf("Instance rank: %d, finished for %lf seconds\n", instance_rank, delta);
    
    //finalizing
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}
