

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
    
    int comm_sz, my_rank;

    double start, finish, local_execution_time, execution_time;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

 
    double xmax;
    double xmin;
    double ymax;
    double ymin;

    int maxiter;


    int xres;
    int yres;

  
    double dx;
    double dy;

    char* filename;

    char* result = NULL;

    if (my_rank == 0)
    {
        if (argc != 8)
        {
            printf("Usage:   %s <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>\n", argv[0]);
            printf("Example: %s 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm\n", argv[0]);
            exit(EXIT_FAILURE);
        }

    
        xmin = atof(argv[1]);
        xmax = atof(argv[2]);
        ymin = atof(argv[3]);
        ymax = atof(argv[4]);

      
        maxiter = atoi(argv[5]);

    
        xres = atoi(argv[6]);
        yres = (xres * (ymax - ymin)) / (xmax - xmin);

       
        filename = argv[7];

        
        dx = (xmax - xmin) / xres;
        dy = (ymax - ymin) / yres;

        
        result = (char*)malloc(6 * yres * xres * sizeof(char));
    }

    start = MPI_Wtime();

    MPI_Bcast(&xres, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&yres, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxiter, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ymax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int chunk_size = yres;
    int total_chunk = 0;

    char* local_result = (char*)malloc(6 * yres * xres * sizeof(char));

   
    for (int i = 0; i < (6 * yres * xres); i++)
        local_result[i] = 'a';

    int local_index = 0;


    if (my_rank == 0)
    {
      
        while (total_chunk < yres)
        {
           
            int stop_signal = 0;

            for (int i = 1; i < comm_sz; i++)
                MPI_Send(&stop_signal, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

         
            int signal_rank[2];
            signal_rank[0] = 0;
            signal_rank[1] = 0;

            int temp_rank[2];
            temp_rank[0] = 0;
            temp_rank[1] = 0;

            int no_work_for_you = 0;

          
            while (signal_rank[0] != 1)
                MPI_Recv(signal_rank, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int start_index = total_chunk;
            total_chunk += chunk_size;

            MPI_Send(&no_work_for_you, 1, MPI_INT, signal_rank[1], 0, MPI_COMM_WORLD);
            MPI_Send(&chunk_size, 1, MPI_INT, signal_rank[1], 2, MPI_COMM_WORLD);
            MPI_Send(&start_index, 1, MPI_INT, signal_rank[1], 3, MPI_COMM_WORLD);

       
            for (int rank = 1; rank < comm_sz; rank++)
            {
               
                no_work_for_you = 1;
                if (rank != signal_rank[1])
                {
                    MPI_Recv(temp_rank, 2, MPI_INT, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&no_work_for_you, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
                }
            }         
        }
    
        int stop_signal = 1;
        for (int i = 1; i < comm_sz; i++)
            MPI_Send(&stop_signal, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
    else	
    {
        
        int stop_signal = 0;
        while (stop_signal != 1)	
        {
            int chunk_size = 0;
            MPI_Recv(&stop_signal, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          
            if (stop_signal == 0)
            {
                int signal_rank[2];
                signal_rank[0] = 1;
                signal_rank[1] = my_rank;

                MPI_Send(signal_rank, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);

                int no_work_for_you;
                MPI_Recv(&no_work_for_you, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

             
                if (no_work_for_you != 1)
                {
                   
                    int start_index;
                    MPI_Recv(&chunk_size, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&start_index, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

             
                   
                    local_index = start_index;

                    double x, y; 
                    int k; 
                    for (int j = start_index; j < start_index + chunk_size; j++) {
                        y = ymax - j * dy;
                        for (int i = 0; i < xres; i++) {
                            double u = 0.0;
                            double v = 0.0;
                            double u2 = u * u;
                            double v2 = v * v;
                            x = xmin + i * dx;
                            
                            for (k = 1; k < maxiter && (u2 + v2 < 4.0); k++) {
                                v = 2 * u * v + y;
                                u = u2 - v2 + x;
                                u2 = u * u;
                                v2 = v * v;
                            };
                          
                            if (k >= maxiter) {
                           
                                for (int a = 0; a < 6; a++)
                                    local_result[local_index++] = 0;
                            }
                            else {
                               
                                local_result[local_index] = k >> 8;
                                local_result[local_index + 1] = k & 255;
                                local_result[local_index + 2] = k >> 8;
                                local_result[local_index + 3] = k & 255;
                                local_result[local_index + 4] = k >> 8;
                                local_result[local_index + 5] = k & 255;
                                local_index += 6;
                            };
                        }
                    }                    
                }
            }
        }
    }
    

    finish = MPI_Wtime();

    if (my_rank != 0)
    {
        MPI_Send(local_result, 6 * yres * xres, MPI_CHAR, 0, my_rank, MPI_COMM_WORLD);
        free(local_result);
    }
    else
    {
        for (int rank = 1; rank < comm_sz; rank++)
        {
            char* temp = (char*)malloc(6 * yres * xres * sizeof(char));
            MPI_Recv(temp, 6 * yres * xres, MPI_CHAR, rank, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < 6 * yres * xres; i++)
            {
                if (temp[i] != 'a')
                    result[i] = temp[i];
            }
            free(temp);
        }  
    }
    
 
    local_execution_time = finish - start;

   
    MPI_Reduce(&local_execution_time, &execution_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
       
        FILE* fp = fopen(filename, "wb");

      
        fprintf(fp,
            "P6\n# Mandelbrot, xmin=%lf, xmax=%lf, ymin=%lf, ymax=%lf, maxiter=%d\n%d\n%d\n%d\n",
            xmin, xmax, ymin, ymax, maxiter, xres, yres, (maxiter < 256 ? 256 : maxiter));

        
        fwrite(result, 6 * yres * xres, 1, fp);

       
        printf("\nExecution time:  %.4f sec\n\n", execution_time / comm_sz);

        fclose(fp);
    }
   
    MPI_Finalize();

    return 0;
}