
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


    int local_yres = yres / comm_sz;

    char* local_result = (char*)malloc(6 * local_yres * xres * sizeof(char));
    int local_index = 0;

    double x, y; 
    int k; 
    for (int j = 0 + (my_rank * local_yres); j < local_yres + (my_rank * local_yres); j++) {
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
   
    int send_count = 6 * local_yres * xres;
    int recv_count = 6 * local_yres * xres;
    MPI_Gather(local_result, send_count, MPI_CHAR, result, recv_count, MPI_CHAR, 0, MPI_COMM_WORLD);


    finish = MPI_Wtime();

    
    local_execution_time = finish - start;

  
    MPI_Reduce(&local_execution_time, &execution_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    /
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