WeChat: cstutorcs
QQ: 749389476
Email: tutorcs@163.com
#include "params.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 	/* For OpenMP & MPI implementations, use their walltime! */

void init(double u[N1][N2][N3]){
  for (int n1=0; n1 < N1; n1++){
    for (int n2=0; n2 < N2; n2++){
      for (int n3=0; n3 < N3; n3++){
        u[n1][n2][n3] = u0(n1,n2,n3);
      }
    }
  }
};

void dudt(const double u[N1][N2][N3], double du[N1][N2][N3]) {
  double sum;
  int count;
  for (int n1 = 0; n1 < N1; n1++) {
    for (int n2 = 0; n2 < N2; n2++) {
      for (int n3 = 0; n3 < N3; n3++) {
        sum = 0.0;
        count = 0;
        for (int l1 = imax(0, n1 - ml); l1 <= imin(n1 + ml, N1-1); l1++) {
          for (int l2 = imax(0, n2 - ml); l2 <= imin(n2 + ml, N2-1); l2++) {
            for (int l3 = imax(0, n3 - ml); l3 <= imin(n3 + ml, N3-1); l3++) {
              sum += u[l1][l2][l3];         // Accumulate the local sum in sum
              count++;                      // Increment the count
            }
          }
        }
        du[n1][n2][n3] = (sum / count);     // store the local mean in du
      }
    }
  }
};

void step(double u[N1][N2][N3], const double du[N1][N2][N3]) {
  for (int n1 = 0; n1 < N1; n1++) { 
    for (int n2 = 0; n2 < N2; n2++) {
      for (int n3 = 0; n3 < N3; n3++) {
        u[n1][n2][n3] = r * du[n1][n2][n3] * (1.0 - du[n1][n2][n3]);
      }
    }
  }
};

void stat(double* stats, const double u[N1][N2][N3]) {
  double mean =    0.0;
  double uvar =    0.0;
  double umin =  100.0; 
  double umax = -100.0;
  for (int n1 = 0; n1 < N1; n1++) {
    for (int n2 = 0; n2 < N2; n2++) {
      for (int n3 = 0; n3 < N3; n3++) {
        mean += u[n1][n2][n3] / (N1 * N2 * N3);
        if (u[n1][n2][n3] > umax){
          umax = u[n1][n2][n3];
        }
        if (u[n1][n2][n3] < umin){
          umin = u[n1][n2][n3];
        }
      }
    }
  }
  for (int n1 = 0; n1 < N1; n1++) {
    for (int n2 = 0; n2 < N2; n2++) {
      for (int n3 = 0; n3 < N3; n3++) {
        uvar += (u[n1][n2][n3] - mean) * (u[n1][n2][n3] - mean) / (N1 * N2 * N3);
      }
    }
  }
  stats[0] = mean;
  stats[1] = umin;
  stats[2] = umax;
  stats[3] = uvar;
};

void write(const double u[N1][N2][N3], const int m) {
  char outfile[80];
  int fileSuccess = sprintf(outfile, "state_%i.txt", m);
  if (fileSuccess > 0) {
    FILE *fptr = fopen(outfile, "w");
    for (int n3 = 0; n3 < N3; n3++) {
      for (int n2 = 0; n2 < N2; n2++) {
        for (int n1 = 0; n1 < N1; n1++) {
          // this segfaults when fptr is null.
          fprintf(fptr, "%2.4f\t", u[n1][n2][n3]);
        }
        fprintf(fptr, "\n");
      }
      fprintf(fptr, "\n");
    }
  } else {
    printf("Failed to write state_%i.txt!\n", m);
  }
};


int main(int argc, char **argv){
  
  double  u[N1][N2][N3];
  double du[N1][N2][N3];
  double stats[M/mm][4];
  int writeInd = 0;
  
  init(u);   
    
  clock_t t0 = clock();                   // for timing serial code
  
  for (int m = 0; m<M; m++){
    dudt(u, du);
    step(u, du);
    if (m%mm == 0){
      writeInd = m/mm;
      stat(&stats[writeInd][0], u);     // Compute statistics and store in stat
    }
    /*
    write(u, m);                        // Slow diagnostic output!
    */
  }
  double t1 = (double)(clock() - t0) / (CLOCKS_PER_SEC);     // for timing serial code
  
  FILE *fptr = fopen("stats.dat","w");
  fprintf(fptr, "iter\t\tmean\t\tmin\t\tmax\t\tvar\n");  // write stats to file
  for (int m = 0; m<(M/mm); m++){
    fprintf(fptr, "%6.0f\t%02.5f\t%02.5f\t%02.5f\t%02.5f\n", 
        (double)(m*mm), stats[m][0], stats[m][1], stats[m][2], stats[m][3]);
  }
  fclose(fptr);
  
  double t2 = (double)(clock() - t0) / (CLOCKS_PER_SEC) - t1; // timing writes
  printf("(%3d,%3d,%3d): average iteration time per element:\t%02.16fs\n", 
    N1,N2,N3,t1 / (N1*N2*N3*M));
  printf("(%5d,%3d,%1d): average write time per element:\t\t%02.16fs\n", 
    M,mm,4,t2 / (4*M/mm));
  
  return 0;
};
