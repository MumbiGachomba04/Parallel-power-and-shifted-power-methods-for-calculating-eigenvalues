#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstring>
#include <mpi.h>
using namespace std;


// void matSub(double **firstMatrix, double **secondMatrix, double **resultMat,int size);
// void matvecMul(double **matrix, double *vec,double *result,int size);
// double power (double **matrix,vector<double>&eigenvec1,int size);
// double shiftedPower(double **matrixA,double **matrixB,double **shiftedMat,vector<double>&eigenvec1,vector<double>&eigenvec2,double eigenval1,int size);
// double twoNorm(double *arr,int size);
void matvecMul(double **matrix, double *vect,double *result,int rows,int size);
double twoNorm(double *arr,int size);

int main(int argc, char *argv[])
{
    double **MatrixA ; // main matrix
    double *eigenvec1;  
    double *resultantvec;
    double** ProcRows;
    double* ProcResult; 
    int size;
    double eigenval ;
    int num_proc, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int *row_cnt = new int (num_proc);

    if (rank == 0){
/****************Start reading file*********/
        string path = argv[1];
        ifstream input;
        input.open(path);
        if(!input)
        {
            cerr << "Failed to open File\n";
            exit(1);
        }
        input >> size ; 
        MatrixA = new double *[size];
        for (int k = 0; k < size; k++)
        {
            MatrixA[k] = new double [size];
        }
        for (int i = 0; i < size; ++i){ 
            for (int j = 0; j < size; ++j){
                    input >> MatrixA[i][j];
                
                }
        }
        input.close();
/****************End of reading file***********/
    for(int i=0;i<num_proc; i++){
		   row_cnt[i] = size / num_proc;	
		   if( i == num_proc-1){
			  row_cnt[i]  +=  (size  % (num_proc)) ; // remainder num of rows go to the last process
		   }
		}
    }
    ProcRows = new double *[row_cnt[rank]]; // block
    for (int i = 0 ; i < row_cnt[rank]; i++ ){
         ProcRows[i] = new double [size];
    }
    ProcResult = new double (row_cnt[rank]); 
    eigenvec1 = new double(size); 
    resultantvec = new double(size);
   

    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(row_cnt ,num_proc ,MPI_INT,0,MPI_COMM_WORLD);
    if (rank == 0){
        for (int i = 0; i < size; i++) {
            eigenvec1[i] = 1;
        }
    }
    for (int p = 0; p < 10; p++) {
            if (rank == 0){
                memcpy(ProcRows,MatrixA ,row_cnt[0] * size * sizeof(double));
                int cumilative_rows = row_cnt[0];
                for(int i=1; i< num_proc;i++){
                MPI_Send(&MatrixA[cumilative_rows * size], row_cnt[i] * size ,MPI_DOUBLE,i,123,MPI_COMM_WORLD);
                cumilative_rows+=row_cnt[i];
                }
            }
            else {
                MPI_Recv(ProcRows, row_cnt[rank] * size, MPI_DOUBLE,0,123,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            
            matvecMul(ProcRows, eigenvec1,ProcResult,row_cnt[rank],size);
            MPI_Barrier( MPI_COMM_WORLD );
            int *disp = new int (num_proc);
            int *recv_cnt = new int (num_proc);
            for(int i=0; i< num_proc;i++){
               recv_cnt[i] = row_cnt[i];
         	   disp[i+1] = disp[i] + row_cnt[i]  ;
         	}
            //MPI_Allgather(ProcResult,row_cnt[rank] , MPI_DOUBLE, resultantvec, row_cnt[rank],  MPI_DOUBLE, MPI_COMM_WORLD); 
        
            MPI_Allgatherv(ProcResult,row_cnt[rank], MPI_DOUBLE,resultantvec,recv_cnt,disp,MPI_DOUBLE,MPI_COMM_WORLD);
            if (rank == 0){
                double norm = twoNorm(resultantvec,size);
                for (int j = 0; j < size; j++) {
                    eigenvec1[j]=resultantvec[j]/norm;
                }
                for (int k = 0; k < size; k++){
                    eigenval += (eigenvec1[k] * resultantvec[k]);
                }
                cout << "eigenval::" << eigenval << endl;
            }
        // send updated eigenvec
            MPI_Bcast(eigenvec1 ,size ,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    MPI_Bcast(&eigenval ,1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (rank == 0){
        cout << "DOMINANT EIGENVALUE: " << eigenval  << endl;
     }

    MPI_Finalize(); 
    return 0;
}

// void matSub(double **firstMatrix, double **secondMatrix, double **resultMat,int size)
// {

//    for (int i = 0; i < size; i++){
//         for (int j = 0; j < size; j++){
//             resultMat[i][j] = firstMatrix[i][j] - secondMatrix[i][j];
//         }
//    }
// }

void matvecMul(double **matrix, double *vect,double *result,int rows,int size)
{
    int num_proc;
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    
    for (int i = 0; i < rows; i++)  {
            result[i] = 0;
        for (int j = 0; j < size; j++) {
            result[i] += ( matrix[i][j] * vect[j]);
        }
    }       
}

double twoNorm(double *arr,int size){
    
   int num_proc;
   MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
   double sum = 0.0;
   for(int i = 0; i < size; i++){
          sum += pow(arr[i],2);
   }  
   return sqrt(sum);
}

// double power (double **matrix,vector<double>&eigenvec1,int size){
//     double eigenval ;
//     for (int i = 0; i < 10; i++)  {
//          eigenval = 0;
//          vector<double> tempeigenvec1(size, 0);
//          matvecMul(matrix, eigenvec1,tempeigenvec1,size);
//          double norm = twoNorm(tempeigenvec1,size);
//          for (int j = 0; j < size; j++) {
//             eigenvec1[j]=tempeigenvec1[j]/norm;
//          }
//          for (int k = 0; k < size; k++){
//              eigenval += (eigenvec1[k] * tempeigenvec1[k]);
//          }
        
//     }
//     return eigenval;
// }

// double shiftedPower(double **matrixA,double **matrixB,double **shiftedMat,vector<double>&eigenvec1,vector<double>&eigenvec2,double eigenval1,int size){  
//       for (int m = 0; m < size; m++) {
//         for (int n = 0; n < size; n++) {
//             matrixB[m][n] = eigenval1 * (eigenvec1[m] * eigenvec1[n]);
//         }
//       }
//       matSub(matrixA,matrixB,shiftedMat,size);       
//       double val = power(shiftedMat,eigenvec2,size);
//       return val;
// }