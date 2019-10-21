  #include<stdio.h>
#include<stdlib.h>

void multiply1DMatrices(double **matrix1, double *matrix2, double *result, int rowCount1, int colCount1);

void multiply1DMatrices(double **matrix1, double *matrix2, double *result, int rowCount1, int colCount1){
	double num = 0;
	int r = 0;
	int i;

	for(i=0; i<colCount1; i++){
		num = num + (matrix1[r][i] * matrix2[i]);
		if(i == colCount1 - 1){
			result[r] = num;
			num = 0;
			r++;
			i = -1;
		}
		if(r == rowCount1){	//rethink
			break;
		}
	}
	return;
}

void multiplyMatrices(double **matrix1, double **matrix2, double **result, int rowCount1, int colCount1, int colCount2);

void multiplyMatrices(double **matrix1, double **matrix2, double **result, int rowCount1, int colCount1, int colCount2){
	double num = 0;
	int r = 0;
	int c = 0;
	int i;

	for(i=0; i<colCount1; i++){
		num = num + (matrix1[r][i] * matrix2[i][c]);
		if(i == colCount1 - 1){
			result[r][c] = num;
			num = 0;
			c++;
			i = -1;
			if(c == colCount2){
				r++;
				c = 0;
			}
		}
		if(r == rowCount1){
			break;
		}
	}
	return;
}

void concatMatrices(double **matrix1, double **matrix2, double **concat, int rows, int cols);

void concatMatrices(double **matrix1, double **matrix2, double **concat, int rows, int cols){
	int halfCols = rows;
	int i;
	int j;
	int k = 0;
	
	for(i=0; i < rows; i++){
		for(j=0; j < halfCols; j++){
			concat[i][j] = matrix1[i][j];
		}
	}
	
	for(i=0; i < rows; i++){
		for(j=halfCols; j < cols; j++){
			concat[i][j] = matrix2[i][k];
			k++;
		}
		k = 0;
	}
	return;
}

void rowDivision(double** matrix, int row1, double constant, int sizeV);

void rowDivision(double** matrix, int row1, double constant, int sizeV){
	int col;
 	for(col = 0; col < (sizeV * 2); col++){
		matrix[row1][col] = matrix[row1][col] / constant;
	}

	return;
		
}

void rowSubMultiply(double** matrix, int row1, int row2, double constant, int sizeV);

void rowSubMultiply(double** matrix, int row1, int row2, double constant, int sizeV){
	int col;
	for(col=0; col < (sizeV * 2); col++){
		matrix[row1][col] = matrix[row1][col] - constant * matrix[row2][col];
	}
	return;
}



int main(int argc, char** argv){
	int N;	// # of examples in a training set
	int K;	// # of attributes in an example

	double **X;	// Matrix representing training set Matrix N x (K + 1)
	double **T;	// Matrix representing transpose of Matrix X
	double *Y;	// Matrix representing prices of the houses N x 1
	double *W;	// Matrix representing weights of the houses (K + 1) x 1

	FILE *fp = fopen(argv[1], "r");
	fscanf(fp, "%d\n", &K);
	fscanf(fp, "%d\n", &N);
	
	//training set Matrix created N (rows) x (K + 1) (cols)	
	int i;
	int j;
	X = (double**)malloc(sizeof(double*) * N);
	for(i=0; i < N; i++){
		X[i] = (double*)malloc(sizeof(double) * (K + 1));
	}

	//transpose Matrix created (K + 1) (rows) x N (cols)
	T = (double**)malloc(sizeof(double*) * (K + 1));
	for(i=0; i < (K + 1); i++){
		T[i] = (double*)malloc(sizeof(double) * N);
	}
	
	//price Matrix created N (rows) x 1 (cols)
	Y = (double*)malloc(sizeof(double) * N);

	//weight Matrix created (K + 1) (rows) x 1 (cols)	
	W = (double*)malloc(sizeof(double) * (K + 1));
	
	//Loads training set Matrix and price Matrix
	double input;
	int lastIndex = 0; //boolean to determine to skip a comma or new-line
	i = 0;
	j = 0;
	while(i < N){
		if(j == 0){
			X[i][j] = 1;
			j++;
			continue;
		}
		fscanf(fp, "%lf%*c" , &input);
		X[i][j] = input;
		if(lastIndex == 1){
			fscanf(fp, "%lf\n", &input);
			Y[i] = input;
			i++;
			j = 0;
			lastIndex = 0;
			continue;
		}
		if(j == (K - 1)){
			lastIndex = 1;
		}
		j++;
	}

	//Loads transpose Matrix of X
	for(i=0; i < N; i++){	//possible here
		for(j=0; j < (K + 1); j++){
			T[j][i] = X[i][j];
		}
	}

	//Creates empty Matrix of X(Matrix) * T(Matrix)-creating K+1(rows) x K+1(cols)	
	double **XT;	//Multiplication of X * T Matrix
	XT = (double**)malloc(sizeof(double*) * (K+1));
	for(i=0; i < (K+1); i++){
		XT[i] = (double*)malloc(sizeof(double) * (K+1));
	}
	//XT is now loaded after multiplyMatrices function runs
	multiplyMatrices(T, X, XT, (K+1), N, (K+1));		//possible overflow here

	//Creates and loads Identity Matrix (K+1)(rows) x (K+1)(cols)
	double **I;
	I = (double**)malloc(sizeof(double*) * (K+1));
	for(i=0; i < (K+1); i++){
		I[i] = (double*)malloc(sizeof(double) * (K+1));
	}

	for(i=0; i < (K+1); i++){
		for(j=0; j < (K+1); j++){
			if(i == j){
				I[i][j] = 1;
			}else{
				I[i][j] = 0;
			}	
		}
	}
	
	//Creates empty Matrix of XT concatenated with I -creating N(rows) x 2*N(cols)
	double **augXTI;
	augXTI = (double**)malloc(sizeof(double*) * (K+1));
	for(i=0; i < (K+1); i++){
		augXTI[i] = (double*)malloc(sizeof(double) * (2*(K+1)));
	}

	//augXTI is now loaded after concatMatrices function runs
	concatMatrices(XT, I, augXTI, (K+1), (2*(K+1)));	//possible overflow here

	for(i=0; i < (K+1); i++){
		for(j=0; j < ((K+1)*2); j++){
			if(i == j){
				if(augXTI[i][j] != 1){
					rowDivision(augXTI, i, augXTI[i][j], (K+1));
					break;		//restart as it hits pivot point
				}
			}else{
				if(augXTI[i][j] != 0){
					rowSubMultiply(augXTI, i, j, augXTI[i][j], (K+1));
				}
			}
		}
	}
	
	int otherHalf = (K+1)-2;
	for(i = (K+1)-2; i > -1; i--){
		for(j = (K+1)-1; j != otherHalf; j--){
			if(augXTI != 0){
				rowSubMultiply(augXTI, i, j, augXTI[i][j], (K+1)); //possible
			}
		}
		otherHalf--;
	}	

	//matrix following W = (Inverse) * T (Part of the equation)
	double **IT;
	IT = (double**)malloc(sizeof(double*) * (K+1));
	for(i=0; i < (K+1); i++){
		IT[i] = (double*)malloc(sizeof(double) * N);
	}
	
	//Identity matrix converted to inverse matrix and multiplied with the transpose matrix
	for(i=0; i < (K+1); i++){
		for(j=0; j < (K+1); j++){
			I[i][j] = augXTI[i][j + (K+1)];
		}
	}

	multiplyMatrices(I, T, IT, (K+1), (K+1), N);
	multiply1DMatrices(IT, Y, W, (K+1), N);	

	//Creating testFile Matrix size V (rows) x K (cols)
	double **testFile;
	int V;		//Represents # of rows in testFile
	FILE* fp2 = fopen(argv[2], "r");
	fscanf(fp2, "%d\n", &V);
	testFile = (double**)malloc(sizeof(double*) * V);
	for(i=0; i < V; i++){
		testFile[i] = (double*)malloc(sizeof(double) * (K+1));
	}
	for(i=0; i < V; i++){
		for(j=0; j < K+1; j++){
			if(j == 0){
				testFile[i][j] = 1;
				continue;
			}
			if(j == K){
				fscanf(fp2, "%lf\n", &input);
				testFile[i][j] = input;
			}else{
				fscanf(fp2, "%lf%*c", &input);
				testFile[i][j] = input;
			}
		}
	}
	double *result;
	result = (double*)malloc(sizeof(double) * V);
	multiply1DMatrices(testFile, W, result, V, (K+1));

	for(i=0; i < V; i++){
		printf("%0.0f\n", result[i]);
	}	
	

	return 0;
}
