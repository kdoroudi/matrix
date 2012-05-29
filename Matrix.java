//import java.lang.*;
public class Matrix {

/************************************************************************************************
A  Simple Matrix Class

Public variables:
   M:                        A two dimensional array the hold the Matrix entries
   numRows:                  An integer to hold the number of rows in the Matrix
   numCols:                  An integer to hold the number of columns in the Matrix
   
   
Public methods (other than constructors) in brief:
   boolean isSquare()            A boolean function that returns true if the Matrix is square, else false
   int sum(A,B,C)                A static integer function that adds Matrices A and B and places the value in C
   int diff(A,B,C)               A static integer function that subtracts Matrices A and B and places the
                                    value in C
   int prod(A,B,C)               A static integer function that multiplies Matrices A and B and places the
                                    value in C
   int solve(A,X,B)              A static integer function that solves AX=B by Gaussian Elimination
                                    Returns 0 if successful; returns 1 if A is singular; returns 2 if dimensions conflict
   int factorlu(A)               A static integer function that factors A without pivoting
                                    Returns 0 if successful; returns 1 if A cannot be factored; returns 2 if A is not square
   int factorlu(A,ipivot,iflag)  A static integer function that factors A with scaled partial pivoting
                                    Returns 0 if successful; returns 1 if A is singular; returns 2 if A is not square
                                    ipivot is an integer array that contains the pivoting information
                                    det(A) = iflag times the product of the diagonals
   int solvelu(A,B)              A static integer function that solves AX=B assuming that A has already been factored by factorlu
                                    Returns 0 if successful, returns 1 if A is singular; returns 2 if dimensions conflict
                                    B contains the right side on passage into solvelu, and X on the exit
   int solvelu(A,B,ipivot)       A static integer function that solves AX=B assuming that A has already been factored by factorlu
                                    Returns 0 if successful, returns 1 if A is singular; returns 2 if dimensions conflict
                                    B contains the righ side on passage into solvelu, and X on the exit
                                    ipivot contains the pivot information from factorlu
   int splitlu(A,L,U)            A static integer function that extracts the lower and upper triangular matrices from a previously
                                    factored matrix
                                    A contains both the lower and upper triangular matrices; L is the lower triangular matrix; 
                                    U the upper triangular matrix

Public methods and constructors in detail:
   Matrix (int nRows, int nCols)
   
   Constructs a Matrix with nRows rows and nCols columns
   The resulting Matrix is initialized to the zero Matrix
   Error message is printed and the program stops if nRows and nCols aren't both bigger than 1

---------------------------------------------------------------------------
   Matrix (int size)
   
   Constructs a square Matrix with size rows and size columns
   The resulting Matrix is initialized to the zero Matrix
   Error message is printed and the program stops if size isn't bigger than 1
   
---------------------------------------------------------------------------
   Matrix (int nRows, int nCols, double[][] A)
   
   Constructs a Matrix with nRows rows and nCols columns
   The resulting Matrix is initialized with the values from the two-dimensional array A
   Error message is printed and the program stops if nRows and nCols aren't both bigger than 1
   No indices on the array A are checked, so the program will stop with an array out of bounds
       exception if indices are incorrect
       
----------------------------------------------------------------------------
   Matrix (int size, A)
   
   Constructs a square Matrix with size rows and size columns
   The resulting Matrix is initialized with the values from the two-dimensional array A
   Error message is printed and the program stops if size isn't bigger than 1
   No indices on the array A are checked, so the program will stop with an array out of bounds
        exception if indices are incorrect
        
----------------------------------------------------------------------------
   Matrix (Matrix A)
   
   Constructs a Matrix with the same dimensions and values as Matrix A
   
----------------------------------------------------------------------------
   boolean isSquare()
   
   Returns true if the Matrix is square; returns false if the Matrix is not square
   
----------------------------------------------------------------------------
   int sum(Matrix A, Matrix B, Matrix C)
   
   A static integer function that adds Matrix A to Matrix B and places the result in Matrix C
   Returns 1 without changing Matrix C if the Matrix dimensions prevent the action
   Returns 0 if the operation is successful
   
--------------------------------------------------------------------------------
   int diff(Matrix A, Matrix B, Matrix C)
   
   A static integer function that subtracts Matrix B from Matrix A and places the result in Matrix C
   Returns 1 without changing Matrix C if the Matrix dimensions prevent the action
   Returns 0 if the operation is successful
   
---------------------------------------------------------------------------------
   int prod(Matrix A, Matrix B, Matrix C)
   
   A static integer function that multiplies Matrix A on the left with Matrix B on the right
   The resulting product is stored in Matrix C
   Returns 1 without changing Matrix C if the Matrix dimensions prevent the action
   Returns 0 if the operation is successful
   
----------------------------------------------------------------------------------
   int solve(Matrix A, Matrix X, Matric B)
   
   A static integer function that solves AX = B by Gaussian Elimination
   Returns 0 if successful; Returns 1 if A is singular; Returns 2 if dimensions conflict

----------------------------------------------------------------------------------
   int factorlu(Matrix A)           

   A static integer function that factors A without pivoting
   Returns 0 if successful; returns 1 if A cannot be factored; returns 2 if A is not square

-----------------------------------------------------------------------------------
   int factorlu(Matrix A, int[] ipivot, int iflag)  

   A static integer function that factors A with scaled partial pivoting
   Returns 0 if successful; returns 1 if A is singular; returns 2 if A is not square
   ipivot is an integer array that contains the pivoting information
   det(A) = iflag times the product of the diagonals

-----------------------------------------------------------------------------------
   int solvelu(Matrix A, Matrix B)              

   A static integer function that solves AX=B assuming that A has already been factored by factorlu
   Returns 0 if successful, returns 1 if A is singular; returns 2 if dimensions conflict
   B contains the right side on passage into solvelu, and X on the exit
   Assumes that no pivoting was done

------------------------------------------------------------------------------------
   int solvelu(Matrix A, Matrix B, int[] ipivot)       

   A static integer function that solves AX=B assuming that A has already been factored by factorlu
   Returns 0 if successful, returns 1 if A is singular; returns 2 if dimensions conflict
   B contains the righ side on passage into solvelu, and X on the exit
   ipivot contains the pivot information from factorlu

------------------------------------------------------------------------------------
   int splitlu(A,L,U)            
   A static integer function that extracts the lower and upper triangular matrices from a previously factored matrix
   A contains both the lower and upper triangular matrices; L is the lower triangular matrix; 
   U the upper triangular matrix

*************************************************************************************************/
   public double[][] M;
   public int numRows;
   public int numCols;
   
   public Matrix(int nRows, int nCols) {
      if ( (nRows < 1) || (nCols < 1) ) {
         System.out.println("The number of rows and columns must be positive.");
         System.exit(-1);
      }
      numRows = nRows;
      numCols = nCols;
      M = new double[numRows][numCols];
      for (int i=0;i<numRows;i++) {
         for (int j=0;j<numCols;j++) {
            M[i][j] = 0.0;
         }
      }
   }
   
   public Matrix(int size) {
      if (size < 1) {
         System.out.println("The number of rows and columns must be positive.");
         System.exit(-1);
      }
      numRows = size;
      numCols = size;
      M = new double[numRows][numCols];
      for (int i=0;i<numRows;i++) {
         for (int j=0;j<numCols;j++) {
            M[i][j] = 0.0;
         }
      }
   }
   
   public Matrix(int nRows, int nCols, double[][] A) {
      if ( (nRows < 1) || (nCols < 1) ) {
         System.out.println("The number of rows and columns must be positive.");
         System.exit(-1);
      }
      numRows = nRows;
      numCols = nCols;
      M = new double[numRows][numCols];
      for (int i=0;i<numRows;i++) {
         for (int j=0;j<numCols;j++) {
            M[i][j] = A[i][j];
         }
      }
   }

   public Matrix(int size, double[][] A) {
      if (size < 1) {
         System.out.println("The number of rows and columns must be positive.");
         System.exit(-1);
      }
      numRows = size;
      numCols = size;
      M = new double[numRows][numCols];
      for (int i=0;i<numRows;i++) {
         for (int j=0;j<numCols;j++) {
            M[i][j] = A[i][j];
         }
      }
   }


   public Matrix(Matrix A) {
      numRows = A.numRows;
      numCols = A.numCols;
      M = new double[numRows][numCols];
      for (int i=0;i<numRows;i++) {
         for (int j=0;j<numCols;j++) {
            M[i][j] = A.M[i][j];
         }
      }
   }
   
   public boolean isSquare() {
      if (numRows == numCols) return true;
      else return false;
   }
         
   public static int sum(Matrix A, Matrix B, Matrix C) {
      if ( (A.numRows != B.numRows ) || (A.numRows != C.numRows ) ||
           (A.numCols != B.numCols ) || (A.numCols != C.numCols ) ) {
         return 1;
      }
      else {
         for (int i = 0; i<A.numRows; i++) {
            for (int j = 0; j<A.numCols; j++) {
               C.M[i][j] =  A.M[i][j] + B.M[i][j];
            }
         }
         return 0;
      }
   }
   
   public static int diff(Matrix A, Matrix B, Matrix C) {
      if ( (A.numRows != B.numRows ) || (A.numRows != C.numRows ) ||
           (A.numCols != B.numCols ) || (A.numCols != C.numCols ) ) {
         return 1;
      }
      else {
         for (int i = 0; i<A.numRows; i++) {
            for (int j = 0; j<A.numCols; j++) {
               C.M[i][j] =  A.M[i][j] - B.M[i][j];
            }
         }
         return 0;
      }
   }
   
   public static int prod(Matrix A, Matrix B, Matrix C) {
      if ( (A.numRows != C.numRows ) || (B.numCols != C.numCols ) ||
           (A.numCols != B.numRows ) ) {
         return 1;
      }
      else {
         for (int i=0; i<A.numRows; i++ ) {
            for (int j=0; j<B.numCols; j++) {
               double value = 0;
               for (int k=0; k<A.numCols; k++) {
                  value = value + A.M[i][k]*B.M[k][j];
               }
               C.M[i][j] = value;
            }
         }
         return 0;
      }
   }
   
   public static int solve(Matrix A, Matrix X, Matrix B) {

// If A is not square, stop and give incorrect dimension indication.
      if (A.isSquare() == false) return 2;

// If A, X, and B have different numbers of rows or if X and B have different numbers
// of columns, stop and give incorrect dimenstion indication.      
      if ( (X.numRows != B.numRows ) || 
           (X.numCols != B.numCols ) ||
           (A.numRows != B.numRows ) ) return 2;

// If A is a 1x1 matrix and its value is zero, the matrix is singular.           
      if ( (A.numRows == 1) && (A.M[0][0] == 0.0) ) return 1;

// If A is a 1x1 matrix and its value is not zero, solution is simple.      
      if ( (A.numRows == 1) && (A.M[0][0] != 0.0) ) {
         for (int i=0;i<B.numCols;i++) {
            X.M[0][i] = B.M[0][i]/A.M[0][0];
         }
         return 0;
      }

// Now we consider the general cases.  Here i indicates the row currently being used
// to eliminate the values below.  We begin by swapping the rows to get the value largest
// in absolute value in position i,i.
      for (int i=0; i<A.numCols-1; i++) {
         double biggie = A.M[i][i];
         int biggiel = i;
         for (int l=i+1; l<A.numCols; l++) {
            if (Math.abs(A.M[l][i]) > Math.abs(A.M[biggiel][i])) {
               biggiel = l;
               biggie = A.M[l][i];
            }
         }

// If all values at i,i and below are zero, the system is singular.  Return 1 and stop.
// Otherwise swap rows.
         if (biggie == 0.0) return 1;
         if (biggiel > i) {
            for (int k=i; k<A.numCols; k++) {
               double tmp = A.M[i][k];
               A.M[i][k] = A.M[biggiel][k];
               A.M[biggiel][k] = tmp;
            }
            for (int k=0; k<B.numCols; k++) {
               double tmp = B.M[i][k];
               B.M[i][k] = B.M[biggiel][k];
               B.M[biggiel][k] = tmp;
            }

         }

// Now we must divide the ith rows through by A_{i,i}.  We divide the right side first
// then the left side, from right to left, so that A_{i,i} is divided last.  Otherwise,
// we would have to save A_{i,i} in a separate variable before the division started.
         for (int j=0; j<B.numCols; j++) {
            B.M[i][j] = B.M[i][j] / A.M[i][i];
         }

         for (int j=A.numCols-1; j>=i; j--) {
            A.M[i][j] = A.M[i][j]/A.M[i][i];
         }

// Now we need get zeros below A_{i,i}.  However, I don't really need to set A_{k,i}, for
// k > i, to zero.  Indeed, if I did, I would need to save the value of A_{k,i} to a separate
// variable before finishing the last step.

// Below, j denote the row that is being changed through the use of row i.  The variable k is used
// to denote the specific column that is being changed.  We begin at k=j+1 for the reasons in the
// paragraph above.
         for (int j=i+1; j<A.numCols; j++) {
            for (int k=i+1; k<A.numCols; k++) {
               A.M[j][k] = A.M[j][k] - A.M[j][i]*A.M[i][k];
            }
            for (int k=0; k<B.numCols; k++) {
               B.M[j][k] = B.M[j][k] - A.M[j][i]*B.M[i][k];
            }
         }
      }

// If A_{n,n} is zero, the matrix is singular.  Return 1 and stop.
      if (A.M[A.numCols-1][A.numCols-1] == 0.0) return 1;

// If A_{n,n} is not zero, first get x_n.
      for (int i=0; i<X.numCols; i++) {
         X.M[A.numCols-1][i] = B.M[A.numCols-1][i] / A.M[A.numCols-1][A.numCols-1];
      }

// Now backsolve to get x.
      for (int i=A.numCols-2; i>=0; i--) {
         for (int l=0; l<B.numCols; l++) {
            double tmp = 0.0;
            for (int k=i+1; k<A.numCols; k++) {
               tmp = tmp + A.M[i][k]*X.M[k][l];
            }
            X.M[i][l] = B.M[i][l] - tmp;
         }
      }
            
      return 0;
   }

   public static int factorlu(Matrix A) {
      if (A.isSquare() == false) return 2;
      if (A.numCols == 1) return 0;
      for (int k=1; k<A.numCols; k++) {
         if (A.M[0][0] == 0.0) return 1;
         A.M[k][0] = A.M[k][0]/A.M[0][0];
         if (k > 1) {
            for (int i=1; i<k; i++) {
               for (int j=0; j<i; j++) {
                  A.M[k][i] = A.M[k][i] - A.M[k][j]*A.M[j][i];
               }
               if (A.M[i][i] == 0.0) return 1;
               A.M[k][i] = A.M[k][i]/A.M[i][i];
            }
         }
         for (int i=k;i<A.numCols;i++) {
            for (int j=0; j<k; j++) {
               A.M[k][i] = A.M[k][i] - A.M[k][j]*A.M[j][i];
            }
         }
      }

      return 0;
   }

   public static int factorlu(Matrix A, int[] ipivot, int iflag) {
      if (A.isSquare() == false) return 2;
      iflag = 1;
      if (A.numCols == 1) return 0;
      double d[];
      d = new double[A.numCols];
      for (int i=0; i<A.numCols; i++) {
         ipivot[i]=i;
         double rowmax = 0.;
         for (int j=0; j<A.numCols; j++) {
            rowmax = Math.max(rowmax,Math.abs(A.M[i][j]));
         }
         if (rowmax == 0.0) {
            iflag = 0;
            rowmax = 1;
         }
         d[i] = rowmax;
      }
      for (int k=0; k<A.numCols-1; k++) {
         double colmax = Math.abs(A.M[k][k])/d[k];
         int istar = k;
         for (int i=k+1; i<A.numCols; i++) {
            double awikod = Math.abs(A.M[i][k])/d[i];
            if (awikod > colmax) {
               colmax = awikod;
               istar = i;
            }
         }
         if (colmax == 0.0) {
            iflag = 0;
         } else {
            if (istar > k) {
               iflag = -1*iflag;
               int i = ipivot[istar];
               ipivot[istar]= ipivot[k];
               ipivot[k]=i;
               double temp = d[istar];
               d[istar]=d[k];
               d[k]=temp;
               for (int j=0; j<A.numCols; j++) {
                  temp = A.M[istar][j];
                  A.M[istar][j]=A.M[k][j];
                  A.M[k][j]=temp;
               }
            }
            for (int i=k+1; i<A.numCols; i++) {
               A.M[i][k] = A.M[i][k]/A.M[k][k];
               double ratio = A.M[i][k];
               for (int j=k+1; j<A.numCols; j++) {
                  A.M[i][j] = A.M[i][j] - ratio*A.M[k][j];
               }
            }
         }
      }
      if (iflag == 0) return 1;
      return 0;
   }

   public static int solvelu(Matrix A, Matrix B) {
      if (A.isSquare() == false) return 2;
      if (A.numRows != B.numRows) return 2;
      if (A.numRows == 1) {
         if (A.M[0][0] == 0.0) return 1;
         for (int i=0; i<B.numCols; i++) {
            B.M[0][i] = B.M[0][i]/A.M[0][0];
         }
         return 0;
      }
      for (int k=1; k<A.numCols; k++) {
         for (int i=0; i<B.numCols; i++) {
            for (int j=0; j<k; j++) {
               B.M[k][i] = B.M[k][i] - A.M[k][j]*B.M[j][i];
            }
         }
      }

      if (A.M[A.numCols-1][A.numCols-1] == 0) return 1;
      for (int i=0; i<B.numCols; i++) {
         B.M[A.numCols-1][i] = B.M[A.numCols-1][i]/A.M[A.numCols-1][A.numCols-1];
      }
      for (int k=A.numCols-2; k>= 0; k--) {
         if (A.M[k][k] == 0.0) return 1;
         for (int i=0; i<B.numCols; i++) {
            for (int j=k+1; j<A.numCols; j++) {
               B.M[k][i] = B.M[k][i] - A.M[k][j]*B.M[j][i];
            }
            B.M[k][i] = B.M[k][i]/A.M[k][k];
         }
      }

      return 0;
   }

   public static int solvelu(Matrix A, Matrix B, int[] ipivot) {
      if (A.isSquare() == false) return 2;
      if (A.numRows != B.numRows) return 2;
      if (A.numRows == 1) {
         if (A.M[0][0] == 0.0) return 1;
         for (int i=0; i<B.numCols; i++) {
            B.M[0][i] = B.M[0][i]/A.M[0][0];
         }
         return 0;
      }
      Matrix X;
      X = new Matrix(B.numRows,B.numCols);
      for (int col=0; col<X.numCols; col++) {
         int ip = ipivot[0];
         X.M[0][col] = B.M[ip][col];
         for (int i=1; i<A.numCols; i++) {
            double sum =0.0;
            for (int j=0; j<i; j++) {
               sum = A.M[i][j]*X.M[j][col] + sum;
            }
            ip = ipivot[i];
            X.M[i][col] = B.M[ip][col] - sum;
         }
         if (A.M[A.numCols-1][A.numCols-1] == 0) return 1;
         X.M[A.numCols-1][col] = X.M[A.numCols-1][col]/A.M[A.numCols-1][A.numCols-1];
         for (int i=A.numCols-2; i>=0; i--) {
            double sum = 0.0;
            for (int j=i+1; j<A.numCols; j++) {
               sum = A.M[i][j]*X.M[j][col] + sum;
            }
            X.M[i][col] = X.M[i][col] - sum;
            if (A.M[i][i] == 0.0) return 1;
            X.M[i][col] = X.M[i][col]/A.M[i][i];
         }
      }
      for (int i=0; i<X.numRows; i++) {
         for (int j=0; j<X.numCols; j++) {
            B.M[i][j] = X.M[i][j];
         }
      }   
      return 0;
   }



   public static int splitlu(Matrix A, Matrix L, Matrix U) {
      if (A.isSquare() == false) return 1;
      if (L.isSquare() == false) return 1;
      if (U.isSquare() == false) return 1;
      if (A.numCols != L.numCols) return 1;
      if (A.numCols != U.numCols) return 1;
      for (int i=0; i<A.numCols; i++) {
         for (int j=0;j<A.numCols; j++) {
            if (i>j) {
               L.M[i][j] = A.M[i][j];
               U.M[i][j] = 0.0;
            }
            if (i == j) {
               L.M[i][i] = 1.0;
               U.M[i][i] = A.M[i][i];
            }
            if (i<j) {
               L.M[i][j] = 0.0;
               U.M[i][j] = A.M[i][j];
            }
         }
      }
      return 0;
   }
}
