import java.util.*;

public class TestMatrix3 {

public static void main(String[] argument) {
   double[][] a = new double[4][4];
   double[][] b = new double[4][3];
   int errCodeF, errCodeP, errCodeD, errCodeSolve;
   double norm=0.0;
   Random r = new Random();
   int rValue;
   Matrix A, X, B, ASaved, BSaved, Res, Err;
   
   for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
         rValue = r.nextInt(20);
         a[i][j] = (double) rValue;
      }
   }
   
   for (int i=0; i<4; i++) {
      for (int j=0; j<3; j++) {
         rValue = r.nextInt(20);
         b[i][j] = (double) rValue;
      }
   }
      
   A = new Matrix(4,4,a);
   X = new Matrix(4,3);
   B = new Matrix(4,3,b);
   Res = new Matrix(4,3);
   Err = new Matrix(4,3);
   ASaved = new Matrix(A);
   BSaved = new Matrix(B);
   
   errCodeF = Matrix.factorlu(A);
   if (errCodeF == 2) System.out.println("System dimensions are inconsistant.");
   if (errCodeF == 1) System.out.println("System is singular.");
   if (errCodeF == 0) {
      errCodeSolve = Matrix.solvelu(A,B);
      if (errCodeSolve == 1) System.out.println("Dimensions inconsistant for some reason.");
      if (errCodeSolve == 0) {
         errCodeP = Matrix.prod(ASaved,B,Res);
         if (errCodeP == 1) System.out.println("Dimensions inconsistant for some reason.");
         if (errCodeP == 0) {
            errCodeD = Matrix.diff(BSaved,Res,Err);
            if (errCodeD == 1) System.out.println("Dimensions inconsistant for some reason.");
            if (errCodeD == 0) {
               for (int i=0; i<4; i++) {
                  for (int j=0; j<3; j++) {
                     norm = norm + (Err.M[i][j])*(Err.M[i][j]);
                  }
               }
            }
         }
      }
   }
   System.out.println("A:");
   for (int i=0;i<4;i++){
      System.out.format("%13.6f %13.6f %13.6f %13.6f%n", ASaved.M[i][0], ASaved.M[i][1],
              ASaved.M[i][2], ASaved.M[i][3]);
   }

   System.out.println("");
   System.out.println("");
   System.out.println("BSaved:");
   for (int i=0;i<4;i++){
      System.out.format("%13.6f %13.6f %13.6f%n", BSaved.M[i][0], BSaved.M[i][1],
              BSaved.M[i][2]);
   }

   System.out.println("");
   System.out.println("");
   System.out.println("Res:");
   for (int i=0;i<4;i++){
      System.out.format("%13.6f %13.6f %13.6f%n", Res.M[i][0], Res.M[i][1],
              Res.M[i][2]);
   }

   System.out.println("");
   System.out.println("");
   System.out.println("X:");
   for (int i=0;i<4;i++){
      System.out.format("%13.6f %13.6f %13.6f%n", B.M[i][0], B.M[i][1],
              B.M[i][2]);
   }


   System.out.println("The square of the norm of the error is " + norm + ".");

}

}
