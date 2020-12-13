package ru.kpekepsalt;

import java.util.Arrays;

public class Main {

    public static final double E = 1E-8;
    public static double[][] J;
    public static double[][] F;

    /*
    //Configuration for system with 2 variables
    public static final int xCount = 2;
    public static double[] h = new double[] {1d, 1d};
    public static double[] x0 = new double[] {2d, 1d};
    public static double[] x1 = new double[] {x0[0]+h[0], x0[1]+h[1]};
    */

    //Configuration for system with 3 variables
    public static final int xCount = 3;
    public static double[] h = new double[] {1d, 1d, 1d};
    public static double[] x0 = new double[] {1d, 2d, 1d};
    public static double[] x1 = new double[] {x0[0]+h[0], x0[1]+h[1], x0[2]+h[2]};


    //3x^2+2y^2-35=0
    //4x^2-3y^2-24=0
    //X = [+-3, +-2]
    //X0 = [1, 1]
    /*
    public static ParamFunctional[] SYSTEM = new ParamFunctional[]{
            (x) -> 3*x[0]*x[0]+2*x[1]*x[1]-35,
            (x) -> 4*x[0]*x[0]-3*x[1]*x[1]-24
    };
     */

    //2x^2-xy-2x+1=0
    //x+3ln(x)-y^2=0
    //X = [1.3734783534098085, -1.5249648363795216]
    //X0 = [2, 1]
    /*
    public static ParamFunctional[] SYSTEM = new ParamFunctional[]{
            (x) -> 2*x[0]*x[0]-x[0]*x[1]-5*x[0]+1,
            (x) -> x[0]+3*Math.log(x[0])-x[1]*x[1]
    };
     */

    //x+y+z-4=0
    //2xy-z^2-16=0
    //z+4=0
    //X = [4.0, 4.0, -4.0]
    //x0 = [1, 2, 1]
    public static ParamFunctional[] SYSTEM = new ParamFunctional[]{
            (x) -> x[0]+x[1]+x[2]-4,
            (x) -> 2*x[0]*x[1]-x[2]*x[2]-16,
            (x) -> x[2]+4
    };

    
    public static void main(String[] args) {
        //Initialize matrices
        F = new double[SYSTEM.length][1];
        J = new double[SYSTEM.length][xCount];

        //Calculate system with first approximation
        F(x0);
        //Calculate derivatives matrix
        J();

        //Inverse J matrix
        double[][] JINV = new double[J.length][J.length];
        double[][] copyJ = new double[J.length][J.length];
        for(int i=0;i<J.length;i++)
            copyJ[i] = Arrays.copyOf(J[i], J.length);

        try {
            inverse(copyJ, JINV);
            for(int i=0;i<J.length;i++)
                J[i] = Arrays.copyOf(JINV[i], J.length);
            root();
        } catch (Exception e) {
            System.out.println(e.getLocalizedMessage());
        }
    }


    /**
     * Calculates derivatives matrix as (f(x0)-f(x0+h))/h
     */
    public static void J() {
        for(int i=0;i<SYSTEM.length;i++) {
            for(int j=0;j<xCount;j++) {
                double[] xh = Arrays.copyOf(x0, x0.length);
                xh[j] += h[j];
                J[i][j] = Math.pow(
                        (SYSTEM[i].function(xh)-SYSTEM[i].function(x0))/h[j],
                        SYSTEM.length
                );
            }
        }
    }

    /**
     * Calculates equations of specified roots
     * @param x - roots array
     */
    public static void F(double[] x) {
        for(int i=0;i< SYSTEM.length;i++) {
            F[i][0] = SYSTEM[i].function(x);
        }
    }

    /**
     * Finding roots as xk = x(k-1)-F(x(k-1))*J^-1
     */
    public static void root() {
        int k = 0;
        double disc = E;
        while(disc>=E && k < 1000000) {
            double[][] mult;
            F(x0);
            mult = multiplication(J, F);
            double[] sub = sub(x0, mult);
            x0 = Arrays.copyOf(x1, x1.length);
            x1 = Arrays.copyOf(sub, sub.length);
            double[] discVector = discrepancy(x0, x1);
            disc = discrepancySum(discVector);
            k++;
        }
        System.out.println("ITERATIONS = "+k);
        System.out.println("X = "+Arrays.toString(x1));
        F(x1);
        Print(F, "SYSTEM WITH ROOT");
    }

    public static double discrepancySum(double[] discrepancy) {
        return Math.abs(Arrays.stream(discrepancy).sum());
    }

    public static double[] discrepancy(double[] x0, double[] x1) {
        double[] discrepancy = new double[x0.length];
        double[] F1 = new double[discrepancy.length];
        double[] F2 = new double[discrepancy.length];
        for(int i=0;i<SYSTEM.length;i++) {
            F1[i] = SYSTEM[i].function(x0);
            F2[i] = SYSTEM[i].function(x1);
        }
        discrepancy = sub(F2, F1);
        return discrepancy;
    }

    public static void Print(double[][] matrix, String title) {
        System.out.println("===="+title+"====");
        for(double[] arr : matrix)
            System.out.println(Arrays.toString(arr));
        System.out.println("==============");
    }

    /**
     * Inverse given matrix
     * @param squareMatrix Matrix for inverse
     * @param inverseMatrix Matrix for storing result
     * @throws Exception
     */
    public static void inverse(final double[][] squareMatrix,
                               final double[][] inverseMatrix)
            throws Exception {
        final int size = squareMatrix.length;
        if (squareMatrix[0].length != size || inverseMatrix.length != size
                || inverseMatrix[0].length != size) {
            throw new Exception(
                    "--- invalid length. column should be 2 times larger than row.");
        }
        for (int i = 0; i < size; ++i) {
            Arrays.fill(inverseMatrix[i], 0.0d);
            inverseMatrix[i][i] = 1.0d;
        }
        for (int i = 0; i < size; ++i) {
            findPivotAndSwapRow(i, squareMatrix, inverseMatrix, size);
            sweep(i, squareMatrix, inverseMatrix, size);
        }
    }

    private static void findPivotAndSwapRow(final int row,
                                            final double[][] squareMatrix0, final double[][] squareMatrix1,
                                            final int size) {
        int ip = row;
        double pivot = Math.abs(squareMatrix0[row][row]);
        for (int i = row + 1; i < size; ++i) {
            if (pivot < Math.abs(squareMatrix0[i][row])) {
                ip = i;
                pivot = Math.abs(squareMatrix0[i][row]);
            }
        }
        if (ip != row) {
            for (int j = 0; j < size; ++j) {
                final double temp0 = squareMatrix0[ip][j];
                squareMatrix0[ip][j] = squareMatrix0[row][j];
                squareMatrix0[row][j] = temp0;
                final double temp1 = squareMatrix1[ip][j];
                squareMatrix1[ip][j] = squareMatrix1[row][j];
                squareMatrix1[row][j] = temp1;
            }
        }
    }

    private static void sweep(final int row, final double[][] squareMatrix0,
                              final double[][] squareMatrix1, final int size)
            throws Exception {
        final double pivot = squareMatrix0[row][row];
        if (pivot == 0) {
            throw new Exception(
                    "Inverse failed. Invalid pivot +" +pivot);
        }
        for (int j = 0; j < size; ++j) {
            squareMatrix0[row][j] /= pivot;
            squareMatrix1[row][j] /= pivot;
        }
        for (int i = 0; i < size; i++) {
            final double sweepTargetValue = squareMatrix0[i][row];
            if (i != row) {
                for (int j = row; j < size; ++j) {
                    squareMatrix0[i][j] -= sweepTargetValue
                            * squareMatrix0[row][j];
                }
                for (int j = 0; j < size; ++j) {
                    squareMatrix1[i][j] -= sweepTargetValue
                            * squareMatrix1[row][j];
                }
            }
        }
    }

    /**
     * Find multiplication of matrices
     */
    public static double[][] multiplication(double[][] A, double[][] B) {

        int aRows = A.length;
        int aColumns = A[0].length;
        int bRows = B.length;
        int bColumns = B[0].length;

        if (aColumns != bRows) {
            throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
        }

        double[][] C = new double[aRows][bColumns];
        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bColumns; j++) {
                C[i][j] = 0.00000;
            }
        }

        for (int i = 0; i < aRows; i++) { // aRow
            for (int j = 0; j < bColumns; j++) { // bColumn
                for (int k = 0; k < aColumns; k++) { // aColumn
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }

    /**
     * Find matrices difference
     */
    public static double[] sub(double[] a, double[][] b) {
        double[] res = new double[a.length];
        for(int i=0;i<res.length;i++) {
            res[i] = a[i]-b[i][0];
        }
        return res;
    }

    /**
     * Find matrices difference
     */
    public static double[] sub(double[] a, double[] b) {
        double[] res = new double[a.length];
        for(int i=0;i<res.length;i++) {
            res[i] = a[i]-b[i];
        }
        return res;
    }

}
