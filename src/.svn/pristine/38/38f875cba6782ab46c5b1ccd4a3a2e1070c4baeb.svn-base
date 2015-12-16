/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.stat;

import cern.colt.list.DoubleArrayList;
import java.io.BufferedReader;
import java.io.FileReader;
import umontreal.iro.lecuyer.functionfit.BSpline;
import umontreal.iro.lecuyer.functionfit.LeastSquares;
import umontreal.iro.lecuyer.functionfit.SmoothingCubicSpline;

/**
 *
 * modified from http://www.codeproject.com/KB/recipes/LinReg.aspx by MX Li
 */
public class LinearRegression {

    double[][] V;            // Least squares and var/covar matrix
    public double[] C;      // Coefficients
    public double[] SEC;    // Std Error of coefficients
    double RYSQ;            // Multiple correlation coefficient
    double SDV;             // Standard deviation of errors
    double FReg;            // Fisher F statistic for regression
    double[] Ycalc;         // Calculated values of Y
    double[] DY;            // Residual values of Y

    public boolean Regress(double[] Y, double[][] X, double[] W) {
        // Y[j]   = j-th observed data point
        // X[i,j] = j-th value of the i-th independent varialble
        // W[j]   = j-th weight value

        int M = Y.length;             // M = Number of data points
        int N = X[0].length;         // N = Number of linear terms
        int NDF = M - N;              // Degrees of freedom
        Ycalc = new double[M];
        DY = new double[M];
        // If not enough data, don't attempt regression
        if (NDF < 1) {
            return false;
        }
        V = new double[N][N];
        C = new double[N];
        SEC = new double[N];
        double[] B = new double[N];   // Vector for LSQ

        // Clear the matrices to start out
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                V[i][j] = 0;
            }
        }

        // Form Least Squares Matrix
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                V[i][j] = 0;
                for (int k = 0; k < M; k++) {
                    V[i][j] = V[i][j] + W[k] * X[k][i] * X[k][j];
                }
            }
            B[i] = 0;
            for (int k = 0; k < M; k++) {
                B[i] = B[i] + W[k] * X[k][i] * Y[k];
            }
        }
        // V now contains the raw least squares matrix
        if (!SymmetricMatrixInvert(V)) {
            return false;
        }
        // V now contains the inverted least square matrix
        // Matrix multpily to get coefficients C = VB
        for (int i = 0; i < N; i++) {
            C[i] = 0;
            for (int j = 0; j < N; j++) {
                C[i] = C[i] + V[i][j] * B[j];
            }
        }

        // Calculate statistics
        double TSS = 0;
        double RSS = 0;
        double YBAR = 0;
        double WSUM = 0;
        for (int k = 0; k < M; k++) {
            YBAR = YBAR + W[k] * Y[k];
            WSUM = WSUM + W[k];
        }
        YBAR = YBAR / WSUM;
        for (int k = 0; k < M; k++) {
            Ycalc[k] = 0;
            for (int i = 0; i < N; i++) {
                Ycalc[k] = Ycalc[k] + C[i] * X[k][i];
            }
            DY[k] = Ycalc[k] - Y[k];
            TSS = TSS + W[k] * (Y[k] - YBAR) * (Y[k] - YBAR);
            RSS = RSS + W[k] * DY[k] * DY[k];
        }
        double SSQ = RSS / NDF;
        RYSQ = 1 - RSS / TSS;
        FReg = 9999999;
        if (RYSQ < 0.9999999) {
            FReg = RYSQ / (1 - RYSQ) * NDF / (N - 1);
        }
        SDV = Math.sqrt(SSQ);

        // Calculate var-covar matrix and std error of coefficients
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                V[i][j] = V[i][j] * SSQ;
            }
            SEC[i] = Math.sqrt(V[i][i]);
        }
        return true;
    }

    public boolean SymmetricMatrixInvert(double[][] V) {
        int N = V.length;
        double[] t = new double[N];
        double[] Q = new double[N];
        double[] R = new double[N];
        double AB;
        int K, L, M;

        // Invert a symetric matrix in V
        for (M = 0; M < N; M++) {
            R[M] = 1;
        }
        K = 0;
        for (M = 0; M < N; M++) {
            double Big = 0;
            for (L = 0; L < N; L++) {
                AB = Math.abs(V[L][L]);
                if ((AB > Big) && (R[L] != 0)) {
                    Big = AB;
                    K = L;
                }
            }
            if (Big == 0) {
                return false;
            }
            R[K] = 0;
            Q[K] = 1 / V[K][K];
            t[K] = 1;
            V[K][K] = 0;
            if (K != 0) {
                for (L = 0; L < K; L++) {
                    t[L] = V[L][K];
                    if (R[L] == 0) {
                        Q[L] = V[L][K] * Q[K];
                    } else {
                        Q[L] = -V[L][K] * Q[K];
                    }
                    V[L][K] = 0;
                }
            }
            if ((K + 1) < N) {
                for (L = K + 1; L < N; L++) {
                    if (R[L] != 0) {
                        t[L] = V[K][L];
                    } else {
                        t[L] = -V[K][L];
                    }
                    Q[L] = -V[K][L] * Q[K];
                    V[K][L] = 0;
                }
            }
            for (L = 0; L < N; L++) {
                for (K = L; K < N; K++) {
                    V[L][K] = V[L][K] + t[L] * Q[K];
                }
            }
        }
        M = N;
        L = N - 1;
        for (K = 1; K < N; K++) {
            M = M - 1;
            L = L - 1;
            for (int J = 0; J <= L; J++) {
                V[M][J] = V[J][M];
            }
        }
        return true;
    }

    /**An example.*/
    public static void main(String[] args) {
        double[][] x = {{1, 38}, {1, 56}, {1, 59}, {1, 64}, {1, 74}};
        double[] y = {41, 63, 70, 72, 84};
        double[] w = {1, 1, 1, 1, 1};
        LinearRegression linReg = new LinearRegression();
        if (linReg.Regress(y, x, w)) {
            System.out.println(linReg.C[0] + " " + linReg.C[1]);
            //  System.out.println(lr.waldTestSlopeP());

        }
        double[] x1 = {38, 56, 59, 64, 74};
        LeastSquares pi = new LeastSquares(x1, y, 1);
        double[] coef = pi.getCoefficients();
        System.out.println(pi.evaluate(56));
        System.out.println(coef[0] + " " + coef[1]);




        try {
            DoubleArrayList xList = new DoubleArrayList();
            DoubleArrayList yList = new DoubleArrayList();
            BufferedReader br = new BufferedReader(new FileReader("lencount.txt"));
            String line = null;

            br.readLine();
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.trim().length() == 0) {
                    continue;
                }
                String[] cells = line.split("\t");
                xList.add(Double.parseDouble(cells[0]));
                yList.add(Double.parseDouble(cells[1]));
            }
            br.close();
            int num = xList.size();
            x1 = new double[num];
            y = new double[num];
            for (int i = 0; i < num; i++) {
                x1[i] = xList.getQuick(i);
                y[i] = yList.getQuick(i) * 10000;
            }
            double rho = 0.6;
            SmoothingCubicSpline fit = new SmoothingCubicSpline(x1, y, rho);
            System.out.println(fit.evaluate(0.6));
            BSpline bfit = new BSpline(x1, y, 100);
            System.out.println(bfit.evaluate(14));
        } catch (Exception nex) {
            nex.printStackTrace();
            //throw new Exception(info);
        }
    }
}
