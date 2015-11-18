/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.stat;

import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.stat.Probability;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.cobi.util.text.Util;

/**
 *
 * @author mxli
 */

public class LogisticRegression {
    
    int nP;
    int nInd;
    double[] coef;
    //RealMatrix covarianceMatrix;
    DoubleMatrix2D covarianceMatrix;
    double[] Y;
    double[][] X;
    boolean toDebug = false;
    BufferedWriter debugOut;
    double priorCorrectionFactor = 1.0;
    
    public LogisticRegression() {
        nP = 0;
        nInd = 0;
        
    }
    
    public void reset() {
        nP = 0;
        nInd = 0;
    }
    
    public double[][] getX() {
        return X;
    }
    
    public void setX(double[][] X) {
        this.X = X;
    }
    
    public double[] getY() {
        return Y;
    }
    
    public void setY(double[] Y) {
        this.Y = Y;
    }
    
    public double getPriorCorrectionFactor() {
        return priorCorrectionFactor;
    }
    
    public void setPriorCorrectionFactor(double priorCorrectionFactor) {
        this.priorCorrectionFactor = priorCorrectionFactor;
    }
    
    public void setVariables(List<double[]> indivList) {
        int n = indivList.size();
        if (n == 0) {
            return;
        }
        DoubleArrayList tmpY = new DoubleArrayList();
        List<double[]> tmpX = new ArrayList<double[]>();
        int varNum = indivList.get(0).length;
        ////  intercept is in the first column
        nP = varNum;
        nInd = indivList.size();
        Y = new double[nInd];
        X = new double[nInd][nP];
        
        for (int i = 0; i < n; i++) {
            double[] cells = indivList.get(i);
            Y[i] = cells[varNum - 1];
            X[i][0] = 1;
            System.arraycopy(cells, 0, X[i], 1, varNum - 1);
            //  double[] tmpXRow = new double[nP];
            // tmpXRow[0] = 1.0;
            //  tmpX.add(tmpXRow);
        }
        
    }
    
    void standardiseIndependent() {
        double[] xMean = new double[nP];           // Attribute means
        double[] xSD = new double[nP];           // Attribute stddev's
        Arrays.fill(xMean, 0.0);
        Arrays.fill(xSD, 0.0);
        for (int i = 0; i < nInd; i++) {
            for (int j = 1; j < nP; j++) {
                xMean[j] += X[i][j];
                xSD[j] += X[i][j] * X[i][j];
            }
        }
        
        xMean[0] = 0;
        xSD[0] = 1;
        for (int j = 1; j < nP; j++) {
            xMean[j] = xMean[j] / nInd;
            if (nInd > 1) {
                xSD[j] = Math.sqrt(Math.abs(xSD[j] - nInd * xMean[j] * xMean[j]) / (nInd - 1));
            } else {
                xSD[j] = 0;
            }
        }

        // Normalise input data
        for (int i = 0; i < nInd; i++) {
            for (int j = 1; j < nP; j++) {
                if (xSD[j] != 0) {
                    X[i][j] = (X[i][j] - xMean[j]) / xSD[j];
                }
            }
        }
    }

    /*
    test results by R
    path<-"D:/home/mxli/MyJava/GenetSimulator/debug.txt";
    path<-"D:/home/mxli/MyJava/GenetSimulator/debug.txt";
    dat<-read.table(path,sep = "\t");
    result1 <- glm(dat$V1 ~ dat$V2 + dat$V3+ dat$V4+ dat$V5 + dat$V6+dat$V7, family=binomial(logit),trace =TRUE)
    summary(result1)
    1-pchisq(result1$null.deviance-result1$deviance,result1$df.null-result1$df.residual)
     */
    public boolean fitLM() throws Exception {
        if (Y == null || Y.length == 0) {
            System.err.println("No dependent variables!!!");
            return false;
        }
        nInd = Y.length;
        
        if (X == null || X.length == 0 || X[0].length == 0) {
            System.err.println("No independent variables!!!");
            return false;
        }
        
        nP = X[0].length;
        coef = new double[nP];
        Arrays.fill(coef, 0.0);
        //toDebug = true;
        if (toDebug) {
            debugOut = new BufferedWriter(new FileWriter("debug.txt"));
            for (int i = 0; i < nInd; i++) {
                debugOut.write(Y[i] + "\t");
                for (int j = 1; j < nP; j++) {
                    debugOut.write(X[i][j] + "\t");
                }
                debugOut.write("\n");
            }
            debugOut.close();
        }


        ///////////////////////////////////////
        // Newton-Raphson to fit logistic model
        boolean converge = false;
        int it = 0;
        double[] p = new double[nInd];
        double[] V = new double[nInd];

        //covarianceMatrix = new Array2DRowRealMatrix(nP, nP);
        covarianceMatrix = new DenseDoubleMatrix2D(nP, nP);
        double[][] T2 = new double[nP][nInd];
        double[] t3 = new double[nInd];
        double[] ncoef = new double[nP];
        //LUDecomposition solver = null;
        Algebra algebra = new Algebra();
        double delta = 0;
        double sum;
        double t;
        int maxIterTime = 20;
        double tolerateDiff = 1e-6;
        
        
        while (!converge) {
            // Determine p and V
            for (int i = 0; i < nInd; i++) {
                t = 0;
                for (int j = 0; j < nP; j++) {
                    t += coef[j] * X[i][j];
                }
                p[i] = 1 / (1 + Math.exp(-t));
                V[i] = p[i] * (1 - p[i]);
            }

            // Update coefficients
            // b <- b +  solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p )
            for (int j = 0; j < nP; j++) {
                for (int k = j; k < nP; k++) {
                    sum = 0;
                    for (int i = 0; i < nInd; i++) {
                        sum += X[i][j] * V[i] * X[i][k];
                    }
                    covarianceMatrix.setQuick(j, k, sum);
                    covarianceMatrix.setQuick(k, j, sum);
                    //covarianceMatrix.setEntry(j, k, sum);
                    //covarianceMatrix.setEntry(k, j, sum);
                }
            }
            // System.out.println(covarianceMatrix.toString());
            //solver = new LUDecompositionImpl(covarianceMatrix);
            //covarianceMatrix = solver.getSolver().getInverse();
            covarianceMatrix = algebra.inverse(covarianceMatrix);
            // Resize and set elements to 0
            for (int i = 0; i < nP; i++) {
                Arrays.fill(T2[i], 0.0);
            }
            // note implicit transpose of X
            for (int i = 0; i < nP; i++) {
                for (int j = 0; j < nInd; j++) {
                    for (int k = 0; k < nP; k++) {
                        //T2[i][j] += (covarianceMatrix.getEntry(i, k) * X[j][k]);
                        T2[i][j] += (covarianceMatrix.getQuick(i, k) * X[j][k]);
                    }
                }
            }
            
            for (int i = 0; i < nInd; i++) {
                t3[i] = Y[i] - p[i];
            }
            
            Arrays.fill(ncoef, 0.0);
            for (int j = 0; j < nP; j++) {
                for (int i = 0; i < nInd; i++) {
                    ncoef[j] += T2[j][i] * t3[i];
                }
            }

            // Update coefficients, and check for
            // convergence
            delta = 0;
            for (int j = 0; j < nP; j++) {
                delta += Math.abs(ncoef[j]);
                coef[j] += ncoef[j];
            }
            if (delta < tolerateDiff) {
                converge = true;
            }
            // Next iteration
            it++;
            if (it > maxIterTime && !converge) {
                System.err.println("Iteration over " + maxIterTime + " times, aborted!");
                return false;
            }
            //System.out.println(coef[0] + "\t" + coef[1] + "\t" + coef[2]);
        }

        /*
        /////////////////////////////////////////
        // Obtain covariance matrix of estimates
        // S <- solve( t(X) %*% V %*% X )
        // Transpose X and multiple by diagonal V
        double[][] Xt = new double[nP][nInd];
        for (int i = 0; i < nInd; i++) {
        for (int j = 0; j < nP; j++) {
        Xt[j][i] = X[i][j] * V[i];
        }
        }
        
        S = multMatrix(Xt, X);
        solver = new LUDecompositionImpl(S);
        S = solver.getSolver().getInverse();
        System.out.print("Std. Error\n");
        for (int j = 0; j < endIndex; j++) {
        System.out.print(Math.sqrt(S.getEntry(j, j)));
        System.out.print("\t");
        }
        System.out.println(Math.sqrt(S.getEntry(endIndex, endIndex)));
         */
        // if ( cluster )      HuberWhite();
        if (toDebug) {
            int endIndex = nP - 1;
            System.out.print("Beta\n");
            for (int j = 0; j < endIndex; j++) {
                System.out.print(coef[j]);
                System.out.print("\t");
            }
            System.out.println(coef[endIndex]);
            
            
            System.out.print("Std. Error\n");
            for (int j = 0; j < endIndex; j++) {
                //System.out.print(Math.sqrt(covarianceMatrix.getEntry(j, j)));
                System.out.print(Math.sqrt(covarianceMatrix.getQuick(j, j)));
                System.out.print("\t");
            }
            //System.out.println(Math.sqrt(covarianceMatrix.getEntry(endIndex, endIndex)));
            System.out.println(Math.sqrt(covarianceMatrix.getQuick(endIndex, endIndex)));
            /*
            System.out.print("Sigma\n");
            for (int i = 0; i < nP; i++) {
            for (int j = 0; j < endIndex; j++) {
            System.out.print(S.getEntry(i, j));
            System.out.print("\t");
            }
            System.out.println(S.getEntry(i, endIndex));
            }
             */
            
            System.out.print(getModelPValue() + "\n");
        }
        return true;
    }
    
    public double getLnLk() {
        // Return  sample log-likelihood
        // We assume the model is fit, and all Y's are either 0 or 1
        double lnlk = 0;
        for (int i = 0; i < nInd; i++) {
            double t = 0;
            for (int j = 0; j < nP; j++) {
                t += coef[j] * X[i][j];
            }
            lnlk += ((Y[i] == 1) ? Math.log(1 / (1 + Math.exp(-t))) : Math.log(1 - (1 / (1 + Math.exp(-t)))));
        }
        return lnlk;
    }
    
    public double getModelPValue() throws Exception {
        // Return p value of the overall model
        // We assume the model is fit, and all Y's are either 0 or 1
        double pi0 = 0.0;   //frequency of classes
        for (int i = 0; i < nInd; i++) {
            pi0 = pi0 + Y[i];
        }
        pi0 = pi0 / nInd;

        //Log-likelihood with constant in the model; Otherwise it is null_LL =totWeights*ln(0.5)=-0.6931472totWeights
        double null_LL = nInd * (pi0 * Math.log(pi0 / (1 - pi0)) + Math.log((1 - pi0)));
        
        
        double chiSqure = 2 * (getLnLk() - null_LL);
        double pValue = 1;
        pValue = Probability.chiSquareComplemented(nP - 1, chiSqure);

        /*
        System.out.println("null -2LL "+(-2*null_LL));
        System.out.println("alter -2LL "+(-2*getLnLk()));
        System.out.println("df "+(nP - 1));
        System.out.println("chiSqure "+(chiSqure));
         * 
         */
        return pValue;
    }
    
    public double getCoefPValue(int testParameter) throws Exception {
        double[] var = getVar();
        double se = Math.sqrt(var[testParameter]);
        double Z = coef[testParameter] / se;        
        return Probability.chiSquareComplemented(1, Z * Z);        
    }
    
    public double getCoef(int testParameter) throws Exception {
        return covarianceMatrix.getQuick(testParameter, testParameter);
    }
    
    public double getNagelkerkeRSquare() {
        double pi0 = 0.0;   //frequency of classes
        for (int i = 0; i < nInd; i++) {
            pi0 = pi0 + Y[i];
        }
        pi0 = pi0 / nInd;

        //Log-likelihood with constant in the model; Otherwise it is null_LL =totWeights*ln(0.5)=-0.6931472totWeights
        double null_LL = nInd * (pi0 * Math.log(pi0 / (1 - pi0)) + Math.log((1 - pi0)));
        //Nagelkerke R squre
        double r = 1 - Math.pow(Math.exp(null_LL - getLnLk()), 2.0 / nInd);
        r = r / (1 - Math.pow(Math.exp(null_LL), 2.0 / nInd));
        return r;
    }
    
    public double[] getVar() {
        double[] var = new double[nP];
        for (int i = 0; i < nP; i++) {
            //var[i] = covarianceMatrix.getEntry(i, i);
            var[i] = covarianceMatrix.getQuick(i, i);
        }
        return var;
    }
    
    public double getGivenProbability(double[] scores) {
        double p = coef[0];
        for (int i = 1; i < nP; i++) {
            p += (coef[i] * scores[i - 1]);
        }
        p = 1 + priorCorrectionFactor * Math.exp(-p);
        
        return 1.0 / p;
    }
    
    public double[] getSE() {
        double[] var = new double[nP];
        for (int i = 0; i < nP; i++) {
            //var[i] = Math.sqrt(covarianceMatrix.getEntry(i, i));
            var[i] = Math.sqrt(covarianceMatrix.getQuick(i, i));
        }
        return var;
    }
    
    public double[] getCoefs() {
        return coef;
    }

    /**
     * Gets a string describing the classifier.
     *
     * @return a string describing the classifer built.
     */
    @Override
    public String toString() {
        double[] coef_SE = getSE();
        int df = nP;
        StringBuffer result = new StringBuffer();
        if (coef == null) {
            return result.append(": No model built yet.").toString();
        }
        
        result.append("\nCoefficients...\n" + "Variable      Coeff.      td. Error       z value      Pr(>|z|)\n");
        
        try {
            for (int j = 1; j < nP; j++) {
                result.append(Util.doubleToString(j, 8, 0));
                result.append(" " + Util.doubleToString(coef[j], 12, 4));
                result.append(" " + Util.doubleToString(coef_SE[j], 12, 4));
                result.append(" " + Util.doubleToString(coef[j] / coef_SE[j], 12, 4));
                result.append(" " + Util.doubleToString(Probability.chiSquareComplemented(1, Math.pow(coef[0] / coef_SE[0], 2)), 12, 4));
                result.append("\n");
            }
            
            result.append("Intercept ");
            
            result.append(" " + Util.doubleToString(coef[0], 10, 4));
            result.append(" " + Util.doubleToString(coef_SE[0], 10, 4));
            result.append(" " + Util.doubleToString((coef[0] / coef_SE[0]), 10, 4));
            result.append(" " + Util.doubleToString(Probability.chiSquareComplemented(1, Math.pow(coef[0] / coef_SE[0], 2)), 10, 4));
            
            result.append("\n");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        result.append("\nOdds Ratios...\n" + "Variable         O.R.\n");
        for (int j = 1; j < nP; j++) {
            result.append(Util.doubleToString(j, 8, 0));
            
            double ORc = Math.exp(coef[j]);
            result.append(" " + ((ORc > 1e10) ? "" + ORc : Util.doubleToString(ORc, 12, 4)));
            
            result.append("\n");
        }
        result.append("Significance of the model with these variable(s) ");
        
        result.append(" p-value ");
        try {
            result.append(getModelPValue());
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        result.append("\n");
        return result.toString();
    }
    /*
    private RealMatrix multMatrix(double[][] a, double[][] b) {
    
    int ar = a.length;
    int br = b.length;
    if (ar == 0 || br == 0) {
    System.err.println("Internal error: multiplying 0-sized matrices");
    }
    
    int ac = a[0].length;
    int bc = b[0].length;
    if (ac != br) {
    System.err.println("Internal error: non-conformable matrices in multMatrix()");
    }
    
    int cr = ar;
    int cc = bc;
    
    
    RealMatrix c = new Array2DRowRealMatrix(cr, cc);
    for (int i = 0; i < ar; i++) {
    for (int j = 0; j < bc; j++) {
    c.setEntry(i, j, 0);
    }
    }
    
    for (int i = 0; i < ar; i++) {
    for (int j = 0; j < bc; j++) {
    for (int k = 0; k < ac; k++) {
    c.setEntry(i, j, c.getEntry(i, j) + a[i][k] * b[k][j]);
    }
    }
    }
    return c;
    }
     */
}
