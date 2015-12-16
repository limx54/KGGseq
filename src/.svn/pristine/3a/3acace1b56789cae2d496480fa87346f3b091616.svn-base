/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.stat;

/**
 *
 * @author mxli
 */
public class LogisticRegressionResultSet {

    double nullm2LnL;
    double convergedm2LnL;
    double overalChisqr;
    double overalP;
    double[] coefficients;
    double[] coefSE;
    double[] coefP;
    //Odds Ratios and 95% Confidence Intervals..
    double[] varOR;
    double[] varLowOR;
    double[] varHighOR;
    double[][] covariances;

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();

        sb.append("Overall Model Fit...");
        sb.append('\n');
        sb.append("-2 Log Likelihood = ");
        sb.append(this.nullm2LnL);
        sb.append(" (Null Model)");
        sb.append('\n');
        sb.append("-2 Log Likelihood = ");
        sb.append(this.convergedm2LnL);
        sb.append(" (Full Model)");
        sb.append('\n');
        sb.append("Chi Square= ");
        sb.append(this.overalChisqr);
        sb.append("\tdf= ");
        sb.append(this.varOR.length);
        sb.append('\n');
        sb.append("p value= ");
        sb.append(this.overalP);
        sb.append('\n');

        sb.append('\n');
        sb.append("Coefficients and Standard Errors...");
        sb.append('\n');
        sb.append("Variable\tCoeff.\tStdErr\tp");
        sb.append('\n');

        for (int i = 0; i < this.coefficients.length; i++) {
            sb.append(i);
            sb.append('\t');
            sb.append(this.coefficients[i]);
            sb.append('\t');
            sb.append(this.coefSE[i]);
            sb.append('\t');
            sb.append(this.coefP[i]);

            sb.append('\n');
        }
        sb.append('\n');
        sb.append("Odds Ratios and 95% Confidence Intervals...");
        sb.append('\n');
        sb.append("Variable\tO.R.\tLow\tHigh");
        sb.append('\n');

        for (int i = 0; i < this.varOR.length; i++) {
            sb.append(i + 1);
            sb.append('\t');
            sb.append(this.varOR[i]);
            sb.append('\t');
            sb.append(this.varLowOR[i]);
            sb.append('\t');
            sb.append(this.varHighOR[i]);
            sb.append('\n');
        }

        /*
        sb.append("Covariances:");
        sb.append('\n');
        for (int i = 0; i < this.covariances.length; i++)
        {
        for (int j = 0; j < this.covariances[i).length; j++)
        sb.append(this.covariances[i)[j)); sb.append('\t');
        sb.append('\n');
        }
         */
        return sb.toString();
    }
}
