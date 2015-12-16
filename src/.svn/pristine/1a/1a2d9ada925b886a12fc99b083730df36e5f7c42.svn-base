/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.stat;

import cern.jet.stat.Probability;
import org.cobi.util.text.Util;

/**
 *
 * @author MX Li
 */
public class SimpleLinearRegression {

    /**An example.*/
    public static void main(String[] args) {
        double[] x = {38, 56, 59, 64, 74};
        double[] y = {41, 63, 70, 72, 84};
        SimpleLinearRegression lr = new SimpleLinearRegression(x, y);
        lr.compute();
        System.out.println(lr.getRoundedModel());
        System.out.println(lr.waldTestSlopeP());
        System.out.println("calculate y given an x of 38 " + lr.calculateY(38));
        System.out.println("calculate x given a y of 41 " + lr.calculateX(41));
    }
    //fields
    private double[] x;
    private double[] y;
    private double slope;
    private double intercept;
    double R2;
    double seSlope;
    double seIntercept;

    public SimpleLinearRegression() {
    }

    public double[] getY() {
        return y;
    }

    public void setY(double[] y) {
        this.y = y;
    }

    public void setX(double[] x) {
        this.x = x;
    }

    //constructor
    public SimpleLinearRegression(double[] x, double[] y) {
        this.x = x;
        this.y = y;
    }

    //methods
    /**Performs linear regression*/
    public void compute() {
        int n = x.length;
        double sumy = 0.0,
                sumx = 0.0,
                sumx2 = 0.0;


        for (int i = 0; i < n; i++) {
            sumx += x[i];
            sumx2 += x[i] * x[i];
            sumy += y[i];
        }
        double xbar = sumx / n;
        double ybar = sumy / n;

        // second pass: compute summary statistics
        double xxbar = 0.0, yybar = 0.0, xybar = 0.0;
        for (int i = 0; i < n; i++) {
            xxbar += (x[i] - xbar) * (x[i] - xbar);
            yybar += (y[i] - ybar) * (y[i] - ybar);
            xybar += (x[i] - xbar) * (y[i] - ybar);
        }
        slope = xybar / xxbar;
        intercept = ybar - slope * xbar;


        // analyze results
        int df = n - 2;
        double rss = 0.0;      // residual sum of squares
        double ssr = 0.0;      // regression sum of squares
        for (int i = 0; i < n; i++) {
            double fit = slope * x[i] + intercept;
            rss += (fit - y[i]) * (fit - y[i]);
            ssr += (fit - ybar) * (fit - ybar);
        }
        R2 = ssr / yybar;
        double svar = rss / df;
        seSlope = svar / xxbar;
        seIntercept = svar / n + xbar * xbar * seSlope;
        seSlope = Math.sqrt(seSlope);
        seIntercept = Math.sqrt(seIntercept);

        /*
        System.out.println("R^2                 = " + R2);
        System.out.println("std error of beta_1 = " + seSlope);
        System.out.println("std error of beta_0 = " + seIntercept);
        seIntercept = svar * sumx2 / (n * xxbar);
        System.out.println("std error of beta_0 = " + Math.sqrt(seIntercept));
        
        System.out.println("SSTO = " + yybar);
        System.out.println("SSE  = " + rss);
        System.out.println("SSR  = " + ssr);
         */
    }

    public double waldTestSlopeP() {
        double s = slope / seSlope;
        return Probability.chiSquareComplemented(1, s * s);
    }

    public double waldTestSlopeZ() {
        double s = slope / seSlope;
        return s;
    }
    //getters

    public double getSlope() {
        return slope;
    }

    public double getIntercept() {
        return intercept;
    }

    public double getRSquared() {

        return R2;
    }

    public double[] getX() {
        return x;
    }

    /**Returns Y=mX+b with full precision, no rounding of numbers.*/
    public String getModel() {
        return "Y= " + slope + "X + " + intercept + " RSqrd=" + getRSquared();
    }

    /**Returns Y=mX+b */
    public String getRoundedModel() {
        return "Y= " + Util.doubleToString(slope, 3) + "X + " + Util.doubleToString(intercept, 3) + " RSqrd=" + Util.doubleToString(getRSquared(), 3);
    }

    /**Calculate Y given X.*/
    public double calculateY(double x) {
        return slope * x + intercept;
    }

    /**Calculate X given Y.*/
    public double calculateX(double y) {
        return (y - intercept) / slope;
    }

    /**Nulls the x and y arrays.  Good to call before saving.*/
    public void nullArrays() {
        x = null;
        y = null;
    }
}
