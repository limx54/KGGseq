// (c) 2009-2011 Miaoxin Li
// This file is distributed as part of the KGG source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
// Permission is granted for you to use this file to compile IGG.
// All computer programs have bugs. Use this file at your own risk.
// Tuesday, March 01, 2011
package org.cobi.util.stat;

import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import cern.jet.stat.Probability;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import umontreal.iro.lecuyer.probdist.BetaDist;
import umontreal.iro.lecuyer.probdist.DiscreteDistributionInt;
import umontreal.iro.lecuyer.probdist.HypergeometricDist;

/**
 *
 * @author mxli
 */
public class MultipleTestingMethod {

    //http://courses.ttu.edu/isqs6348-westfall/images/6348/BonHolmBenHoch.htm
    //Problem:  The ordinary Bonferroni method is conservative:  That is, the true FWE is less than .05 because of the Bonferroni inequality.  This means that the method does not reject hypotheses as often as it should and therefore lacks power.
    //The Bonferroni-Holm method allows more rejections, and is therefore less conservative and more powerful than the Bonferroni method.
    //, which is called a “step-down method,”
    public static double bonferroniHolmFWE(double fdrThreshold, DoubleArrayList pValues) throws Exception {
        int snpSize = pValues.size();
        if (snpSize == 0) {
            return 1;
        }
        pValues.quickSort();
        int i = 0;

        while (i < snpSize) {
            if (pValues.getQuick(i) > fdrThreshold / (snpSize - i)) {
                return fdrThreshold / (snpSize - i);
            }
            i++;
        }
        //it must be less than or equal to this value
        return fdrThreshold;

    }
//The Benjamini-Hochberg method aims to control the FDR, or “False Discovery Rate,” rather than the FWE or “FamilyWise Error rate.”
//FDR is defined as
//FDR = E{ (false rejections)/(total rejections) }.   
// Benjamini-Hochberg method is a “step-up” metho
    //make sure the pValues have been sorted by ascending

    /*
     public static double benjaminiHochbergFDR(double fdrThreshold, DoubleArrayList pValues) throws Exception {
     int i;
     int snpSize = pValues.size();
     if (snpSize == 0) {
     return 1;
     }
     pValues.quickSort();

     //fdrThreshold = fdrThreshold / (snpSize*(Math.log(snpSize)+0.577721));
     fdrThreshold = fdrThreshold / (snpSize);
     i = snpSize - 1;
     while (i >= 0) {
     if (pValues.getQuick(i) <= (i + 1) * fdrThreshold) {
     return fdrThreshold * (i + 1);
     }
     i--;
     }

     //it must be less than or equal to this value
     return fdrThreshold;

     }
     */
    //This is derived from PLINK
    //Benjamini & Hochberg (1995)
    public static double benjaminiHochbergFDR(double fdrThreshold, DoubleArrayList sp) {

        int ti = sp.size();
        if (ti == 0) {
            return fdrThreshold;
        }
        // BH 

        double[] pv_BH = new double[ti];
        double t = (double) ti;

        pv_BH[ti - 1] = sp.getQuick(ti - 1);
        double x = 0;
        for (int i = ti - 2; i >= 0; i--) {
            x = (t / (double) (i + 1)) * sp.getQuick(i) < 1 ? (t / (double) (i + 1)) * sp.getQuick(i) : 1;
            pv_BH[i] = pv_BH[i + 1] < x ? pv_BH[i + 1] : x;
        }
        if (pv_BH[0] <= fdrThreshold) {
            for (int i = 1; i < ti; i++) {
                if (pv_BH[i] >= fdrThreshold) {
                    return sp.getQuick(i - 1);
                }
            }
        }
        return fdrThreshold / ti;
    }

    //This is derived from PLINK
    //Benjamini & Yekutieli (2001) ("BY")
    public static double BenjaminiYekutieliFDR(double fdrThreshold, DoubleArrayList sp) {
        int ti = sp.size();
        if (ti == 0) {

            return fdrThreshold;
        }
        // BY 

        double[] pv_BY = new double[ti];
        double a = 0;
        double t = (double) ti;
        for (double i = 1; i <= t; i++) {
            a += 1 / i;
        }
        pv_BY[ti - 1] = a * sp.getQuick(ti - 1) < 1 ? a * sp.getQuick(ti - 1) : 1;

        for (int i = ti - 2; i >= 0; i--) {
            double x = ((t * a) / (double) (i + 1)) * sp.getQuick(i) < 1 ? ((t * a) / (double) (i + 1)) * sp.getQuick(i) : 1;
            pv_BY[i] = pv_BY[i + 1] < x ? pv_BY[i + 1] : x;
        }
        if (pv_BY[0] <= fdrThreshold) {
            for (int i = 1; i < ti; i++) {
                if (pv_BY[i] >= fdrThreshold) {
                    return sp.getQuick(i - 1);
                }
            }
        }
        return fdrThreshold / ti;
    }

    //make sure the pValues have been sorted by ascending

    public static double adaptiveBenjaminiFDR(double fdrThreshold, DoubleArrayList pValues) {
        int pvSize = pValues.size();
        int i;
        pValues.quickSort();
        // stage I
        fdrThreshold = fdrThreshold / (1 + fdrThreshold);
        fdrThreshold = fdrThreshold / pvSize;

        i = pvSize - 1;
        while ((i >= 0) && (pValues.getQuick(i) > (i + 1) * fdrThreshold)) {
            i--;
        }

        if (i <= 0) {
            return fdrThreshold;
        } else if (i >= pvSize) {
            return fdrThreshold * (pvSize);
        }

        // stage II
        int m0 = pvSize - i;
        fdrThreshold = fdrThreshold * pvSize / m0;

        i = pvSize - 1;
        while ((i >= 0) && (pValues.getQuick(i) > (i + 1) * fdrThreshold)) {
            i--;
        }
        //it must be less than or equal to this value
        return fdrThreshold * (i + 1);
    }

    //make sure the pValues have been sorted by ascending
    public static double storeyFDRTest(double alpha, DoubleArrayList pValues) {
        int snpSize = pValues.size();

        double a, tmpDouble;
        int i;

        tmpDouble = 0;
        a = 0;
        i = 0;

        while (pValues.getQuick(i) <= 0.5) {
            i++;
            if (i == snpSize) {
                break;
            }
        }
        a = i;
        a /= snpSize;
        a = (a - 0.5) / (1 - 0.5);
        //adjustedAlpha = alpha / (1 - a);
        tmpDouble = alpha / ((1 - a) * snpSize);

        i = 0;
        while ((i < snpSize) && (pValues.getQuick(i) <= (i + 1) * tmpDouble)) {
            i++;
        }
        //it must be less than this value
        return tmpDouble * (i + 1);
    }

    public static void zScores(final DoubleArrayList pValues) {
        double q = 1;

        int size = pValues.size();

        // assume they are two-tailed I2-values
        for (int i = 0; i < size; i++) {
            if (Double.isNaN(pValues.getQuick(i))) {
                continue;
            }
            q = 1 - pValues.getQuick(i);
            //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
            if (q > 0.5) {
                q = 1 - q;
                if (q < 1E-323) {
                    q = 1E-323;
                }
                q = Probability.normalInverse(q);
                pValues.setQuick(i, -q);

            } else {
                if (q < 1E-323) {
                    q = 1E-323;
                }
                q = Probability.normalInverse(q);
                pValues.setQuick(i, q);
            }
        }
    }

    public static double zScore(final double pValue) {
        double q = 0;
        // assume they are two-tailed I2-values 
        if (Double.isNaN(pValue)) {
            return q;
        }
        q = 1 - pValue;
        //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
        if (q > 0.5) {
            q = 1 - q;
            if (q < 1E-323) {
                q = 1E-323;
            }
            q = Probability.normalInverse(q);
            return -q;

        } else {
            if (q < 1E-323) {
                q = 1E-323;
            }
            q = Probability.normalInverse(q);
            return q;
        }

    }

    //http://www.unm.edu/~marcusj/WilcoxonSR.pdf
    //this is not accurate
    public static double wilcoxonSignedRankTest(double median, DoubleArrayList pValues) {
        int snpSize = pValues.size();
        class PairData {

            double value;
            int sign;

            public PairData(double value, int sign) {
                this.value = value;
                this.sign = sign;
            }
        }
        class PairComparator implements Comparator<PairData> {

            @Override
            public int compare(PairData arg0, PairData arg1) {
                return Double.compare(arg0.value, arg1.value);
            }
        }
        //after the transversion to z scores 
        zScores(pValues);
        System.out.println(pValues.toString());
        List<PairData> diff = new ArrayList<PairData>();

        for (int i = 0; i < snpSize; i++) {
            diff.add(new PairData(Math.abs(pValues.getQuick(i) - median), (pValues.getQuick(i) - median) > 0 ? 1 : -1));
        }
        Collections.sort(diff, new PairComparator());

        double ranks = 0;
        for (int i = 0; i < snpSize; i++) {
            if (diff.get(i).sign > 0) {
                ranks += (i + 1);
            }
        }
        if (snpSize > 9) {
            double expectRanks = 0.25 * snpSize * (snpSize + 1);
            // double se = Math.sqrt(snpSize * (snpSize + 1) * (2 * snpSize + 1) / 6.0);
            // double sdiff = (ranks - 0.5) / se;

            double se = Math.sqrt(snpSize * (snpSize + 1) * (2 * snpSize + 1) / 24);
            double sdiff = (ranks - expectRanks) / se;
//two side test
            double prob = Probability.normal(sdiff);
            if (prob > 0.5) {
                return 2 - 2 * prob;
            } else {
                return 2 * prob;
            }

        } else {
            //Critical Values of ±W for Small Samples
            //http://vassarstats.net/textbook/ch12a.html
            /*
             * Non-Directional Test
             N	--	.05	.02	.01
             5	  15  	  --  	  --  	  --  
             6	  17  	  21  	  --  	  --  
             7	  22  	  24  	  28  	  --  
             8	  26  	  30  	  34  	  36  
             9	  29  	  35  	  39  	  43  
             */
            if (snpSize < 4) {
                return 1;
            } else if (snpSize == 5) {
                if (ranks < 15) {
                    return 1;
                } else {
                    return 0.05;
                }
            } else if (snpSize == 6) {
                if (ranks < 21) {
                    return 1;
                } else {
                    return 0.05;
                }
            } else if (snpSize == 7) {
                if (ranks < 24) {
                    return 1;
                } else if (ranks < 28) {
                    return 0.05;
                } else {
                    return 0.02;
                }
            } else if (snpSize == 8) {
                if (ranks < 30) {
                    return 1;
                } else if (ranks < 34) {
                    return 0.05;
                } else if (ranks < 36) {
                    return 0.02;
                } else {
                    return 0.01;
                }
            } else if (snpSize == 9) {
                if (ranks < 35) {
                    return 1;
                } else if (ranks < 39) {
                    return 0.05;
                } else if (ranks < 43) {
                    return 0.02;
                } else {
                    return 0.01;
                }
            }
            return 1;
        }
    }

    private static double[][] getPoints() {
// The density of the standard normal probability distribution
// points contains one array for X values and one array for Y values
        double[][] points = new double[2][25];
        final double CPI = Math.sqrt(2 * Math.PI);
        for (int i = 0; i < points[0].length; ++i) {
            double x = -3.5 + i * 7.0 / (points[0].length - 1);
            points[0][i] = x;
            points[1][i] = Math.exp(-x * x / 2.0) / CPI;
        }
        return points;
    }

    public static void main(String[] args) {
        try {
            BetaDist betaDist = new BetaDist(300, 383700);
            //Beta betaDist = new Beta(3000, 383700);
            double x = betaDist.cdf(0.331);
            System.out.println(x);

            //System.out.println(Probability.chiSquareComplemented(1, 15.62));
            //MultipleTestingMethod.simulationSimesTest();
            //MultipleTestingMethod.simulationModifiedSidakTest();
            //MultipleTestingMethod.simulationTestEstimateNullHypothesisProportion();
            System.out.println(MultipleTestingMethod.hypergeometricEnrichmentTest(19043, 733, 40, 4));
            System.out.println(MultipleTestingMethod.hypergeometricEnrichmentTest(19043, 1744, 40, 8));

            System.out.println(MultipleTestingMethod.hypergeometricEnrichmentTest(18988, 40, 1000, 8));
            System.out.println(MultipleTestingMethod.hypergeometricEnrichmentTest(18988, 1000, 40, 8));
            //double[] pValues = new double[]{17, 50, 45, 59.8, 21.74, 16, 9, 15.43, 5.12, 40, 35, 13.35, 13.4};
            // double [] pValues=new double[]{8.97E-02, 1.17E-14, 2.84E-06, 8.52E-02, 1.77E-01, 1.00E+00, 2.03E-03, 1.17E-01};
            double[] pValues = new double[]{0.08970206984892933, 1.1721659248909853E-14, 2.8415038031821715E-6, 0.08522038516202218, 0.17709282017413752, 1.0, 0.002032425164707917, 0.11738624009898738};
            DoubleArrayList dList = new DoubleArrayList(pValues);
            System.out.println(MultipleTestingMethod.wilcoxonSignedRankTest(7.38, dList));
            double[] ps = new double[]{5.89E-08, 0.00000336, 0.00000412, 0.0000106, 0.000011, 0.0000139, 0.0000198, 0.0000204, 0.0000219, 0.0000242, 0.0000492, 0.0000495, 0.0000659, 0.000071, 0.0000732, 0.0000892, 0.0000974, 0.000106, 0.000115, 0.000116, 0.000143, 0.000224, 0.00024, 0.000298, 0.000332, 0.000344, 0.000353, 0.000367, 0.000377, 0.000422, 0.000467, 0.000474, 0.000493, 0.000496, 0.000518, 0.000543, 0.00055, 0.000574, 0.000578, 0.00058, 0.000613, 0.361, 0.361, 0.361, 0.361, 0.361, 0.361, 0.361, 0.361, 0.361, 0.361, 0.361, 0.361, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.362, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.363, 0.364, 0.364, 0.364, 0.364, 0.364, 0.364};
            //double[] ps = new double[]{0.1, 0.3, 0.5, 0.8};

            DoubleArrayList plist = new DoubleArrayList(ps);

            double s = MultipleTestingMethod.benjaminiHochbergFDR(0.05, plist);
            System.out.println(s);

            /*
             WilcoxonTest wt = new WilcoxonTest(pValues, 0.5);
             System.out.println(wt.getSP());
             System.out.println(wt.approxSP() + " " + wt.exactSP() + " " + wt.getSP() + " " + wt.getZ());
             // WilcoxonSignedRankTest wst = new WilcoxonSignedRankTest();

             //wst.test(pValues,  alternative = "two.sided", mu = 0, paired = TRUE, exact = TRUE, correct = FALSE);
             // DecimalFormat df = new DecimalFormat("0.00E0", new DecimalFormatSymbols(Locale.US));
             // MultipleTestingMethod.simulationTestWeightedCombinationProbabilities();

             /*
             double[][] points = getPoints();
             XYLineChart chart = new XYLineChart(null, "X", null, points);
             chart.setAutoRange00(true, true); // Axes pass through (0,0)
             chart.toLatexFile("NormalChart.tex", 12, 8);
             * 
             */
            /*
             if (true) {
             MultipleTestingMethod.simulationFDR();
             return;
             }
            
             ChiSquaredDistribution chiDis = new ChiSquaredDistributionImpl(1);
            
            
             double I2 = 1E-311;
             I2 = Probability.normalInverse(0.5 * I2);// two tails to be one tail
             I2 = I2 * I2;
             System.out.println(I2);
             I2 = Probability.chiSquareComplemented(1, I2);
            
             String formattedNumber = df.format(I2, new StringBuffer(), new FieldPosition(NumberFormat.INTEGER_FIELD)).toString();
             System.out.println(formattedNumber);
             double Q = -2;
             I2 = 1 - Probability.normal(Q);
             formattedNumber = df.format(I2, new StringBuffer(), new FieldPosition(NumberFormat.INTEGER_FIELD)).toString();
             System.out.println(formattedNumber);
            
             I2 = Probability.chiSquareComplemented(1, Q * Q) / 2;
             if (Q < 0) {
             I2 = 1 - I2;
             }
            
             double[] weights = {0.9, 0.1};
             double[] pValues = new double[2];
             pValues[0] = 1E-30;
             pValues[1] = 1E-3;
             double x = -2E-12;
             double I2 = (weightedCombinationProbabilities(weights, pValues));
            
            
             int n = 2;
             double p1 = 0;
            
             double s = 1;
             for (int k = n; k >= 1; k--) {
             s *= ((n - k + 1) * 1.0 / k);
             }
            
             for (int i = n; i >= 1; i--) {
             p1 += (s * Math.pow(x, i));
             String formattedNumber = df.format(p1, new StringBuffer(), new FieldPosition(NumberFormat.INTEGER_FIELD)).toString();
             System.out.println(formattedNumber);
             s *= (i / ((n - i + 1) * 1.0));
             }
             String formattedNumber = df.format(-p1, new StringBuffer(), new FieldPosition(NumberFormat.INTEGER_FIELD)).toString();
             System.out.println(formattedNumber);
            
             I2 = 1 - Math.pow(1 + x, n);
             formattedNumber = df.format(I2, new StringBuffer(), new FieldPosition(NumberFormat.INTEGER_FIELD)).toString();
             System.out.println(formattedNumber);
             *
            
             List<double[]> pVs = new ArrayList<double[]>();
             pVs.add(new double[]{2.46E-06, 6.58E-05});
             pVs.add(new double[]{7.97E-06, 1.38E-04});
             pVs.add(new double[]{2.11E-06, 1.21E-03});
             pVs.add(new double[]{3.37E-06, 1});
             pVs.add(new double[]{1.85E-05, 1.35E-02});
             pVs.add(new double[]{1.12E-05, 1});
             pVs.add(new double[]{1.80E-04, 1});
             pVs.add(new double[]{3.96E-06, 4.30E-04});
             pVs.add(new double[]{2.60E-05, 9.15E-08});
             pVs.add(new double[]{4.27E-04, 1});
             pVs.add(new double[]{3.23E-04, 1});
             pVs.add(new double[]{1.00E-04, 1});
             double[] weights = {1, 1};
             for (int i = 0; i < pVs.size(); i++) {
             double p = (weightedCombinationProbabilities(weights, pVs.get(i)));
             String formattedNumber = df.format(p, new StringBuffer(), new FieldPosition(NumberFormat.INTEGER_FIELD)).toString();
             System.out.println(formattedNumber);
             }
             * 
             */
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    //Algorithm proposed by  Statistics & Probability Letters 73 (2005) 179??87; A simple approximation for the distribution of the weighted combination ofnon-independent or independent probabilities
    // and Liptak-Stouffer method  http://en.wikipedia.org/wiki/Fisher%27s_method
    public static double weightedCombinationProbabilities(double[] weights, double[] pValues) {
        double p = 1;
        double w2 = 0;
        double w = 0;
        int size = weights.length;
        double X = 0;
        double Z = 0;
        double q = 0;

        // weights[0] = 1;
        //  weights[1] = 1;
        for (int i = 0; i < size; i++) {
            if (pValues[i] > 0.5) {
                if (pValues[i] >= 1) {
                    //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
                    q = Probability.normalInverse(1E-323);
                } else {
                    q = Probability.normalInverse(1 - pValues[i]);
                }
            } else {
                if (pValues[i] < 1E-320) {
                    pValues[i] = 1E-320;
                }
                q = -Probability.normalInverse(pValues[i]);
            }

            // q =  Probability.normalInverse(1 - pValues[i]);
            w2 += weights[i] * weights[i];
            Z += weights[i] * q;
            //w += weights[i];
            //X += weights[i] * Math.log(pValues[i]);
        }
        Z = Z / Math.sqrt(w2);

        /*
         //Statistics & Probability Letters 73 (2005) 179??87
         double c = w2 / w;
         double f = 2 * w * w / w2;
         I2 = Probability.chiSquareComplemented(f, -2 * X / c);
         */
        //Liptak-Stouffer method  http://en.wikipedia.org/wiki/Fisher%27s_method
        //I2 = 1 - Probability.normal(Q);
        p = Probability.chiSquareComplemented(1, Z * Z) / 2;
        if (Z < 0) {
            p = 1 - p;
        }
        return p;
    }

    public static double combinePValuebyFisherCombinationTest(DoubleArrayList pValueArray) throws Exception {
        int snpSize = pValueArray.size();
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray.getQuick(0);
        }

        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            double p = pValueArray.getQuick(i);
            Y += (-2 * Math.log(p));
        }
        if (Y < 0) {
            Y = 0;
        }
        double p = Probability.chiSquareComplemented(snpSize * 2, Y);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;

    }

    public static double combinationHeterogeneityCochranQTest(final double[] pValues) {
        double p = 1;
        double Q = 0;
        double q = 1;
        double meanQ = 0;
        int size = pValues.length;
        double[] qValues = new double[size];
        // assume they are two-tailed I2-values
        for (int i = 0; i < size; i++) {
            //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
            q = Probability.normalInverse(pValues[i] / 2);
            qValues[i] = -q;
            //Cochran's Q statistic  http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0000841
            // meanQ += qValues[i];
        }

        Q = qValues[0] - qValues[1];
        Q = Q * Q / 2;
        p = Probability.chiSquareComplemented(1, Q);
        //Quantifying heterogeneity in a meta-analysis Julian P. T. Higgins?? ??and Simon G. Thompson  I 2 = 100%?(Q ??df)/Q http://onlinelibrary.wiley.com/doi/10.1002/sim.1186/pdf       
        if (p < 0) {
            p = 0;
        }
        return p;
    }

    public static double combinationHeterogeneityI2(final double[] pValues) {
        double I2 = 1;
        double Q = 0;
        double q = 1;
        double meanQ = 0;
        int size = pValues.length;
        double[] qValues = new double[size];
        // assume they are two-tailed I2-values
        for (int i = 0; i < size; i++) {
            //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
            q = Probability.normalInverse(pValues[i] / 2);
            qValues[i] = -q;
            //Cochran's Q statistic  http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0000841
            // meanQ += qValues[i];
        }

        Q = qValues[0] - qValues[1];
        Q = Q * Q / 2;

        //Quantifying heterogeneity in a meta-analysis Julian P. T. Higgins?? ??and Simon G. Thompson  I 2 = 100%?(Q ??df)/Q http://onlinelibrary.wiley.com/doi/10.1002/sim.1186/pdf
        I2 = 1 - 1 / Q;
        if (I2 < 0) {
            I2 = 0;
        }
        return I2;
    }

    /*
     BufferedWriter bw = new BufferedWriter(new FileWriter("test.txt"));
     for (int i = 0; i < pValues.size(); i++) {
     bw.write(String.valueOf(pValues.getQuick(i)));
     bw.write("\n");
     }
     bw.close();
    
     pp <- scan("D:/home/mxli/MyJava/KGG2Eclipse/test.txt", quiet= TRUE);
     #hist(pp, breaks=1000);
    
     observed <- sort(pp)
     lobs <- -(log10(observed))
    
     expected <- c(1:length(observed))
     lexp <- -(log10(expected / (length(expected)+1)))
    
    
    
     #pdf("qqplot.pdf", width=6, height=6)
     plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
     points(lexp, lobs, pch=23, cex=.4, bg="black")
    
     *
     */
    //Algorithm proposed by Storey, J. D. and R. Tibshirani (2003). "Statistical significance for genomewide studies." Proc Natl Acad Sci U S A 100(16): 9440-5.
    public static double estimateNullHypothesisProportion(DoubleArrayList pValues) throws Exception {
        double step = 0.005;
        int pValueNum = pValues.size();
        int num = (int) (0.95 / step) + 1; // I modified this from 0.95 to be 1
        DoubleArrayList pi0List = new DoubleArrayList();
        DoubleArrayList lambdaList = new DoubleArrayList();
        int counter = 0;
        double currentValue = 0.0;
        double pi0;

        for (int i = 0; i < num; i++) {
            counter = 0;
            //note: the I2-value must be ordered
            while (counter < pValueNum && pValues.get(counter) <= currentValue) {
                counter++;
            }
            counter = pValueNum - counter;
            pi0 = (double) counter / (pValueNum * (1 - currentValue));
            pi0List.add(pi0);
            lambdaList.add(currentValue);

            currentValue += step;
        }

        pi0List.quickSort();
        boolean useMode = true;
        if (useMode) {
            //note the mode often works a little better than the median
            // but their performances are very close in my simulation
            double cutoffLen = 0.00125;
            Map<Double, Integer> modeCounter = new HashMap<Double, Integer>();
            for (int i = 0; i < num; i++) {
                double mid = pi0List.getQuick(i);
                int count = 0;
                int minIndex = i - 1;
                while (minIndex >= 0 && pi0List.getQuick(minIndex) >= (mid - cutoffLen)) {
                    minIndex--;
                }
                int maxIndex = i + 1;
                while (maxIndex < num && pi0List.getQuick(maxIndex) <= (mid + cutoffLen)) {
                    maxIndex++;
                }
                count = maxIndex - minIndex - 1;
                Integer goupNum = modeCounter.get(mid);
                if (goupNum == null) {
                    modeCounter.put(mid, count);
                } else {
                    modeCounter.put(mid, goupNum + count);
                }
            }

            double model = pi0List.getQuick(0);
            int maxCounter = modeCounter.get(model);
            for (Map.Entry<Double, Integer> m : modeCounter.entrySet()) {
                if (m.getValue() > maxCounter) {
                    model = m.getKey();
                    maxCounter = m.getValue();
                }
            }
            if (model > 1) {
                model = 1;
            }
            return model;
        } else {
            //note the median often works better than natural cubic spline
            // because the  pi0 becomes quite unstable when is close to currentValue
            double median = Descriptive.median(pi0List);
            if (median > 1) {
                median = 1;
            }
            return median;
        }

        /*        
         // System.out.println(median);
         num = lambdaList.size();
         if (num == 0) {
         return 0;
         } else if (num < 3) {
         return pi0List.get(num - 1);
         }
        
         double[] m0Array = new double[num];
         double[] lambdaArray = new double[num];
        
         for (int i = 0; i < num; i++) {
         m0Array[i] = pi0List.getQuick(i);
         lambdaArray[i] = lambdaList.getQuick(i);
         }
        
         UnivariateRealInterpolator interpolator = new SplineInterpolator();
         UnivariateRealFunction function = interpolator.interpolate(lambdaArray, m0Array);
         double interpolationX = lambdaArray[num - 1];
         double interpolatedY = function.value(interpolationX);
         //System.out.println(interpolatedY);
         return interpolatedY;
         */

        /*
         BufferedWriter bw = new BufferedWriter(new FileWriter("test.txt"));
         for (int i = 0; i < pValues.size(); i++) {
         bw.write(String.valueOf(pValues.getQuick(i)));
         bw.write("\n");
         }
         bw.close();
        
         pp <- scan("D:/home/mxli/MyJava/KGG2/test.txt", quiet= TRUE);
         hist(pp);
         *
         */
    }

    //note: try to make sampleSize smaller than subPopulationSize
    //for example (19061, 89, 307, 8) will not work, but (19061, 307, 89, 8) wiil although they have the same p value
    public static double hypergeometricEnrichmentTest(int populationSize, int subPopulationSize, int sampleSize,
            int subSampleSize) throws Exception {
        if (sampleSize == 0) {
            return 1;
        }
        double p = -9;
        if (subPopulationSize < sampleSize) {
            int tmp = sampleSize;
            sampleSize = subPopulationSize;
            subPopulationSize = tmp;
        }

        /*
         //calcluate pvalues       
         HyperGeometric hp = new HyperGeometric(populationSize, subPopulationSize, sampleSize, null);
         double I2 = 0;
         for (int k = 0; k < subSampleSize; k++) {
         I2 += hp.pdf(k);
         }
         p = 1 - I2;
        
        
         //I found the above colt function is more robust 
         //For example, this function does not work for  (19061, 307, 89, 8))
         HypergeometricDistribution hpd = new HypergeometricDistribution(populationSize, subPopulationSize, sampleSize);
        
         p = 1 - hpd.cumulativeProbability(0, subSampleSize - 1);
         if (p < 0) {
         p = 1.0E-40;
         }
         */
        //Later, I found this function is the most robust function among the three and the results are identical to R 1-phyper(7, 289,19061-289, 407)
        DiscreteDistributionInt hpd1 = new HypergeometricDist(subPopulationSize, populationSize, sampleSize);
        p = 1 - hpd1.cdf(subSampleSize - 1);

        if (Double.isNaN(p)) {
            p = -9;
        }
        return p;
    }
}
