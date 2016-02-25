/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import java.util.concurrent.Callable;
import org.cobi.util.thread.Task;

/**
 *
 * @author mxli
 */
public class LDCalcTask extends Task implements Callable<String> {

    int totalPedSubjectNum;
    int startRowIndex;
    int endRowIndex;//not inclusive
    int startColIndex;
    int endColIndex;//not inclusive
    boolean isPhased;
    //for LD calculation 
    int[][] bits1;
    int[][] bits2;
    int[][] bits3;
    int[][] bits4;
    double[] mean1;
    int[] sum12;
    int totalVarNum;
    double[] rArray;

    int unitNum;
    boolean isSmallRectangle = true;
    boolean[] hasMissingGty;
    double adjIndivSize;

    public boolean isIsSmallRectangle() {
        return isSmallRectangle;
    }

    public void setIsSmallRectangle(boolean isSmallRectangle) {
        this.isSmallRectangle = isSmallRectangle;
    }

    public LDCalcTask(int[][] bits1, int[][] bits2, int[][] bits3, int[][] bits4, double[] mean1, int[] sum12, boolean[] hasMissingGty,
            int totalPedSubjectNum, int startRowIndex, int endRowIndex, int startColIndex, int endColIndex, double[] rArray, boolean isPhased) {
        this.totalPedSubjectNum = totalPedSubjectNum;
        this.startRowIndex = startRowIndex;
        this.endRowIndex = endRowIndex;
        this.startColIndex = startColIndex;
        this.endColIndex = endColIndex;
        this.rArray = rArray;
        this.isPhased = isPhased;

        // 
        this.bits1 = bits1;
        this.bits2 = bits2;
        this.bits3 = bits3;
        this.bits4 = bits4;

        this.mean1 = mean1;
        this.sum12 = sum12;
        // this.totalVarNum =varSize;
        // this.totalIndiv =indivSize;
        this.unitNum = this.bits1[0].length;
        this.hasMissingGty = hasMissingGty;

    }

    public void fillSmallRectangle() {
        double r = 0;
        int x;
        double douSize = totalPedSubjectNum * 2;
        totalVarNum = endColIndex - startColIndex;
        int countAB, tmpNum = 0;
        double tmpD1, tmpD2, tmpD3;
        int sum1i, sum1j, sum12i, sum12j;
        int nij = 0;
        int exsit;
        int tempInt;
        if (isPhased) {
            if (startRowIndex == startColIndex) {
                for (int i = startRowIndex; i < endRowIndex; i++) {
                    for (int j = i + 1; j < endColIndex; j++) {
                        countAB = 0;

                        //refecen http://blog.csdn.net/hitwhylz/article/details/10122617
                        for (int k = 0; k < unitNum; k++) {
                            x = bits1[i][k] & bits1[j][k];
                            countAB += (Integer.bitCount(x));  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
                            x = bits2[i][k] & bits2[j][k];
                            countAB += (Integer.bitCount(x));
                        }
                        if (hasMissingGty[i] || hasMissingGty[j]) {
                            sum1i = 0;
                            sum1j = 0;
                            nij = 0;

                            for (int k = 0; k < unitNum; k++) {
                                exsit = bits3[i][k] & bits3[j][k];
                                x = exsit;
                                nij += (Integer.bitCount(x));

                                x = bits1[i][k] & exsit;
                                tempInt = Integer.bitCount(x);
                                sum1i += tempInt;
                                x = bits2[i][k] & exsit;
                                tempInt = Integer.bitCount(x);
                                sum1i += tempInt;

                                x = bits1[j][k] & exsit;
                                tempInt = Integer.bitCount(x);
                                sum1j += tempInt;
                                x = bits2[j][k] & exsit;
                                tempInt = Integer.bitCount(x);
                                sum1j += tempInt;
                            }
                            nij = nij << 1;
                            tmpD1 = ((double) sum1i) * sum1j / nij;
                            tmpD1 = (countAB - tmpD1);
                            tmpD2 = ((double) sum1i) / nij;
                            tmpD3 = ((double) sum1j) / nij;
                            tmpD2 = tmpD2 * (1 - tmpD2) * tmpD3 * (1 - tmpD3);
                            r = tmpD1 / (nij * Math.sqrt(tmpD2));
                        } else {
                            tmpD1 = (countAB / douSize - mean1[i] * mean1[j]);
                            tmpD2 = mean1[i] * (1 - mean1[i]) * mean1[j] * (1 - mean1[j]);
                            r = tmpD1 / Math.sqrt(tmpD2);
                            if (Double.isNaN(r)) {
                                int ssss = 0;
                            }
                        }
                        //correction
                        // r = (douSize * r - 1) / (douSize - 3);
                        //johny's correction                
                        //r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1.0 + 2 * (1 - r) / (adjIndivSize - 3.3));
                        rArray[(i - startRowIndex) * totalVarNum + j - startColIndex] = r;
                        // System.out.println();
                    }
                }
            } else {
                for (int i = startRowIndex; i < endRowIndex; i++) {
                    for (int j = startColIndex; j < endColIndex; j++) {
                        countAB = 0;
                        //refecen http://blog.csdn.net/hitwhylz/article/details/10122617
                        for (int k = 0; k < unitNum; k++) {
                            x = bits1[i][k] & bits1[j][k];
                            countAB += (Integer.bitCount(x));  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
                            x = bits2[i][k] & bits2[j][k];
                            countAB += (Integer.bitCount(x));
                        }
                        //refecen http://blog.csdn.net/hitwhylz/article/details/10122617
                        if (hasMissingGty[i] || hasMissingGty[j]) {
                            sum1i = 0;
                            sum1j = 0;
                            nij = 0;

                            for (int k = 0; k < unitNum; k++) {
                                exsit = bits3[i][k] & bits3[j][k];
                                x = exsit;
                                nij += (Integer.bitCount(x));

                                x = bits1[i][k] & exsit;
                                tempInt = Integer.bitCount(x);
                                sum1i += tempInt;
                                x = bits2[i][k] & exsit;
                                tempInt = Integer.bitCount(x);
                                sum1i += tempInt;

                                x = bits1[j][k] & exsit;
                                tempInt = Integer.bitCount(x);
                                sum1j += tempInt;
                                x = bits2[j][k] & exsit;
                                tempInt = Integer.bitCount(x);
                                sum1j += tempInt;
                            }
                            nij = nij << 1;
                            tmpD1 = ((double) sum1i) * sum1j / nij;
                            tmpD1 = (countAB - tmpD1);
                            tmpD2 = ((double) sum1i) / nij;
                            tmpD3 = ((double) sum1j) / nij;
                            tmpD2 = tmpD2 * (1 - tmpD2) * tmpD3 * (1 - tmpD3);
                            r = tmpD1 / (nij * Math.sqrt(tmpD2));
                        } else {
                            tmpD1 = (countAB / douSize - mean1[i] * mean1[j]);
                            tmpD2 = mean1[i] * (1 - mean1[i]) * mean1[j] * (1 - mean1[j]);
                            r = tmpD1 / Math.sqrt(tmpD2);
                        }
                        //correction
                        // r = (douSize * r - 1) / (douSize - 3);
                        //johny's correction                
                        // r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1 + 2.0 * (1 - r) / (adjIndivSize - 3.3));
                        rArray[(i - startRowIndex) * totalVarNum + j - startColIndex] = r;
                    }
                }
            }
        } else if (startRowIndex == startColIndex) {
            for (int i = startRowIndex; i < endRowIndex; i++) {
                for (int j = i + 1; j < endColIndex; j++) {
                    countAB = 0;
                    for (int k = 0; k < unitNum; k++) {
                        x = bits1[i][k] & bits1[j][k];
                        countAB += (Integer.bitCount(x) << 1);
                        x = bits2[i][k] & bits2[j][k];
                        countAB += (Integer.bitCount(x) << 1);
                        x = bits3[i][k] & bits3[j][k];
                        countAB -= Integer.bitCount(x);
                    }

                    if (hasMissingGty[i] || hasMissingGty[j]) {
                        sum1i = 0;
                        sum1j = 0;
                        sum12i = 0;
                        sum12j = 0;
                        nij = 0;

                        for (int k = 0; k < unitNum; k++) {
                            exsit = bits4[i][k] & bits4[j][k];
                            x = exsit;
                            nij += (Integer.bitCount(x));  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 

                            x = bits1[i][k] & exsit;
                            tempInt = Integer.bitCount(x);
                            sum1i += (tempInt << 1);  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
                            sum12i += (tempInt << 2);  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
                            x = bits3[i][k] & exsit;
                            tempInt = Integer.bitCount(x);
                            sum1i += (tempInt);
                            sum12i += (tempInt);

                            x = bits1[j][k] & exsit;
                            tempInt = Integer.bitCount(x);
                            sum1j += (tempInt << 1);  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
                            sum12j += (tempInt << 2);
                            x = bits3[j][k] & exsit;
                            tempInt = Integer.bitCount(x);
                            sum1j += (tempInt);
                            sum12j += (tempInt);
                        }

                        countAB = nij * countAB - sum1i * sum1j;
                        sum12i = nij * sum12i - sum1i * sum1i;
                        sum12j = nij * sum12j - sum1j * sum1j;
                        r = (countAB) / Math.sqrt(sum12i * sum12j);
                        //r = (freqAB - mean1[i] * mean1[j]) / (sum12[i] * sum12[j]);
                    } else {
                        //Strange! this formula is even faster
                        r = (totalPedSubjectNum * countAB - mean1[i] * mean1[j]) / Math.sqrt((totalPedSubjectNum * sum12[i] - mean1[i] * mean1[i]) * (totalPedSubjectNum * sum12[j] - mean1[j] * mean1[j]));
                        //r = (freqAB - mean1[i] * mean1[j]) / (sum12[i] * sum12[j]);
                    }

                    //r = r * r;
                    //correction
                    // r = (douSize * r - 1) / (douSize - 3);
                    //johny's correction                
                    // r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1 + 2.0 * (1 - r) / (adjIndivSize - 3.3));
                     rArray[(i - startRowIndex) * totalVarNum + j - startColIndex] = r;
                   
                    //   System.out.print((i - startRowIndex) + "," + (j - startColIndex) + ":" + r + " ");
                }
                //  System.out.println();

            }
        } else {
            for (int i = startRowIndex; i < endRowIndex; i++) {
                for (int j = startColIndex; j < endColIndex; j++) {
                    countAB = 0;
                    for (int k = 0; k < unitNum; k++) {
                        x = bits1[i][k] & bits1[j][k];
                        countAB += (Integer.bitCount(x) << 1);
                        x = bits2[i][k] & bits2[j][k];
                        countAB += (Integer.bitCount(x) << 1);
                        x = bits3[i][k] & bits3[j][k];
                        countAB -= Integer.bitCount(x);
                    }

                    if (hasMissingGty[i] || hasMissingGty[j]) {
                        sum1i = 0;
                        sum1j = 0;
                        sum12i = 0;
                        sum12j = 0;
                        nij = 0;

                        for (int k = 0; k < unitNum; k++) {
                            exsit = bits4[i][k] & bits4[j][k];
                            x = exsit;
                            nij += (Integer.bitCount(x));  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 

                            x = bits1[i][k] & exsit;
                            tempInt = Integer.bitCount(x);
                            sum1i += (tempInt << 1);  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
                            sum12i += (tempInt << 2);  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
                            x = bits3[i][k] & exsit;
                            tempInt = Integer.bitCount(x);
                            sum1i += (tempInt);
                            sum12i += (tempInt);

                            x = bits1[j][k] & exsit;
                            tempInt = Integer.bitCount(x);
                            sum1j += (tempInt << 1);  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
                            sum12j += (tempInt << 2);
                            x = bits3[j][k] & exsit;
                            tempInt = Integer.bitCount(x);
                            sum1j += (tempInt);
                            sum12j += (tempInt);
                        }
                        countAB = nij * countAB - sum1i * sum1j;
                        sum12i = nij * sum12i - sum1i * sum1i;
                        sum12j = nij * sum12j - sum1j * sum1j;
                        r = (countAB) / Math.sqrt(sum12i * sum12j);
                    } else {
                        // the mean1 are acuually are sum for unphased genotypes 
                        r = (totalPedSubjectNum * countAB - mean1[i] * mean1[j]) / Math.sqrt((totalPedSubjectNum * sum12[i] - mean1[i] * mean1[i]) * (totalPedSubjectNum * sum12[j] - mean1[j] * mean1[j]));
                        //r = (freqAB - mean1[i] * mean1[j]) / (sum12[i] * sum12[j]);
                    }
                    // r = r * r;

                    //correction
                    // r = (douSize * r - 1) / (douSize - 3);
                    //johny's correction                
                    // r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1 + 2.0 * (1 - r) / (adjIndivSize - 3.3));
                      rArray[(i - startRowIndex) * totalVarNum + j - startColIndex] = r;
                   
                }
            }
        }
    }

    @Override
    public String call() throws Exception {
        long startTime = System.currentTimeMillis();
        fillSmallRectangle();
        fireTaskComplete();
        String info = "  LD vas estimated on block [" + startRowIndex + " , " + startRowIndex + "). Elapsed time: " + (System.currentTimeMillis() - startTime) / 1000 + " seconds.";
        //  System.out.println(info);
        //return info;
        return info;
        // 

    }
}
