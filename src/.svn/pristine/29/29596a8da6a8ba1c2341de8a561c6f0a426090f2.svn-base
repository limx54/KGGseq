/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.thread;

import java.util.Arrays;

/**
 *
 * @author mxli
 */
public class Util {

    //the block number cannot exceed threadNum

    static public int[] partitionEvenBlock(int threadNum, int startIndex, int endIndex) throws Exception {
        int totalSnpSize = endIndex - startIndex;
        int intervalLen = totalSnpSize / threadNum;
        int[] bigBlockIndexes = null;

        if (intervalLen == 0) {
            //no need to block
            bigBlockIndexes = new int[2];
            bigBlockIndexes[0] = startIndex;
            bigBlockIndexes[1] = endIndex;

        } else {
            bigBlockIndexes = new int[threadNum + 1];
            Arrays.fill(bigBlockIndexes, startIndex);
            for (int s = 1; s < threadNum; s++) {
                bigBlockIndexes[s] = startIndex + s * intervalLen;
            }
            bigBlockIndexes[threadNum] = endIndex;
        }

        return bigBlockIndexes;
    }

}
