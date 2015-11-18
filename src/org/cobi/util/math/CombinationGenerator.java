// (c) 2009-2011 Miaoxin Li
// This file is distributed as part of the KGG source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
// Permission is granted for you to use this file to compile IGG.
// All computer programs have bugs. Use this file at your own risk.
// Tuesday, March 01, 2011
package org.cobi.util.math;

import java.math.BigInteger;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * copied from http://www.merriampark.com/comb.htm
 */
public class CombinationGenerator {

    private int[] a;
    private int n;
    private int r;
    private BigInteger numLeft;
    private BigInteger total;

    //------------
    // Constructor
    //------------
    public CombinationGenerator(int n, int r) {
        if (r > n) {
            throw new IllegalArgumentException();
        }
        if (n < 1) {
            throw new IllegalArgumentException();
        }
        this.n = n;
        this.r = r;
        a = new int[r];
        BigInteger nFact = getFactorial(n);
        BigInteger rFact = getFactorial(r);
        BigInteger nminusrFact = getFactorial(n - r);
        total = nFact.divide(rFact.multiply(nminusrFact));
        reset();
    }

    //------
    // Reset
    //------
    public void reset() {
        for (int i = 0; i < a.length; i++) {
            a[i] = i;
        }
        numLeft = new BigInteger(total.toString());
    }

    //------------------------------------------------
    // Return number of combinations not yet generated
    //------------------------------------------------
    public BigInteger getNumLeft() {
        return numLeft;
    }

    //-----------------------------
    // Are there more combinations?
    //-----------------------------
    public boolean hasMore() {
        return numLeft.compareTo(BigInteger.ZERO) == 1;
    }

    //------------------------------------
    // Return total number of combinations
    //------------------------------------
    public BigInteger getTotal() {
        return total;
    }

    //------------------
    // Compute factorial
    //------------------
    private static BigInteger getFactorial(int n) {
        BigInteger fact = BigInteger.ONE;
        for (int i = n; i > 1; i--) {
            fact = fact.multiply(new BigInteger(Integer.toString(i)));
        }
        return fact;
    }

    //--------------------------------------------------------
    // Generate next combination (algorithm from Rosen p. 286)
    //--------------------------------------------------------
    public int[] getNext() {

        if (numLeft.equals(total)) {
            numLeft = numLeft.subtract(BigInteger.ONE);
            return a;
        }

        int i = r - 1;
        while (a[i] == n - r + i) {
            i--;
        }
        a[i] = a[i] + 1;
        for (int j = i + 1; j < r; j++) {
            a[j] = a[i] + j - i;
        }

        numLeft = numLeft.subtract(BigInteger.ONE);
        return a;

    }

    public static boolean[] toBinaryArray(int number, int base) {
        final boolean[] ret = new boolean[base];
        for (int i = 0; i < base; i++) {
            ret[base - 1 - i] = (1 << i & number) != 0;
        }
        return ret;
    }

    public static String toBinaryString(int n, int base) {
        StringBuilder binary = new StringBuilder(Integer.toBinaryString(n));
        int bit = binary.length();
        StringBuilder zeroes = new StringBuilder();
        for (int i = bit; i < base; i++) {
            zeroes.append('0');
        }
        return zeroes.append(binary).toString();
    }

    public static void main(String[] args) {
        String[] elements = {"a", "b", "c", "d", "e", "f", "g"};
        int[] indices = null;
        CombinationGenerator x = new CombinationGenerator(elements.length, 2);
        StringBuffer combination;
        while (x.hasMore()) {
            combination = new StringBuffer();
            indices = x.getNext();
            for (int i = 0; i < indices.length; i++) {
                combination.append(elements[indices[i]]);
            }
            //   System.out.println(combination.toString());
        }
        Map<String, boolean[]> unphasedGtyCodingMap = new HashMap<String, boolean[]>();
        Map<String, String> codingUnphasedGtyCodingMap = new HashMap<String, String>();
        for (int alleleNum = 2; alleleNum <= 10; alleleNum++) {
            String[] gtys = new String[alleleNum * (alleleNum + 1) / 2];
            String bits = Integer.toBinaryString(gtys.length);
            int base = bits.length();
            x = new CombinationGenerator(alleleNum, 2);
            int s = 0;

            while (x.hasMore()) {
                indices = x.getNext();
                if (indices[1] - indices[0] == 1) {
                    gtys[s] = indices[0] + "/" + indices[0];
                    s++;
                }
                gtys[s] = indices[0] + "/" + indices[1];
                s++;
            }
            gtys[s] = indices[1] + "/" + indices[1];
            s = 0;
            for (s = 0; s < gtys.length; s++) {
                String str = gtys[s];
                String bitS = toBinaryString(s + 1, base);

                unphasedGtyCodingMap.put(str + ":" + alleleNum, toBinaryArray(s + 1, base));
                codingUnphasedGtyCodingMap.put(bitS + ":" + alleleNum, str);

                System.out.println(str + "\t" + bitS);
            }
        }
        String str = codingUnphasedGtyCodingMap.get("01:2");
        System.out.println(str);

        Map<String, boolean[]> phasedGtyCodingMap = new HashMap<String, boolean[]>();
        Map<String, String> codingPhasedGtyCodingMap = new HashMap<String, String>();
        for (int alleleNum = 2; alleleNum <= 10; alleleNum++) {
            String[] gtys = new String[alleleNum * (alleleNum)];
            String bits = Integer.toBinaryString(gtys.length);
            int base = bits.length();
            int s = 0;
            for (int i = 0; i < alleleNum; i++) {
                for (int j = 0; j < alleleNum; j++) {
                    gtys[s] = i + "|" + j;
                    s++;
                }
            }
            s = 0;
            for (s = 0; s < gtys.length; s++) {
                str = gtys[s];
                String bitS = toBinaryString(s + 1, base);

                phasedGtyCodingMap.put(str + ":" + alleleNum, toBinaryArray(s + 1, base));
                codingPhasedGtyCodingMap.put(bitS + ":" + alleleNum, str);

                System.out.println(str + "\t" + bitS);
            }
        }
        str = codingPhasedGtyCodingMap.get("100:2");
        System.out.println(str);
    }

}
