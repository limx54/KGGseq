/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.util.Comparator;

/**
 *
 * @author mxli
 */
public class DoubleArrayListComparator implements Comparator<double[]> {

    int index;

    public DoubleArrayListComparator(int index) {
        this.index = index;
    }
    
    @Override
    public int compare(final double[] arg0, final double[] arg1) {
        return Double.compare(arg0[index], arg1[index]);
    }
}
