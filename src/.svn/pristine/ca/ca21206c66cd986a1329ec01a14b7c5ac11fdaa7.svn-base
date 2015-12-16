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
public class IntListComparator implements Comparator<int[]> {

    int index;

    public IntListComparator(int index) {
        this.index = index;
    }

    @Override
    public int compare(final int[] arg0, final int[] arg1) {
        return Integer.compare(arg0[index], arg1[index]);
    }
}
