/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

import java.util.Comparator;

/**
 *
 * @author mxli
 */
public class StringArrayDoubleComparator implements Comparator<String[]> {

    private int index = 0;

    public StringArrayDoubleComparator(int ind) {
        this.index = ind;
    }

    @Override
    public int compare(String[] o1, String[] o2) {
        if (o1[index] == null || o2[index] == null || o1[index].equals(".") || o2[index].equals(".")) {
            return 0;
        }
        return Double.compare(Double.parseDouble(o1[index]), Double.parseDouble(o2[index]));
    }
}
