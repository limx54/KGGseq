/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

import java.util.Comparator;

/**
 *
 * @author mxli
 */
public class StringArrayIntegerComparator implements Comparator<String[]> {

    private int index = 0;

    public StringArrayIntegerComparator(int ind) {
        this.index = ind;
    }

    @Override
    public int compare(String[] o1, String[] o2) {
        return Util.parseInt(o1[index]) - Util.parseInt(o2[index]);
    }
}