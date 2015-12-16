/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

import java.util.Comparator;

/**
 *
 * @author Miaoxin Li
 */
public class StringArrayStringComparator implements Comparator<String[]> {

    private int index = 0;

    public StringArrayStringComparator(int ind) {
        this.index = ind;
    }

    @Override
    public int compare(String[] o1, String[] o2) {
        return o1[index].compareTo(o2[index]);
    }
}