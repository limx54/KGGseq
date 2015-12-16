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
public class StringArrayStringFinder implements Comparator<Object> {

    private int index = 0;

    public StringArrayStringFinder(int ind) {
        this.index = ind;
    }

    public int compare(Object o1, Object o2) {
        String[] str1 = (String[]) o1;
        String str2 = (String) o2;
        if ((str1[index] == null) || (str2 == null)) {
            return -1;
        } else {
            return str1[index].compareTo(str2);
        }
    }
}