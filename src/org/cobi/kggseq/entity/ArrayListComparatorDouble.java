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
public class ArrayListComparatorDouble implements Comparator<String[]> {

    int index;

    public ArrayListComparatorDouble(int index) {
        this.index = index;
    }
    
    

    @Override
    public int compare(String[] arg0, String[] arg1) {
        return Double.compare(Double.parseDouble(arg0[index]), Double.parseDouble(arg1[index]));
    }
}
