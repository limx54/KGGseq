/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

/**
 *
 * @author mxli
 */
public class IntSet {

    public int index1;
    public int index2;
    public byte len;

    public IntSet(int index1, int index2, byte len) {
        this.index1 = index1;
        this.index2 = index2;
        this.len = len;
    }
}
