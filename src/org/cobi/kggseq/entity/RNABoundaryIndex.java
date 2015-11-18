/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import cern.colt.list.IntArrayList;

/**
 *
 * @author mxli
 */
public class RNABoundaryIndex {

    protected int position;
    protected IntArrayList mRNAIndexList;
    
    public RNABoundaryIndex(int position) {
        this.position = position;
        mRNAIndexList = new IntArrayList(3);
    }

    public void addIndexes(int index) {
        mRNAIndexList.add(index);
    }
}
