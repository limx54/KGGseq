/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.util.Comparator;

/**
 *
 * @author mxli
 */
public class GeneSetWilcoxPValueComparator implements Comparator<GeneSet> {

    @Override
    public int compare(final GeneSet arg0, final GeneSet arg1) { 
        return Double.compare(arg0.getWilcoxonPValue(), arg1.getWilcoxonPValue());
    }
}