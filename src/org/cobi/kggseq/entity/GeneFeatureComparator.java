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
public class GeneFeatureComparator implements Comparator<GeneFeature> {

    @Override
    public int compare(GeneFeature arg0, GeneFeature arg1) {
        return arg0.id - arg1.id;
    }
}