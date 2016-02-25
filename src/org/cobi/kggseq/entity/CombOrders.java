/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.util.HashSet;
import java.util.Set;
import org.cobi.randomforests.MyRandomForest;
import org.cobi.util.text.Util;

/**
 *
 * @author mxli
 */
public class CombOrders {

    public  Set<Integer> indexes = null;
    public double auc = 0;
    public RegressionParams rp;
    public MyRandomForest rf;

    public CombOrders(String index, double v, RegressionParams r) {
        indexes = new HashSet<Integer>();
        String[] indexesStr = index.split("-");
        for (int i = 0; i < indexesStr.length; i++) {
            indexes.add(Util.parseInt(indexesStr[i]));
        }
        auc = v;
        rp = r;
    }

    public CombOrders(String index, double v, MyRandomForest r) {
        indexes = new HashSet<Integer>();
        String[] indexesStr = index.split("-");
        for (int i = 0; i < indexesStr.length; i++) {
            indexes.add(Util.parseInt(indexesStr[i]));
        }
        auc = v;
        rf = r;
    }
}
