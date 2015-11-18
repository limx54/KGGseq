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
public class CombOrderComparator implements Comparator<CombOrders> {

        @Override
        public int compare(CombOrders arg0, CombOrders arg1) {
            return Double.compare(arg0.auc, arg1.auc);
        }
    }
