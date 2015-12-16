/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.util.Comparator;

/**
 *
 * @author MX Li
 */
public class SeqSegmentComparator implements Comparator<SeqSegment> {

    @Override
    public int compare(SeqSegment arg0, SeqSegment arg1) {
        int result = -1; 
        if (arg0.start == arg1.start) {
            result = arg0.end - arg1.end;
        } else {
            result = arg0.start - arg1.start;
        }
        return result; 
    }
}
