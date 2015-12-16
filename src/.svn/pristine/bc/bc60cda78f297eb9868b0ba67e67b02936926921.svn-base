/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

/**
 *
 * @author mxli
 */
public class RefCNV extends SeqSegment {

 
    String variantSubtype;
    int sampleSize = 0;
    int observedGains = 0;
    int observedLosses = 0;

    public RefCNV(String refID, int start, int end) {
        super(refID, start, end);
    }

    public RefCNV(String refID, int start, int end, String variantSubtype, int sampleSize, int observedGains, int observedLosses) {
        super(refID, start, end);
        this.variantSubtype = variantSubtype;
        this.sampleSize = sampleSize;
        this.observedGains = observedGains;
        this.observedLosses = observedLosses;
    }

    public int getObservedGains() {
        return observedGains;
    }

    public int getObservedLosses() {
        return observedLosses;
    }

    public int getSampleSize() {
        return sampleSize;
    }
    
    
}
