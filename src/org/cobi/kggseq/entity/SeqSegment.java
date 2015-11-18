/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

/**
 *
 * @author mxli
 */
public class SeqSegment {

    protected int start = -1;
    protected int end = -1;
    protected String description = null;

    public SeqSegment(int start, int end) {
        this.start = start;
        this.end = end;
    }
    public SeqSegment() {
      
    }

    public SeqSegment(String des, int start, int end) {
        this.description = des;
        this.start = start;
        this.end = end;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String feature) {
        this.description = feature;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }
}
