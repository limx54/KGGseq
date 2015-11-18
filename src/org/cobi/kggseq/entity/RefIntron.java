/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.io.Serializable;

/**
 *
 * @author MX Li
 */
public class RefIntron  implements Cloneable, Serializable{

    private int start;
    private int end;
    private int length;

    public RefIntron(int start, int end, int length) {
        this.start = start;
        this.end = end;
        this.length = length;
    }


    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

}
