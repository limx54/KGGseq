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
public class RefExon  implements Cloneable, Serializable{

    private int start;
    private int end;
    private int length;
    private int codingStart;
    private int codingEnd;
    private int codingLength;

    public RefExon(int start, int end, int length, int codingStart, int codingEnd, int codingLength) {
        this.start = start;
        this.end = end;
        this.length = length;
        this.codingStart = codingStart;
        this.codingEnd = codingEnd;
        this.codingLength = codingLength;
    }

    public RefExon(int start, int end, int length) {
        this.start = start;
        this.end = end;
        this.length = length;
        this.codingStart = 0;
        this.codingEnd = 0;
        this.codingLength = 0;
    }

    public int getCodingEnd() {
        return codingEnd;
    }

    public void setCodingEnd(int codingEnd) {
        this.codingEnd = codingEnd;
    }

    public int getCodingLength() {
        return codingLength;
    }

    public void setCodingLength(int codingLength) {
        this.codingLength = codingLength;
    }

    public int getCodingStart() {
        return codingStart;
    }

    public void setCodingStart(int codingStart) {
        this.codingStart = codingStart;
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
