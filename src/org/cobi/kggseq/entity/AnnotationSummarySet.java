/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.io.BufferedReader;

/**
 *
 * @author mxli
 */
public class AnnotationSummarySet {

    String name;
    BufferedReader br;
    StringBuilder lastLine;
    int availableFeatureIndex = -1;
    int annotNum;
    int leftNum;
    int totalNum;

    public AnnotationSummarySet(String name, BufferedReader br, StringBuilder lastLine, int annotNum, int leftNum, int totalNum, int availableFeatureIndex) {
        this.name = name;
        this.br = br;
        this.lastLine = lastLine;
        this.annotNum = annotNum;
        this.leftNum = leftNum;
        this.totalNum = totalNum;
        this.availableFeatureIndex = availableFeatureIndex;
    }

    public int getAvailableFeatureIndex() {
        return availableFeatureIndex;
    }

    public void setAvailableFeatureIndex(int availableFeatureIndex) {
        this.availableFeatureIndex = availableFeatureIndex;
    }

    public int getTotalNum() {
        return totalNum;
    }

    public void setTotalNum(int totalNum) {
        this.totalNum = totalNum;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public BufferedReader getBr() {
        return br;
    }

    public void setBr(BufferedReader br) {
        this.br = br;
    }

    public StringBuilder getLastLine() {
        return lastLine;
    }

    public void setLastLine(StringBuilder lastLine) {
        this.lastLine = lastLine;
    }

    public int getAnnotNum() {
        return annotNum;
    }

    public void setAnnotNum(int annotNum) {
        this.annotNum = annotNum;
    }

    public int getLeftNum() {
        return leftNum;
    }

    public void setLeftNum(int leftNum) {
        this.leftNum = leftNum;
    }

}
