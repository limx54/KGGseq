/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.util.ArrayList;
import java.util.List;

/**
 * this class is used to store the gene-based association analysis
 *
 * @author mxli
 */
public class Gene extends SeqSegment {

    public String geneSymb;
    public List<Double> testValues;
    public List<String> featureValues;
    public List<String> transcripts;

    public Gene(String refID, int start, int end) {
        super(start, end);
        geneSymb = refID;
        featureValues = new ArrayList<String>();
        testValues = new ArrayList<Double>();
        transcripts = new ArrayList<String>();
    }

    public Gene(String refID) {
        geneSymb = refID;
        featureValues = new ArrayList<String>();
        testValues = new ArrayList<Double>();
        transcripts = new ArrayList<String>();
    }

    public void addFeatureValue(String val) {
        featureValues.add(val);
    }

    public void addFeatureValueAt(int pos, String val) {
        if (pos >= featureValues.size()) {
            for (int i = featureValues.size(); i <= pos; i++) {
                featureValues.add(null);
            }
        }
        featureValues.set(pos, val);
    }

    public void addTestValue(double val) {
        testValues.add(val);
    }

    public void addTestValueAt(int pos, double val) {
        if (pos >= testValues.size()) {
            for (int i = testValues.size(); i <= pos; i++) {
                testValues.add(null);
            }
        }
        testValues.set(pos, val);
    }

    public String getSymbol() {
        return geneSymb;
    }

    public void setGeneSymb(String stringID) {
        this.geneSymb = stringID;
    }
}
