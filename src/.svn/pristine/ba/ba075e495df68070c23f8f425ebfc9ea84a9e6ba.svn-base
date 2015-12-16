/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.util.ArrayList;
import java.util.List;
import org.cobi.kggseq.Constants;

/**
 *
 * @author mxli
 */
public class Individual implements Constants {

    private static final long serialVersionUID = 2L;
    private String familyID;
    private String individualID;
    private String momID;
    private String dadID;
    private int gender;
    private int affectedStatus;
    private int liability; //optional
    private boolean hasGenotypes=false;
    private double[] traitValues;
    private String labelInChip;

    public Individual() {
       
    }

    public boolean isHasGenotypes() {
        return hasGenotypes;
    }

    public void setHasGenotypes(boolean hasGenotypes) {
        this.hasGenotypes = hasGenotypes;
    }

    public Individual(String individualID) {
        this.individualID = individualID;     
    }

    /*
     @Override
     public Object clone() throws CloneNotSupportedException {
     Individual o = null;
     o = (Individual) super.clone();
     o.traitValues = new ArrayList<String>();
     o.traitValues.addAll(this.traitValues);
    
     return o;
     }
     * 
     */
    public String getLabelInChip() {
        return labelInChip;
    }

    public void setLabelInChip(String labelInChip) {
        this.labelInChip = labelInChip;
    }

    /**
     * Get the value of traitValues
     *
     * @return the value of traitValues
     */
    public double[] getTraits() {
        return traitValues;
    }

    /**
     * Set the value of traitValues
     *
     * @param traits 
     */
    public void setTraits(double[] traits) {
        this.traitValues = traits;
    }
 
    public int getAffectedStatus() {
        return affectedStatus;
    }

    public void setAffectedStatus(int affectedStatus) {
        this.affectedStatus = affectedStatus;
    }

    public String getDadID() {
        return dadID;
    }

    public void setDadID(String dadID) {
        this.dadID = dadID;
    }

    public String getFamilyID() {
        return familyID;
    }

    public void setFamilyID(String familyID) {
        this.familyID = familyID;
    }

    public int getGender() {
        return gender;
    }

    public void setGender(int gender) {
        this.gender = gender;
    }

    public String getIndividualID() {
        return individualID;
    }

    public void setIndividualID(String individualID) {
        this.individualID = individualID;
    }

    public int getLiability() {
        return liability;
    }

    public void setLiability(int liability) {
        this.liability = liability;
    }

    public String getMomID() {
        return momID;
    }

    public void setMomID(String momID) {
        this.momID = momID;
    }
}
