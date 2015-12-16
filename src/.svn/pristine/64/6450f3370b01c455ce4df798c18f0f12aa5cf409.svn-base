/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author MX Li
 */
public class RefRefmRNA implements Cloneable, Serializable {

    private String refID;
    private int length;
    private String proteinRefID;
    private List<RefExon> exons;
    private List<RefIntron> introns;

    public RefRefmRNA(String refID) {
        this.refID = refID;
    }

    public void writeObject(ObjectOutputStream outputStream) throws IOException {
        outputStream.defaultWriteObject();
        //outputStream.writeObject(exons);
       // outputStream.writeObject(introns);
    }

    private void readObject(ObjectInputStream inputStream) throws IOException, ClassNotFoundException {
        inputStream.defaultReadObject();
        //exons = (List<Exon>) inputStream.readObject();
       // introns = (List<Intron>) inputStream.readObject();
    }

    public RefRefmRNA(String refID, int length, String proteinRefID) {
        this.refID = refID;
        this.length = length;
        this.proteinRefID = proteinRefID;
        this.exons = new ArrayList<RefExon>();
        this.introns = new ArrayList<RefIntron>();
    }

    public List<RefExon> getExons() {
        return exons;
    }

    public void setExons(List<RefExon> exons) {
        this.exons = exons;
    }

    public void addExon(RefExon exon) {
        exons.add(exon);
    }

    public int getExonsLength() {
        int len = 0;
        for (int i = 0; i < exons.size(); i++) {
            len += exons.get(i).getLength();
        }
        return len;
    }

    public int getIntronsLength() {
        int len = 0;
        for (int i = 0; i < introns.size(); i++) {
            len += introns.get(i).getLength();
        }
        return len;
    }

    public int getExonsNum() {
        return exons.size();
    }

    public int getIntronsNum() {
        return introns.size();
    }

    public List<RefIntron> getIntrons() {
        return introns;
    }

    public void setIntrons(List<RefIntron> introns) {
        this.introns = introns;
    }

    public void addIntron(RefIntron intron) {
        this.introns.add(intron);
    }

    public int getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public String getProteinRefID() {
        return proteinRefID;
    }

    public void setProteinRefID(String proteinRefID) {
        this.proteinRefID = proteinRefID;
    }

    public String getRefID() {
        return refID;
    }

    public void setRefID(String refID) {
        this.refID = refID;
    }
}
