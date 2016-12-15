/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import cern.colt.list.IntArrayList;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author mxli
 */
public class FiltrationSummarySet {

    String name;
    IntArrayList counts;
    List<String> messages;
    int availableFeatureIndex = -1;

    public FiltrationSummarySet(String name, int availableFeatureIndex) {
        this.name = name;
        counts = new IntArrayList();
        messages = new ArrayList<String>();
        this.availableFeatureIndex = availableFeatureIndex;
    }

    public int getAvailableFeatureIndex() {
        return availableFeatureIndex;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void initiateAMessage(int count, String msg) {
        counts.add(count);
        messages.add(msg);
    }

    public void increaseCount(int index, int count) {
        int t = counts.getQuick(index) + count;
        counts.setQuick(index, t);
    }

    @Override
    public String toString() {
        int size = counts.size();
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < size; i++) {
            if (counts.getQuick(i) == 0) {
                continue;
            }
            sb.append(counts.getQuick(i));
            sb.append(" ");
            sb.append(messages.get(i));
            sb.append("\n");
        }
        if (sb.length() > 0) {
            return sb.substring(0, sb.length() - 1);
        } else {
            return "";
        }
    }

    public String toString(int startIndex, int endIndex, String delimeter) {
        int size = counts.size();
        if (size < endIndex) {
            endIndex = size;
        }
        StringBuilder sb = new StringBuilder();
        for (int i = startIndex; i < endIndex; i++) {
            if (counts.getQuick(i) == 0) {
                continue;
            }
            sb.append(counts.getQuick(i));
            sb.append(" ");
            sb.append(messages.get(i));
            sb.append(delimeter);
        }
        return sb.toString();
    }
}
