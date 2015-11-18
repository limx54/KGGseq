/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

/**
 *
 * @author mxli
 */
public class ProteinDomain {

    String type;
    String description;
    String status;
    int start;
    int end;

    public ProteinDomain(String type, String description, String status, int start, int end) {
        this.type = type;
        this.description = description;
        this.status = status;
        this.start = start;
        this.end = end;
    }

    @Override
    public String toString(){
        return type+'|'+description+'|'+status+'|'+start+':'+end;
    }
    
    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public String getStatus() {
        return status;
    }

    public void setStatus(String status) {
        this.status = status;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }
    
    
}
