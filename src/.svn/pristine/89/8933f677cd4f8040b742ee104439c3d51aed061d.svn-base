// (c) 2009-2011 Miaoxin Li
// This file is distributed as part of the KGG source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
// Permission is granted for you to use this file to compile IGG.
// All computer programs have bugs. Use this file at your own risk.
// Tuesday, March 01, 2011
package org.cobi.kggseq.entity;

import java.io.Serializable;

/**
 *
 * @author MX Li
 */
public class PPIEdge implements Serializable {

    private static final long serialVersionUID = 3L; 
    int id;
    double score;
   

    
    public PPIEdge(int id, double score) {
        this.id = id;
        this.score = score;
    }

    
    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    @Override
    public String toString() { // Always good for debugging
        return "E" + id;
    }
}
