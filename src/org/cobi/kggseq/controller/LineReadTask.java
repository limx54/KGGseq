/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import java.io.BufferedReader;
import java.util.concurrent.Callable;
import org.cobi.util.thread.Task;

/**
 *
 * @author mxli
 */
public class LineReadTask extends Task implements Callable<String> {

    String[] rows;
    BufferedReader br;
    int effectiveIndex;
    int maxRowNum;
    String currentLine;

    public LineReadTask(BufferedReader br, int maxRowNum) {
        this.rows = new String[maxRowNum];
        this.br = br;
        this.maxRowNum = maxRowNum;
    }

    public String[] getRows() {
        return rows;
    }

    public int getEffectiveIndex() {
        return effectiveIndex;
    }

    @Override
    public String call() throws Exception {
        // long startTime = System.currentTimeMillis();
//by default
        effectiveIndex = 0;
        while ((currentLine = br.readLine()) != null) {
            rows[effectiveIndex++] = currentLine;
            if (effectiveIndex >= maxRowNum) {
                break;
            }
        }

        // fireTaskComplete();
        //  String info = "Elapsed time: " + (System.currentTimeMillis() - startTime) / 1000 + " seconds.";
        //  System.out.println(info);
        //return info;
        return "";
        // 

    }

}
