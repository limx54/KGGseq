/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.coding;

import cern.colt.list.ObjectArrayList;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author mxli
 */
public class LoopPerformance {

    static int size = 10000000;
    private static List<String> list = new ArrayList<String>(size);
    private static String[] arry = new String[size];
    private static ObjectArrayList oal = new ObjectArrayList(size);

    static {
        for (int i = 0; i < size; i++) {
            list.add(String.valueOf(i));
            arry[i] = String.valueOf(i);
            oal.add(String.valueOf(i));
        }
    }

    public void testLoopPerformanceColt() {
        long startTime;
        long endTime;
        //style 1  
        String s = null;
        startTime = System.nanoTime();

        //style 2  
        startTime = System.nanoTime();
        for (int j = 0; j < oal.size(); j++) {
            s = (String) oal.getQuick(j);
            //  
        }
        endTime = System.nanoTime();
        System.out.println("Using collection.size() :: " + (endTime - startTime) / 1000000.0 + " ms");

        //style 3  
        startTime = System.nanoTime();
        int size = oal.size();
        for (int j = 0; j < size; j++) {
            s = (String) oal.getQuick(j);
        }
        endTime = System.nanoTime();
        System.out.println("Using [int size = list.size(); int j = 0; j < size ; j++] :: " + (endTime - startTime) / 1000000.0 + " ms");

        //style 4  
        startTime = System.nanoTime();

        for (int j = oal.size() - 1; j > 0; j--) {
            s = (String) oal.getQuick(j);
            //System.out.println(j);  
        }
        endTime = System.nanoTime();
        System.out.println("Using [int j = list.size(); j > 0 ; j--] :: " + (endTime - startTime) / 1000000.0 + " ms");
    }

    public void testLoopPerformance() {
        long startTime;
        long endTime;
        //style 1  
        String s = null;
        startTime = System.nanoTime();
        for (String i : list) {
            s = i;
//  
        }
        endTime = System.nanoTime();
        System.out.println("For each loop :: " + (endTime - startTime) / 1000000.0 + " ms");

        //style 2  
        startTime = System.nanoTime();
        for (int j = 0; j < list.size(); j++) {
            s = list.get(j);
            //  
        }
        endTime = System.nanoTime();
        System.out.println("Using collection.size() :: " + (endTime - startTime) / 1000000.0 + " ms");

        //style 3  
        startTime = System.nanoTime();
        int size = list.size();
        for (int j = 0; j < size; j++) {
            s = list.get(j);
        }
        endTime = System.nanoTime();
        System.out.println("Using [int size = list.size(); int j = 0; j < size ; j++] :: " + (endTime - startTime) / 1000000.0 + " ms");

        //style 4  
        startTime = System.nanoTime();

        for (int j = list.size() - 1; j > 0; j--) {
            s = list.get(j);
            //System.out.println(j);  
        }
        endTime = System.nanoTime();
        System.out.println("Using [int j = list.size(); j > 0 ; j--] :: " + (endTime - startTime) / 1000000.0 + " ms");
    }

    public void testLoopPerformanceArray() {
        long startTime;
        long endTime;
        //style 1  
        String s = null;
        startTime = System.nanoTime();
        for (String i : arry) {
            //   s = i;
//  
        }
        endTime = System.nanoTime();
        System.out.println("For each loop :: " + (endTime - startTime) + " ms");

        //style 2  
        startTime = System.nanoTime();
        for (int j = 0; j < arry.length; j++) {
            s = arry[j];
            //  
        }
        endTime = System.nanoTime();
        System.out.println("Using collection.size() :: " + (endTime - startTime) + " ms");

        //style 3  
        startTime = System.nanoTime();
        int size = arry.length;
        for (int j = 0; j < size; j++) {
            s = arry[j];
        }
        endTime = System.nanoTime();
        System.out.println("Using [int size = list.size(); int j = 0; j < size ; j++] :: " + (endTime - startTime) + " ms");

        //style 4  
        startTime = System.nanoTime();

        for (int j = arry.length - 1; j > 0; j--) {
            s = arry[j];
            //System.out.println(j);  
        }
        endTime = System.nanoTime();
        System.out.println("Using [int j = list.size(); j > 0 ; j--] :: " + (endTime - startTime) + " ms");
    }

    public static void main(String[] args) {
        new LoopPerformance().testLoopPerformance();
        System.out.println("Array");
        new LoopPerformance().testLoopPerformanceArray();
        System.out.println("Colt");
        new LoopPerformance().testLoopPerformanceColt();

    }
}
