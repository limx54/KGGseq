/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.math;

/**
 *
 * @author mxli
 */
public class ArrayQuickSort {

    private int keyIndex = 0;
    private int array[];
    private int arrays[][];
    private int length;

    public void sort(int[] inputArr) {

        if (inputArr == null || inputArr.length == 0) {
            return;
        }
        this.array = inputArr;
        length = inputArr.length;
        quickSort(0, length - 1);
    }

    public void sort(int[][] inputArr, int index) {
        keyIndex = index;
        if (inputArr == null || inputArr.length == 0) {
            return;
        }
        this.arrays = inputArr;
        length = inputArr.length;
        quickSort2D(0, length - 1);
    }

    private void quickSort(int lowerIndex, int higherIndex) {
        int i = lowerIndex;
        int j = higherIndex;
        // calculate pivot number, I am taking pivot as middle index number
        int pivot = array[lowerIndex + (higherIndex - lowerIndex) / 2];
        // Divide into two arrays
        while (i <= j) {
            /**
             * In each iteration, we will identify a number from left side which
             * is greater then the pivot value, and also we will identify a
             * number from right side which is less then the pivot value. Once
             * the search is done, then we exchange both numbers.
             */
            while (array[i] < pivot) {
                i++;
            }
            while (array[j] > pivot) {
                j--;
            }
            if (i <= j) {
                exchangeNumbers(i, j);
                //move index to next position on both sides
                i++;
                j--;
            }
        }
        // call quickSort() method recursively
        if (lowerIndex < j) {
            quickSort(lowerIndex, j);
        }
        if (i < higherIndex) {
            quickSort(i, higherIndex);
        }
    }

    public void quickSort2D(int lowerIndex, int higherIndex) {
        int i = lowerIndex;
        int j = higherIndex;
        int[] temp;
        // calculate pivot number, I am taking pivot as middle index number
        int pivot = arrays[lowerIndex + (higherIndex - lowerIndex) / 2][keyIndex];
        // Divide into two arrays
        while (i <= j) {
            /**
             * In each iteration, we will identify a number from left side which
             * is greater then the pivot value, and also we will identify a
             * number from right side which is less then the pivot value. Once
             * the search is done, then we exchange both numbers.
             */
            while (arrays[i][keyIndex] < pivot) {
                i++;
            }
            while (arrays[j][keyIndex] > pivot) {
                j--;
            }
            if (i <= j) {
                temp = arrays[i];
                arrays[i] = arrays[j];
                arrays[j] = temp;
                //move index to next position on both sides
                i++;
                j--;
            }
        }
        // call quickSort() method recursively
        if (lowerIndex < j) {
            quickSort2D(lowerIndex, j);
        }
        if (i < higherIndex) {
            quickSort2D(i, higherIndex);
        }
    }

    private void exchangeNumbers(int i, int j) {
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }

    public static void main(String a[]) {

        ArrayQuickSort sorter = new ArrayQuickSort();
        int[] input = {24, 2, 45, 20, 56, 75, 2, 56, 99, 53, 12};
        sorter.sort(input);
        for (int i : input) {
            System.out.print(i);
            System.out.print(" ");
        }
    }
    // - See more at   : http://java2novice.com/java-sorting-algorithms/quick-sort/#sthash.691c9Fka.dpuf
}
