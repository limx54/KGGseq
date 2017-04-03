/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import cern.colt.list.IntArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import org.cobi.kggseq.entity.AnnotationSummarySet;
import org.cobi.kggseq.entity.Chromosome;
import org.cobi.kggseq.entity.FiltrationSummarySet;
import org.cobi.kggseq.entity.Variant;
import org.cobi.util.text.BZPartReader;
import org.cobi.util.text.Util;
import org.cobi.util.thread.Task;

/**
 *
 * @author mxli
 */
public class VarAnnotTask extends Task implements Callable<String> {

    BZPartReader br;
    List<Variant> variantList;
    Chromosome chromosome;
    final AnnotationSummarySet ass;
    final FiltrationSummarySet fss;
    int indexCHROM = -1;
    int indexPOS;
    int[] predicIndex;
    int[] scoreIndexes;
    int[] altFreqColIndexes;

    int varAssignedScore = 0;
    int annotationType;

    public void setAltFreqColIndexes(int[] altFreqColIndexes) {
        this.altFreqColIndexes = altFreqColIndexes;
    }

    public VarAnnotTask(BZPartReader br, List<Variant> variantList, int indexCHROM, int indexPOS, int[] predicIndex, int[] scoreIndexes, final FiltrationSummarySet fss, int annotationType) {
        this.br = br;
        this.variantList = variantList;
        this.indexCHROM = indexCHROM;
        this.indexPOS = indexPOS;
        this.predicIndex = predicIndex;
        this.scoreIndexes = scoreIndexes;
        this.annotationType = annotationType;
        this.fss = fss;
        ass = null;
    }

    public VarAnnotTask(BZPartReader br, int indexCHROM, int indexPOS, int[] scoreIndexes) {
        this.br = br;
        this.indexCHROM = indexCHROM;
        this.indexPOS = indexPOS;
        this.scoreIndexes = scoreIndexes;
        ass = null;
        fss = null;
    }

    //assume the resource data have no identical coordinates
    public VarAnnotTask(BZPartReader br, List<Variant> variantList, int indexCHROM, int indexPOS, final AnnotationSummarySet ass, int annotationType) {
        this.br = br;
        this.variantList = variantList;
        this.indexCHROM = indexCHROM;
        this.indexPOS = indexPOS;
        this.annotationType = annotationType;
        this.ass = ass;
        fss = null;
    }

    public VarAnnotTask(BZPartReader br, int indexCHROM, int indexPOS, Chromosome chromosome, final AnnotationSummarySet ass, int annotationType) {
        this.br = br;
        this.variantList = chromosome.variantList;
        this.indexCHROM = indexCHROM;
        this.indexPOS = indexPOS;
        this.chromosome = chromosome;
        this.ass = ass;
        this.annotationType = annotationType;
        fss = null;
    }

    public void dbNSFPAnnot() throws Exception {
        int varListSize = variantList.size();
        int varPos = -1;
        boolean needNewRow = true;
        int varIndex = 0;
        int varIndex1 = 0;
        int varPos1 = 0;
        int indexREF = 1;
        int indexALT = 2;
        boolean fullMatch;
        int unmatchedNum = 0;
        char ref, alt;
        Variant var = variantList.get(varIndex);
        varPos = var.refStartPosition;
        byte[] currentLine = null;
        int currFilePostion = 0;
        byte[] sByte;
        String currChr;

        boolean hasNoScore = false;
        byte[] missingLabel = ".".getBytes();
        final byte semicolonByte = (byte) ';';
        int len;
        int varFeatureNum = fss.getAvailableFeatureIndex();
        int maxColNum = Math.max(indexCHROM, indexPOS);
        for (int i = 0; i < predicIndex.length; i++) {
            if (maxColNum < predicIndex[i]) {
                maxColNum = predicIndex[i];
            }
        }
        for (int i = 0; i < scoreIndexes.length; i++) {
            if (maxColNum < scoreIndexes[i]) {
                maxColNum = scoreIndexes[i];
            }
        }

        float[] scores = new float[scoreIndexes.length];
        String[] predicResults = new String[predicIndex.length];
        int[] cellsBT = new int[20];
        int[] cellDelimts = new int[maxColNum + 2];
        int shouldbeLen = maxColNum + 1;
        while (varIndex < varListSize) {
            // System.out.println(currentLine);    

            if (needNewRow) {
                // StringTokenizer st = new
                // StringTokenizer(currentLine.trim(), "\t");
                //Important: This function assumes each compressed block has at least one complete line.
                currentLine = br.readLine(cellDelimts);
                if (currentLine == null) {
                    break;
                }
                //this line may be trancated
                if (cellDelimts[0] < shouldbeLen) {
                    continue;
                }
                // System.out.println(new String(currentLine));
                //indexCHROM==0;
                if (indexCHROM >= 0) {
                    sByte = Arrays.copyOfRange(currentLine, 0, cellDelimts[indexCHROM + 1]);
                }

                //indexPOS=1;
                currFilePostion = parseInt(currentLine, indexPOS == 0 ? 0 : cellDelimts[indexPOS] + 1, cellDelimts[indexPOS + 1]);

                needNewRow = false;
            }

            if (varPos > currFilePostion) {
                needNewRow = true;
                continue;
            } else if (varPos < currFilePostion) {
                varIndex++;
                if (varIndex >= varListSize) {
                    break;
                }
                var = variantList.get(varIndex);
                varPos = var.refStartPosition;
                continue;
            }

            //ingore indel
            if (var.isIndel) {
                varIndex++;
                if (varIndex >= varListSize) {
                    break;
                }
                var = variantList.get(varIndex);
                varPos = var.refStartPosition;
                continue;
            }

            // in the variant genome, there is only one unique variant
            varIndex1 = varIndex;

            //this alogrithm does not work for indels
            //it must be equal and there may be multiple positions equal
            //indexREF=3
            sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexREF] + 1, cellDelimts[indexREF + 1]);
            ref = (char) sByte[0];
            //alt = currentLine.substring(cellDelimts[3] + 1, cellDelimts[4]);
            //indexALT=4
            //for reference data, sometimes we do not have alternative alleles
            sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexALT] + 1, cellDelimts[indexALT + 1]);
            alt = (char) sByte[0];
            unmatchedNum = 0;

            //The list may have variants with the same coordinates 
            do {
                var = variantList.get(varIndex1);
                varPos1 = var.refStartPosition;

                if (varPos1 > currFilePostion) {
                    break;
                }
                fullMatch = false;
                String[] altAlleles = var.getAltAlleles();
                for (String altA : altAlleles) {
                    if (altA == null) {
                        continue;
                    }
                    // assum there is ony one alternative allele Note
                    // the second condition is not safe
                    if (var.getRefAllele().charAt(0) == ref && altA.charAt(0) == alt
                            || var.getRefAllele().charAt(0) == alt && altA.charAt(0) == ref) {
                        fullMatch = true;
                        break;
                    }
                }
                if (fullMatch && var.scores1 == null) {
                    Arrays.fill(scores, Float.NaN);
                    Arrays.fill(predicResults, null);

                    for (int iCol = 0; iCol < scoreIndexes.length; iCol++) {
                        if (equal(currentLine, cellDelimts[scoreIndexes[iCol]] + 1, cellDelimts[scoreIndexes[iCol] + 1], missingLabel)) {
                            scores[iCol] = Float.NaN;
                        } else if (indexof(currentLine, cellDelimts[scoreIndexes[iCol]] + 1, cellDelimts[scoreIndexes[iCol] + 1], semicolonByte) < 0) {
                            scores[iCol] = parseFloat(currentLine, cellDelimts[scoreIndexes[iCol]] + 1, cellDelimts[scoreIndexes[iCol] + 1]);
                        } else {
                            // just use the first one
                            // System.out.println(tmpBuffer.toString());
                            Arrays.fill(cellsBT, -1);
                            len = tokenizeDelimiter(currentLine, cellDelimts[scoreIndexes[iCol]] + 1, cellDelimts[scoreIndexes[iCol] + 1], (byte) ';', cellsBT);

                            hasNoScore = true;
                            for (int s = 0; s < len; s++) {
                                if (!equal(currentLine, s == 0 ? cellsBT[s] : cellsBT[s] + 1, cellsBT[s + 1], missingLabel)) {
                                    scores[iCol] = parseFloat(currentLine, s == 0 ? cellsBT[s] : cellsBT[s] + 1, cellsBT[s + 1]);
                                    hasNoScore = false;
                                    break;
                                }
                            }

                            if (hasNoScore) {
                                scores[iCol] = Float.NaN;
                            }
                        }
                    }
                    for (int iCol = 0; iCol < predicIndex.length; iCol++) {
                        if (equal(currentLine, cellDelimts[predicIndex[iCol]] + 1, cellDelimts[predicIndex[iCol] + 1], missingLabel)) {
                            predicResults[iCol] = ".";
                        } else {
                            predicResults[iCol] = new String(currentLine, cellDelimts[predicIndex[iCol]] + 1, cellDelimts[predicIndex[iCol] + 1] - cellDelimts[predicIndex[iCol]] - 1);
                        }
                    }
                    var.scores1 = new float[scores.length];
                    System.arraycopy(scores, 0, var.scores1, 0, scores.length);
                    varAssignedScore++;
                    for (int k = 0; k < predicIndex.length; k++) {
                        var.setFeatureValue(varFeatureNum + k, predicResults[k]);
                    }

                } else {
                    unmatchedNum++;
                }
                varIndex1++;
            } while (varIndex1 < varListSize);

//The same coordinate but unmatched
            if (unmatchedNum == 0) {
                varIndex++;
                if (varIndex >= varListSize) {
                    break;
                }
                var = variantList.get(varIndex);
                varPos = var.refStartPosition;
            }
            //otherwise only move file's index
            needNewRow = true;
        }
        synchronized (fss) {
            fss.increaseCount(0, varAssignedScore);
            // System.out.println( ": " + varAssignedScore);
        }
    }

    public void dbNSFPAnnotSimple(IntArrayList positions, List<char[]> alleleList, List<float[]> scoreList, int[] regions) throws Exception {
        int varIndex1 = 0;
        int indexREF = 1;
        int indexALT = 2;
        char ref, alt;

        byte[] currentLine = null;
        int currFilePostion = 0;
        byte[] sByte;

        boolean hasNoScore = false;
        boolean allNoScore = false;
        byte[] missingLabel = ".".getBytes();
        final byte semicolonByte = (byte) ';';
        int len;

        int maxColNum = Math.max(indexCHROM, indexPOS);
        for (int i = 0; i < scoreIndexes.length; i++) {
            if (maxColNum < scoreIndexes[i]) {
                maxColNum = scoreIndexes[i];
            }
        }

        int[] cellsBT = new int[20];
        int[] cellDelimts = new int[maxColNum + 2];
        int shouldbeLen = maxColNum + 1;

        while ((currentLine = br.readLine(cellDelimts)) != null) {
            //this line may be trancated
            if (cellDelimts[0] < shouldbeLen) {
                continue;
            }
            // System.out.println(new String(currentLine));
            //indexCHROM==0;
            if (indexCHROM >= 0) {
                sByte = Arrays.copyOfRange(currentLine, 0, cellDelimts[indexCHROM + 1]);
            }

            //indexPOS=1;
            currFilePostion = parseInt(currentLine, indexPOS == 0 ? 0 : cellDelimts[indexPOS] + 1, cellDelimts[indexPOS + 1]);
            varIndex1 = Arrays.binarySearch(regions, currFilePostion);
            if (varIndex1 < 0) {
                varIndex1 = -varIndex1 - 1;
                //assume the regions are sorted in desendent order and not overlapped
                if (varIndex1 >= regions.length) {
//beyond all regions
                    break;
                } else if (varIndex1 % 2 == 0) {
                    //beyond a region
                    continue;
                }
            }

            //this alogrithm does not work for indels
            //it must be equal and there may be multiple positions equal
            //indexREF=3
            sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexREF] + 1, cellDelimts[indexREF + 1]);
            ref = (char) sByte[0];
            //alt = currentLine.substring(cellDelimts[3] + 1, cellDelimts[4]);
            //indexALT=4
            //for reference data, sometimes we do not have alternative alleles
            sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexALT] + 1, cellDelimts[indexALT + 1]);
            alt = (char) sByte[0];

            allNoScore = true;
            float[] scores = new float[scoreIndexes.length];
            Arrays.fill(scores, Float.NaN);

            for (int iCol = 0; iCol < scoreIndexes.length; iCol++) {
                if (equal(currentLine, cellDelimts[scoreIndexes[iCol]] + 1, cellDelimts[scoreIndexes[iCol] + 1], missingLabel)) {
                    scores[iCol] = Float.NaN;
                } else if (indexof(currentLine, cellDelimts[scoreIndexes[iCol]] + 1, cellDelimts[scoreIndexes[iCol] + 1], semicolonByte) < 0) {
                    scores[iCol] = parseFloat(currentLine, cellDelimts[scoreIndexes[iCol]] + 1, cellDelimts[scoreIndexes[iCol] + 1]);
                    allNoScore = false;
                } else {
                    // just use the first one
                    // System.out.println(tmpBuffer.toString());
                    Arrays.fill(cellsBT, -1);
                    len = tokenizeDelimiter(currentLine, cellDelimts[scoreIndexes[iCol]] + 1, cellDelimts[scoreIndexes[iCol] + 1], (byte) ';', cellsBT);

                    hasNoScore = true;
                    for (int s = 0; s < len; s++) {
                        if (!equal(currentLine, s == 0 ? cellsBT[s] : cellsBT[s] + 1, cellsBT[s + 1], missingLabel)) {
                            scores[iCol] = parseFloat(currentLine, s == 0 ? cellsBT[s] : cellsBT[s] + 1, cellsBT[s + 1]);
                            hasNoScore = false;
                            break;
                        }
                    }
                    allNoScore = false;
                    if (hasNoScore) {
                        scores[iCol] = Float.NaN;
                    }
                }
            }
            if (allNoScore) {
                continue;
            }
            positions.add(currFilePostion);
            alleleList.add(new char[]{ref, alt});
            scoreList.add(scores);
        }

    }

    //It only works when there is no indentifcal coordinates in the resource data
    public void dbNCFPAnnot() throws Exception {
        int varListSize = variantList.size();
        int varPos = -1;
        boolean needNewRow = true;
        int varIndex = 0;
        int varIndex1 = 0;
        int varPos1 = 0;
        int indexREF = 3;
        int indexALT = 4;
        boolean fullMatch;

        char ref, alt;
        byte[] sByte;
        byte[] alts;
        Variant var = variantList.get(varIndex);
        varPos = var.refStartPosition;
        byte[] currentLine = null;
        int currFilePostion = 0;

        String currChr;

        byte[] missingLabel = ".".getBytes();
        final byte semicolonByte = (byte) ';';

        int maxColNum = Math.max(indexCHROM, indexPOS);

        maxColNum = 5;
        int scoreNum = 10;
        float[] scores = new float[scoreNum];
        //the first column is Allele label
        scoreNum++;
        int[] cellsBT = new int[20];
        int[] cellDelimts = new int[maxColNum + 2];
        int shouldbeLen = maxColNum + 1;
        int indexA = 0;
        while (varIndex < varListSize) {
            // System.out.println(currentLine);            
            if (needNewRow) {
                // StringTokenizer st = new
                // StringTokenizer(currentLine.trim(), "\t");
                //Important: This function assumes each compressed block has at least one complete line.
                currentLine = br.readLine(cellDelimts);
                if (currentLine == null) {
                    break;
                }
                //this line may be trancated
                if (cellDelimts[0] < shouldbeLen) {
                    continue;
                }
                // System.out.println(new String(currentLine));
                //indexCHROM==0;
                if (indexCHROM >= 0) {
                    sByte = Arrays.copyOfRange(currentLine, 0, cellDelimts[indexCHROM + 1]);
                }

                //indexPOS=1;
                currFilePostion = parseInt(currentLine, indexPOS == 0 ? 0 : cellDelimts[indexPOS] + 1, cellDelimts[indexPOS + 1]);
                needNewRow = false;
            }

            if (varPos > currFilePostion) {
                needNewRow = true;
                continue;
            } else if (varPos < currFilePostion) {
                varIndex++;
                if (varIndex >= varListSize) {
                    break;
                }
                var = variantList.get(varIndex);
                varPos = var.refStartPosition;
                continue;
            }

            //ingore indel
            if (var.isIndel) {
                varIndex++;
                if (varIndex >= varListSize) {
                    break;
                }
                var = variantList.get(varIndex);
                varPos = var.refStartPosition;
                continue;
            }

            // in the variant genome, there is only one unique variant
            varIndex1 = varIndex;
            //  System.out.println(varIndex1);

            //this alogrithm does not work for indels
            //it must be equal and there may be multiple positions equal
            //indexREF=3
            sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexREF] + 1, cellDelimts[indexREF + 1]);
            ref = (char) sByte[0];
            //alt = currentLine.substring(cellDelimts[3] + 1, cellDelimts[4]);
            //indexALT=4
            //for reference data, sometimes we do not have alternative alleles
            sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexALT] + 1, cellDelimts[indexALT + 1]);
            alts = sByte;

            //The list may have variants with the same coordinates 
            do {
                var = variantList.get(varIndex1);
                varPos1 = var.refStartPosition;

                if (varPos1 > currFilePostion) {
                    break;
                }
                fullMatch = false;
                String[] altAlleles = var.getAltAlleles();
                indexA = 0;
                Arrays.fill(scores, Float.NaN);
                for (String altA : altAlleles) {
                    if (altA == null) {
                        continue;
                    }
                    // assum there is ony one alternative allele Note
                    // the second condition is not safe
                    for (int i = 0; i < alts.length; i += 2) {
                        if (var.getRefAllele().charAt(0) == ref && altA.charAt(0) == alts[i]
                                || var.getRefAllele().charAt(0) == alts[i] && altA.charAt(0) == ref) {
                            fullMatch = true;
                            indexA = i / 2;
                            break;
                        }
                    }

                    if (fullMatch) {
                        //System.out.println(new String(currentLine));
                        Arrays.fill(cellsBT, -1);
                        tokenizeDelimiter(currentLine, cellDelimts[5] + 1, cellDelimts[6], (byte) ',', cellsBT);
                        tokenizeDelimiter(currentLine, (indexA == 0 ? cellsBT[0] : cellsBT[indexA] + 1), cellsBT[indexA + 1], (byte) '|', cellsBT);
                        for (int iCol = 1; iCol < scoreNum; iCol++) {
                            if (equal(currentLine, cellsBT[iCol] + 1, cellsBT[iCol + 1], missingLabel)) {
                                scores[iCol - 1] = Float.NaN;
                            } else if (cellsBT[iCol + 1] > cellsBT[iCol] + 1) {
                                scores[iCol - 1] = parseFloat(currentLine, cellsBT[iCol] + 1, cellsBT[iCol + 1], iCol == 8);
                            } else // just use the first one
                            //System.out.println(tmpBuffer.toString()); 
                            {
                                //empty
                                scores[iCol - 1] = Float.NaN;
                            }
                        }
                        var.scores2 = new float[scores.length];
                        System.arraycopy(scores, 0, var.scores2, 0, scores.length);

                        varAssignedScore++;
                    }
                }
                varIndex1++;
            } while (varIndex1 < varListSize);
//The same coordinate but unmatched

            //otherwise only move file's index
            needNewRow = true;
        }
        for (varIndex = 0; varIndex < varListSize; varIndex++) {
            var = variantList.get(varIndex);
            if (var.scores2 == null) {
                var.scores2 = new float[scores.length];
                Arrays.fill(var.scores2, Float.NaN);
            }
        }
        //System.out.println(" : " + varAssignedScore);
        synchronized (br) {
            // fss.increaseCount(0, varAssignedScore);
            // System.out.println( ": " + varAssignedScore);
        }
    }

    public void markByANNOVARefFormat() {
        indexCHROM = 0;
        indexPOS = 1;
        int indexREF = 2;
        int indexALT = 3;
        int indexMAF = 4;

        if (variantList.isEmpty()) {
            return;
        }
        int feautreNum = ass.getAvailableFeatureIndex();

        byte[] currentLine = null;
        try {
            int maxColNum = indexCHROM;
            maxColNum = Math.max(maxColNum, indexPOS);
            maxColNum = Math.max(maxColNum, indexREF);
            maxColNum = Math.max(maxColNum, indexALT);
            maxColNum = Math.max(maxColNum, indexMAF);

            int lineCounter = 0;

            int filePosition = -1;

            String ref;
            String alt;

            float maf;
            String[] alts;
            int[] mafIndex;
            int[] varIndex = null;
            boolean isDeleion = false;
            boolean isInsertion = false;

            char[] backSpaces = null;
            int delNum = 0;
            int existVarNum = 0;
            StringBuilder sb = new StringBuilder();
            boolean hitOnce = false;
            int[] cellDelimts = new int[maxColNum + 2];
            int shouldbeLen = maxColNum + 1;
            byte[] sByte;
            byte[] currChr = chromosome.getName().getBytes();
            int indexTmp = 0;
            int tmpPos = 0;
            while ((currentLine = br.readLine(cellDelimts)) != null) {
                lineCounter++;

                //this line may be trancated
                if (cellDelimts[0] < shouldbeLen) {
                    continue;
                }

                //StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles
                //sByte = Arrays.copyOfRange(currentLine, indexCHROM == 0 ? 0 : cellDelimts[indexCHROM] + 1, cellDelimts[indexCHROM + 1]);
                //  System.err.println(new String(currentLine));
                if (!equal(currentLine, indexCHROM == 0 ? 0 : cellDelimts[indexCHROM] + 1, cellDelimts[indexCHROM + 1], currChr)) {
                    continue;
                }
                //indexPOS=1;
                filePosition = parseInt(currentLine, indexPOS == 0 ? 0 : cellDelimts[indexPOS] + 1, cellDelimts[indexPOS + 1]);

                //indexREF=3
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexREF] + 1, cellDelimts[indexREF + 1]);
                ref = new String(sByte);
                //alt = currentLine.substring(cellDelimts[3] + 1, cellDelimts[4]);
                //indexALT=4
                //for reference data, sometimes we do not have alternative alleles
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexALT] + 1, cellDelimts[indexALT + 1]);
                alt = new String(sByte);

                alts = alt.split(",");
                mafIndex = tokenizeDelimiter(currentLine, cellDelimts[indexMAF] + 1, cellDelimts[indexMAF + 1], (byte) ',');
                if (mafIndex[0] == mafIndex[mafIndex.length - 1]) {
                    continue;
                }
                hitOnce = false;
                tmpPos = 0;
                //once the variant is in db, it at least has a zero freq
                for (int s = 0; s < alts.length; s++) {
                    if (alts[s] == null || alts[s].isEmpty()) {
                        continue;
                    }
                    maf = Float.NaN;
                    if (mafIndex != null && s < mafIndex.length - 1) {
                        //this missing score is denoted by .
                        indexTmp = (s == 0 ? mafIndex[s] : mafIndex[s] + 1);

                        if (indexTmp >= currentLine.length) {
                            continue;
                        }
                        if (currentLine[indexTmp] != '.') {
                            maf = parseFloat(currentLine, s == 0 ? mafIndex[s] : mafIndex[s] + 1, mafIndex[s + 1]);
                        }
                    }

                    tmpPos = filePosition;
                    isDeleion = false;
                    isInsertion = false;

                    alt = alts[s];
                    //1 13211293 13211294 TC - comments: rs59770105, a 2-bp deletion
                    //deletion
                    //format:1	45113	-	0TATGG	0.715732
///1	53599	CTA	3	0.890916
//1	223450	CT	1	0.207385
                    //1	229450	C	T,G	0.207385,0.1
                    //1	229450	C	T,G	0.207385,

                    if (alt.charAt(0) == '0') {
                        isInsertion = true;
                    } else if (alt.charAt(0) == '-' || (alt.charAt(0) - '0' <= 9 && alt.charAt(0) - '0' > 0)) {
                        isDeleion = true;
                        tmpPos--;
                    }

                    varIndex = chromosome.lookupVariantIndexes(tmpPos);
                    if (varIndex == null) {
                        continue;
                    }

                    // System.out.println(fileChr);
                    for (int index : varIndex) {
                        Variant var = variantList.get(index);
                        //althoug it is not thread-safe; it should be OK because it is will given the same value by different threads
                        if (var.hasBeenAcced) {
                            continue;
                        }

                        if (isDeleion || isInsertion) {
                            if (var.isIndel) {
                                String varRef = var.getRefAllele();
                                String[] altAlleles = var.getAltAlleles();
                                //keep variants with score less than minAlleleFreq
                                for (String varAlt : altAlleles) {
                                    //insertion in 1KG
                                    if (isInsertion) {
                                        if (varAlt.substring(1).equals(alt.substring(1))) {
                                            //record the maximal allele frequencies
                                            if (var.altAF == -1 || Float.isNaN(var.altAF) || (maf > var.altAF)) {
                                                var.altAF = maf;
                                            }
                                            hitOnce = true;
                                            if (Float.isNaN(maf)) {
                                                var.setFeatureValue(feautreNum, ".");
                                            } else {
                                                var.setFeatureValue(feautreNum, String.valueOf(maf));
                                            }
                                            break;
                                        }
                                    } else if (alt.charAt(0) != '0') {
                                        //deletion in 1KG
                                        sb.delete(0, sb.length());
                                        for (int t = 0; t < varAlt.length(); t++) {
                                            if (varAlt.charAt(t) == '-') {
                                                sb.append(varRef.charAt(t));
                                            }
                                        }
                                        if (alt.charAt(0) == '-') {
                                            delNum = ref.length();
                                        } else {
                                            delNum = Util.parseInt(alt);
                                        }
                                        if (sb.toString().equals(ref.substring(ref.length() - delNum))) {
                                            //record the maximal allele frequencies
                                            if (var.altAF == -1 || Float.isNaN(var.altAF) || (maf > var.altAF)) {
                                                // if (Float.isNaN(var.altAF) || (!Float.isNaN(score) && score > var.altAF)) {
                                                var.altAF = maf;
                                            }
                                            hitOnce = true;
                                            if (Float.isNaN(maf)) {
                                                var.setFeatureValue(feautreNum, ".");
                                            } else {
                                                var.setFeatureValue(feautreNum, String.valueOf(maf));
                                            }
                                            var.hasBeenAcced = true;
                                            break;
                                        }
                                    }
                                }
                            } else {
                                continue;
                            }
                        } else if (var.isIndel) {
                            continue;
                        } else {
                            String[] altAlleles = var.getAltAlleles();
                            for (String str : altAlleles) {
                                if (str.charAt(0) == alt.charAt(0)) {
                                    //record treadVariantsInFileOnlyFastTokenhe maximal allele frequencies
                                    if (var.altAF == -1 || Float.isNaN(var.altAF) || (maf > var.altAF)) {
                                        //if (Float.isNaN(var.altAF) || (!Float.isNaN(score) && score > var.altAF)) {
                                        var.altAF = maf;
                                    }
                                    hitOnce = true;
                                    if (Float.isNaN(maf)) {
                                        var.setFeatureValue(feautreNum, ".");
                                    } else {
                                        var.setFeatureValue(feautreNum, String.valueOf(maf));
                                    }
                                    var.hasBeenAcced = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (hitOnce) {
                    existVarNum++;
                }
            }
            synchronized (ass) {
                ass.setLeftNum(existVarNum + ass.getLeftNum());
                ass.setAnnotNum(existVarNum + ass.getAnnotNum());
                ass.setTotalNum(lineCounter + ass.getTotalNum());
                // System.out.println(chromosome.getName()+": "+lineCounter);
            }

        } catch (Exception ex) {
            if (currentLine != null) {
                String ssss = new String(currentLine);
                System.err.println("Errors in a row: " + ssss);
            }
            ex.printStackTrace();
        }
    }

    public void markByVCFAllelesFormat() {
        indexCHROM = 0;
        indexPOS = 1;
        int indexREF = 2;
        int indexALT = 3;

        if (variantList.isEmpty()) {
            return;
        }
        int feautreNum = ass.getAvailableFeatureIndex();
        int featureIndex;
        byte[] currentLine = null;
        try {
            int maxColNum = indexCHROM;
            maxColNum = Math.max(maxColNum, indexPOS);
            maxColNum = Math.max(maxColNum, indexREF);
            maxColNum = Math.max(maxColNum, indexALT);
            for (int i = 0; i < altFreqColIndexes.length; i++) {
                maxColNum = Math.max(maxColNum, altFreqColIndexes[i]);
            }

            int lineCounter = 0;

            int filePosition = -1;

            String ref;
            String alt;

            double maf;
            String[] alts;
            int[] mafIndexInCol;
            int[] varIndex = null;
            boolean isDeleion = false;
            boolean isInsertion = false;

            char[] backSpaces = null;
            int delNum = 0;
            int existVarNum = 0;
            StringBuilder sb = new StringBuilder();
            boolean hitOnce = false;
            int[] cellDelimts = new int[maxColNum + 2];
            int shouldbeLen = maxColNum + 1;
            byte[] sByte;
            byte[] currChr = chromosome.getName().getBytes();
            int indexTmp = 0;
            int tmpPos = 0;
            int indexMAF;
            while ((currentLine = br.readLine(cellDelimts)) != null) {
                lineCounter++;

                //this line may be trancated
                if (cellDelimts[0] < shouldbeLen) {
                    continue;
                }

                //StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles
                //sByte = Arrays.copyOfRange(currentLine, indexCHROM == 0 ? 0 : cellDelimts[indexCHROM] + 1, cellDelimts[indexCHROM + 1]);
                //  System.err.println(new String(currentLine));
                if (!equal(currentLine, indexCHROM == 0 ? 0 : cellDelimts[indexCHROM] + 1, cellDelimts[indexCHROM + 1], currChr)) {
                    continue;
                }
                //indexPOS=1;
                filePosition = parseInt(currentLine, indexPOS == 0 ? 0 : cellDelimts[indexPOS] + 1, cellDelimts[indexPOS + 1]);

                //indexREF=3
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexREF] + 1, cellDelimts[indexREF + 1]);
                ref = new String(sByte);
                //alt = currentLine.substring(cellDelimts[3] + 1, cellDelimts[4]);
                //indexALT=4
                //for reference data, sometimes we do not have alternative alleles
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexALT] + 1, cellDelimts[indexALT + 1]);
                alt = new String(sByte);

                alts = alt.split(",");

                hitOnce = false;
                tmpPos = 0;

                //once the variant is in db, it at least has a zero freq
                for (int s = 0; s < alts.length; s++) {
                    if (alts[s] == null || alts[s].isEmpty()) {
                        continue;
                    }
                    tmpPos = filePosition;
                    isDeleion = false;
                    isInsertion = false;

                    alt = alts[s];
                    //1 13211293 13211294 TC - comments: rs59770105, a 2-bp deletion
                    //deletion
                    //format:1	3739822	AAG	A  0.715732
                    //insertion 
///1	8064650	AAAAC	A,AAAACAAAC	0.890916
                    //1	229450	C	T,G	0.207385,0.1
                    //1	229450	C	T,G	0.207385,

                    if (ref.length() > alt.length()) {
                        //deletiion 
                        isDeleion = true;
                        //Assume on
                        tmpPos += (alt.length() - 1);
                    } else if (ref.length() < alt.length()) {
                        isInsertion = true;
                        tmpPos += (ref.length() - 1);
                    }

                    varIndex = chromosome.lookupVariantIndexes(tmpPos);
                    if (varIndex == null) {
                        continue;
                    }
                    //System.err.println(new String(currentLine));
                   
                    featureIndex = feautreNum - 1;
                    for (int i = 0; i < altFreqColIndexes.length; i++) {
                        featureIndex++;
                        indexMAF = altFreqColIndexes[i];
                        mafIndexInCol = tokenizeDelimiter(currentLine, cellDelimts[indexMAF] + 1, cellDelimts[indexMAF + 1], (byte) ',');
                        if (mafIndexInCol[0] == mafIndexInCol[mafIndexInCol.length - 1]) {
                            continue;
                        }

                        maf = Double.NaN;
                        if (s < mafIndexInCol.length - 1) {
                            //this missing score is denoted by .
                            indexTmp = (s == 0 ? mafIndexInCol[s] : mafIndexInCol[s] + 1);

                            if (indexTmp >= currentLine.length) {
                                continue;
                            }
                            if (currentLine[indexTmp] != '.') {
                                maf = parseDouble(currentLine, s == 0 ? mafIndexInCol[s] : mafIndexInCol[s] + 1, mafIndexInCol[s + 1]);
                            }
                        }

                        // System.out.println(fileChr);
                        for (int index : varIndex) {
                            Variant var = variantList.get(index);
                            //althoug it is not thread-safe; it should be OK because it is will given the same value by different threads
                            if (var.hasBeenAcced) {
                                //continue;
                            }

                            if (isDeleion || isInsertion) {
                                if (var.isIndel) {
                                    String varRef = var.getRefAllele();
                                    String[] altAlleles = var.getAltAlleles();
                                    //keep variants with score less than minAlleleFreq
                                    for (String varAlt : altAlleles) {
                                        //insertion in 1KG
                                        if (isInsertion) {
                                            if (varAlt.substring(1).equals(alt.substring(1))) {
                                                //record the maximal allele frequencies
                                                if (var.altAF == -1 || Float.isNaN(var.altAF) || (maf > var.altAF)) {
                                                    var.altAF = (float) maf;
                                                }
                                                hitOnce = true;
                                                if (Double.isNaN(maf)) {
                                                    var.setFeatureValue(featureIndex, ".");
                                                } else {
                                                    var.setFeatureValue(featureIndex, String.valueOf(maf));
                                                }
                                                break;
                                            }
                                        } else if (alt.charAt(0) != '0') {
                                            //deletion in 1KG
                                            sb.delete(0, sb.length());
                                            for (int t = 0; t < varAlt.length(); t++) {
                                                if (varAlt.charAt(t) == '-') {
                                                    sb.append(varRef.charAt(t));
                                                }
                                            }
                                            if (alt.charAt(0) == '-') {
                                                delNum = ref.length();
                                            } else {
                                                delNum = Util.parseInt(alt);
                                            }
                                            if (sb.toString().equals(ref.substring(ref.length() - delNum))) {
                                                //record the maximal allele frequencies
                                                if (var.altAF == -1 || Float.isNaN(var.altAF) || (maf > var.altAF)) {
                                                    // if (Float.isNaN(var.altAF) || (!Float.isNaN(score) && score > var.altAF)) {
                                                    var.altAF = (float) maf;
                                                }
                                                hitOnce = true;
                                                if (Double.isNaN(maf)) {
                                                    var.setFeatureValue(featureIndex, ".");
                                                } else {
                                                    var.setFeatureValue(featureIndex, String.valueOf(maf));
                                                }
                                                var.hasBeenAcced = true;
                                                break;
                                            }
                                        }
                                    }
                                } else {
                                    continue;
                                }
                            } else if (var.isIndel) {
                                continue;
                            } else {
                                String[] altAlleles = var.getAltAlleles();
                                for (String str : altAlleles) {
                                    if (str.charAt(0) == alt.charAt(0)) {
                                        //record treadVariantsInFileOnlyFastTokenhe maximal allele frequencies
                                        if (var.altAF == -1 || Float.isNaN(var.altAF) || (maf > var.altAF)) {
                                            //if (Float.isNaN(var.altAF) || (!Float.isNaN(score) && score > var.altAF)) {
                                            var.altAF = (float) maf;
                                        }
                                        hitOnce = true;
                                        if (Double.isNaN(maf)) {
                                            var.setFeatureValue(featureIndex, ".");
                                        } else {
                                            var.setFeatureValue(featureIndex, String.valueOf(maf));
                                        }
                                        var.hasBeenAcced = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }

                }
                if (hitOnce) {
                    existVarNum++;
                }

            }
            synchronized (ass) {
                ass.setLeftNum(existVarNum + ass.getLeftNum());
                ass.setAnnotNum(existVarNum + ass.getAnnotNum());
                ass.setTotalNum(lineCounter + ass.getTotalNum());
                // System.out.println(chromosome.getName()+": "+lineCounter);
            }

        } catch (Exception ex) {
            if (currentLine != null) {
                String ssss = new String(currentLine);
                System.err.println("Errors in a row: " + ssss);
            }
            ex.printStackTrace();
        }
    }

    public void markByANNOVARefFormatSimpleVar(IntArrayList positions, List<char[]> alleleList, List<double[]> frequnceList) {
        indexCHROM = 0;
        indexPOS = 1;
        int indexREF = 2;
        int indexALT = 3;
        int indexMAF = 4;

        byte[] currentLine = null;
        try {
            int maxColNum = indexCHROM;
            maxColNum = Math.max(maxColNum, indexPOS);
            maxColNum = Math.max(maxColNum, indexREF);
            maxColNum = Math.max(maxColNum, indexALT);
            maxColNum = Math.max(maxColNum, indexMAF);

            int lineCounter = 0;

            int filePosition = -1;

            String ref;
            String alt;

            double maf;
            String[] alts;
            int[] mafIndex;

            int varIndex, varIndex0, varIndex1;

            boolean isDeleion = false;
            boolean isInsertion = false;

            int existVarNum = 0;

            boolean hitOnce = false;
            int[] cellDelimts = new int[maxColNum + 2];
            int shouldbeLen = maxColNum + 1;
            byte[] sByte;
            byte[] currChr = chromosome.getName().getBytes();
            int indexTmp = 0;
            int tmpPos = 0;
            while ((currentLine = br.readLine(cellDelimts)) != null) {
                lineCounter++;

                //this line may be trancated
                if (cellDelimts[0] < shouldbeLen) {
                    continue;
                }

                //StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles
                //sByte = Arrays.copyOfRange(currentLine, indexCHROM == 0 ? 0 : cellDelimts[indexCHROM] + 1, cellDelimts[indexCHROM + 1]);
                //  System.err.println(new String(currentLine));
                if (!equal(currentLine, indexCHROM == 0 ? 0 : cellDelimts[indexCHROM] + 1, cellDelimts[indexCHROM + 1], currChr)) {
                    continue;
                }
                //indexPOS=1;
                filePosition = parseInt(currentLine, indexPOS == 0 ? 0 : cellDelimts[indexPOS] + 1, cellDelimts[indexPOS + 1]);

                //indexREF=3
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexREF] + 1, cellDelimts[indexREF + 1]);
                ref = new String(sByte);
                //alt = currentLine.substring(cellDelimts[3] + 1, cellDelimts[4]);
                //indexALT=4
                //for reference data, sometimes we do not have alternative alleles
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexALT] + 1, cellDelimts[indexALT + 1]);
                alt = new String(sByte);

                alts = alt.split(",");
                mafIndex = tokenizeDelimiter(currentLine, cellDelimts[indexMAF] + 1, cellDelimts[indexMAF + 1], (byte) ',');
                if (mafIndex[0] == mafIndex[mafIndex.length - 1]) {
                    continue;
                }
                hitOnce = false;
                tmpPos = 0;
                //once the variant is in db, it at least has a zero freq
                for (int s = 0; s < alts.length; s++) {
                    if (alts[s] == null || alts[s].isEmpty()) {
                        continue;
                    }
                    maf = Double.NaN;
                    if (s < mafIndex.length - 1) {
                        //this missing score is denoted by .
                        indexTmp = (s == 0 ? mafIndex[s] : mafIndex[s] + 1);

                        if (indexTmp >= currentLine.length) {
                            continue;
                        }
                        if (currentLine[indexTmp] != '.') {
                            maf = parseDouble(currentLine, s == 0 ? mafIndex[s] : mafIndex[s] + 1, mafIndex[s + 1]);
                        }
                    }

                    tmpPos = filePosition;
                    isDeleion = false;
                    isInsertion = false;

                    alt = alts[s];
                    //1 13211293 13211294 TC - comments: rs59770105, a 2-bp deletion
                    //deletion
                    //format:1	45113	-	0TATGG	0.715732
///1	53599	CTA	3	0.890916
//1	223450	CT	1	0.207385
                    //1	229450	C	T,G	0.207385,0.1
                    //1	229450	C	T,G	0.207385,

                    if (alt.charAt(0) == '0') {
                        isInsertion = true;
                    } else if (alt.charAt(0) == '-' || (alt.charAt(0) - '0' <= 9 && alt.charAt(0) - '0' > 0)) {
                        isDeleion = true;
                        tmpPos--;
                    }

                    varIndex = positions.binarySearch(tmpPos);
                    if (varIndex < 0) {
                        continue;
                    }
                    varIndex0 = varIndex;
                    while (positions.getQuick(varIndex0) == tmpPos) {
                        varIndex0--;
                    }
                    varIndex0++;
                    varIndex1 = varIndex;
                    while (positions.getQuick(varIndex1) == tmpPos) {
                        varIndex1++;
                    }

                    // System.out.println(fileChr);
                    for (int index = varIndex0; index < varIndex1; index++) {
                        char[] allels = alleleList.get(index);
                        if (allels[1] == alt.charAt(0) && ref.charAt(0) == allels[0]) {
                            frequnceList.get(index)[0] = maf;
                        }

                    }
                }
                if (hitOnce) {
                    existVarNum++;
                }
            }

        } catch (Exception ex) {
            if (currentLine != null) {
                String ssss = new String(currentLine);
                System.err.println("Errors in a row: " + ssss);
            }
            ex.printStackTrace();
        }
    }

    //I stop writing this new merge function; because it will not work when both variant list and file list have identifcal variants
    public void markByANNOVARefFormatMerge() {
        indexCHROM = 0;
        indexPOS = 1;
        int indexREF = 2;
        int indexALT = 3;
        int indexMAF = 4;

        if (variantList.isEmpty()) {
            return;
        }
        int feautreNum = ass.getAvailableFeatureIndex();

        byte[] currentLine = null;
        try {
            int maxColNum = indexCHROM;
            maxColNum = Math.max(maxColNum, indexPOS);
            maxColNum = Math.max(maxColNum, indexREF);
            maxColNum = Math.max(maxColNum, indexALT);
            maxColNum = Math.max(maxColNum, indexMAF);

            int lineCounter = 0;

            int filePosition = -1;

            String ref;
            String alt;

            float maf = Float.NaN;
            String[] alts;
            int[] mafIndex;

            boolean isDeleion = false;
            boolean isInsertion = false;

            char[] backSpaces = null;
            int delNum = 0;
            int existVarNum = 0;
            StringBuilder sb = new StringBuilder();
            boolean hitOnce = false;
            int[] cellDelimts = new int[maxColNum + 2];
            int shouldbeLen = maxColNum + 1;
            byte[] sByte;
            byte[] currChr = chromosome.getName().getBytes();
            int varListSize = variantList.size();
            boolean needNewRow = true;
            int varI = 0, varIndex1;
            int varPos1, varPos;
            Variant var = variantList.get(varI);
            varPos = var.refStartPosition;
            int tmpPos = 0;
            while (varI < varListSize) {
                // System.out.println(currentLine);       
                if (needNewRow) {
                    // StringTokenizer st = new
                    // StringTokenizer(currentLine.trim(), "\t");
                    //Important: This function assumes each compressed block has at least one complete line.
                    currentLine = br.readLine(cellDelimts);
                    if (currentLine == null) {
                        break;
                    }
                    //this line may be trancated
                    if (cellDelimts[0] < shouldbeLen) {
                        continue;
                    }
                    //StringTokenizer st = new StringTokenizer(currentLine.trim());
                    //initialize varaibles
                    if (indexCHROM >= 0) {
                        sByte = Arrays.copyOfRange(currentLine, indexCHROM == 0 ? 0 : cellDelimts[indexCHROM] + 1, cellDelimts[indexCHROM + 1]);
                        if (!equal(currentLine, indexCHROM == 0 ? 0 : cellDelimts[indexCHROM] + 1, cellDelimts[indexCHROM + 1], currChr)) {
                            continue;
                        }

                    }

                    //indexPOS=1;
                    filePosition = parseInt(currentLine, indexPOS == 0 ? 0 : cellDelimts[indexPOS] + 1, cellDelimts[indexPOS + 1]);
                    needNewRow = false;
                }

                //indexREF=3
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexREF] + 1, cellDelimts[indexREF + 1]);
                ref = new String(sByte);
                //alt = currentLine.substring(cellDelimts[3] + 1, cellDelimts[4]);
                //indexALT=4
                //for reference data, sometimes we do not have alternative alleles
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexALT] + 1, cellDelimts[indexALT + 1]);
                alt = new String(sByte);

                alts = alt.split(",");

                mafIndex = tokenizeDelimiter(currentLine, cellDelimts[indexMAF] + 1, cellDelimts[indexMAF + 1], (byte) ',');

                // in the variant genome, there is only one unique variant
                varIndex1 = varI;
                hitOnce = false;
                //The list may have variants with the same coordinates 
                do {
                    var = variantList.get(varIndex1);
                    varPos1 = var.refStartPosition;

                    if (varPos1 > filePosition) {
                        break;
                    }

                    for (int s = 0; s < alts.length; s++) {
                        if (alts[s] == null || alts[s].isEmpty()) {
                            continue;
                        }
                        maf = Float.NaN;
                        if (mafIndex != null && s < mafIndex.length - 1) {
                            //this missing score is denoted by .
                            if (currentLine[s == 0 ? mafIndex[s] : mafIndex[s] + 1] != '.') {
                                maf = parseFloat(currentLine, s == 0 ? mafIndex[s] : mafIndex[s] + 1, mafIndex[s + 1]);
                            }
                        }

                        tmpPos = filePosition;
                        isDeleion = false;
                        isInsertion = false;

                        alt = alts[s];
                        //deletion
                        //format:1	45113	-	0TATGG	0.715732
///1	53599	CTA	3	0.890916
//1	223450	CT	1	0.207385
                        //1	229450	C	T,G	0.207385,0.1
                        //1	229450	C	T,G	0.207385,

                        if (alt.charAt(0) == '0') {
                            isInsertion = true;
                        } else if (alt.charAt(0) - '0' <= 9 && alt.charAt(0) - '0' > 0) {
                            isDeleion = true;
                            tmpPos--;
                        }

                        if (isDeleion || isInsertion) {
                            if (var.isIndel) {
                                String varRef = var.getRefAllele();
                                String[] altAlleles = var.getAltAlleles();
                                //keep variants with score less than minAlleleFreq
                                for (String varAlt : altAlleles) {
                                    //insertion in 1KG
                                    if (isInsertion) {
                                        if (varAlt.substring(1).equals(alt.substring(1))) {
                                            //record the maximal allele frequencies
                                            if (var.altAF == -1 || Float.isNaN(var.altAF) || (maf > var.altAF)) {
                                                var.altAF = maf;
                                            }
                                            hitOnce = true;
                                            if (Float.isNaN(maf)) {
                                                var.setFeatureValue(feautreNum, ".");
                                            } else {
                                                var.setFeatureValue(feautreNum, String.valueOf(maf));
                                            }
                                            break;
                                        }
                                    } else if (alt.charAt(0) != '0') {
                                        //deletion in 1KG
                                        sb.delete(0, sb.length());
                                        for (int t = 0; t < varAlt.length(); t++) {
                                            if (varAlt.charAt(t) == '-') {
                                                sb.append(varRef.charAt(t));
                                            }
                                        }

                                        delNum = Util.parseInt(alt);
                                        if (sb.toString().equals(ref.substring(ref.length() - delNum))) {
                                            //record the maximal allele frequencies
                                            if (var.altAF == -1 || Float.isNaN(var.altAF) || (maf > var.altAF)) {
                                                // if (Float.isNaN(var.altAF) || (!Float.isNaN(score) && score > var.altAF)) {
                                                var.altAF = maf;
                                            }
                                            hitOnce = true;
                                            if (Float.isNaN(maf)) {
                                                var.setFeatureValue(feautreNum, ".");
                                            } else {
                                                var.setFeatureValue(feautreNum, String.valueOf(maf));
                                            }
                                            break;
                                        }
                                    }
                                }
                            } else {
                                continue;
                            }
                        } else if (var.isIndel) {
                            continue;
                        } else {
                            String[] altAlleles = var.getAltAlleles();
                            for (String str : altAlleles) {
                                if (str.charAt(0) == alt.charAt(0)) {
                                    //record treadVariantsInFileOnlyFastTokenhe maximal allele frequencies
                                    if (var.altAF == -1 || Float.isNaN(var.altAF) || (maf > var.altAF)) {
                                        //if (Float.isNaN(var.altAF) || (!Float.isNaN(score) && score > var.altAF)) {
                                        var.altAF = maf;
                                    }
                                    hitOnce = true;
                                    if (Float.isNaN(maf)) {
                                        var.setFeatureValue(feautreNum, ".");
                                    } else {
                                        var.setFeatureValue(feautreNum, String.valueOf(maf));
                                    }
                                    break;
                                }
                            }
                        }

                    }
                    varIndex1++;
                } while (varIndex1 < varListSize);
//The same coordinate but unmatched
                if (!hitOnce) {
                    varI++;
                    if (varI >= varListSize) {
                        break;
                    }
                    var = variantList.get(varI);
                    varPos = var.refStartPosition;
                } else {
                    existVarNum++;
                }
                //otherwise only move file's index
                needNewRow = true;

                synchronized (ass) {
                    ass.setLeftNum(existVarNum + ass.getLeftNum());
                    ass.setAnnotNum(existVarNum + ass.getAnnotNum());
                    ass.setTotalNum(lineCounter + ass.getTotalNum());
                }
            }
        } catch (Exception ex) {
            if (currentLine != null) {
                String ssss = new String(currentLine);
                System.err.println("Errors in a row: " + ssss);
            }
            ex.printStackTrace();
        }
    }

    @Override
    public String call() throws Exception {
        if (annotationType == 0) {
            dbNSFPAnnot();
        } else if (annotationType == 1) {
            markByANNOVARefFormat();
        } else if (annotationType == 2) {
            dbNCFPAnnot();
        } else if (annotationType == 3) {
            markByVCFAllelesFormat();
        }

        fireTaskComplete();
        return "";
    }

    private int parseInt(final byte[] s, int start, int end) {
        // Check for a sign.
        int num = 0;
        int sign = -1;
        int i = start;
        //ACSII
        //'0' : 48 
        //'9': 57
        //'-' 45
        //'.' 46
        //'e':101
        //'E':69
        //' ': 32
        while (s[i] == 32) {
            i++;
        }

        final byte ch = s[i++];
        if (ch == 45) {
            sign = 1;
        } else {
            num = 48 - ch;
        }

        // Build the number. 
        while (i < end) {
            if (s[i] == 46) {
                return sign * num;
            } else if (s[i] < 48 || s[i] > 57) {
                i++;
            } else {
                num = num * 10 + 48 - s[i++];
            }
        }
        return sign * num;
    }

    private boolean equal(byte[] src, int start, int end, byte[] tar) {
        if (end - start != tar.length) {
            return false;
        }
        for (int i = start; i < end; i++) {
            if (src[i] != tar[i - start]) {
                return false;
            }
        }
        return true;
    }

    private int indexof(byte[] src, int start, int end, byte tar) {
        if (end - start < 1) {
            return -1;
        }
        for (int i = start; i < end; i++) {
            if (tar == src[i]) {
                return i;
            }
        }
        return -1;
    }

    private float parseFloat(byte[] f, int start, int end) {
        double ret = 0f;         // return value
        int pos = start;          // read pointer position
        int part = 0;          // the current part (int, float and sci parts of the number)
        boolean neg = false;      // true if part is a negative number
        // the max long is 2147483647
        final int MAX_INT_BIT = 9;

        //ACSII
        //'0' : 48 
        //'9': 57
        //'-' 45
        //'.' 46
        //'e':101
        //'E':69
        while (f[pos] == ' ') {
            pos++;
        }
        // find start
        while (pos < end && (f[pos] < 48 || f[pos] > 57) && f[pos] != 45 && f[pos] != 46) {
            pos++;
        }

        // sign
        if (f[pos] == 45) {
            neg = true;
            pos++;
        }

        // integer part
        while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
            part = part * 10 + (f[pos++] - 48);
        }
        ret = neg ? (double) (part * -1) : (double) part;

        // float part
        if (pos < end && f[pos] == 46) {
            pos++;
            int mul = 1;
            part = 0;
            int num = 0;
            while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
                num++;
                if (num <= MAX_INT_BIT) {
                    part = part * 10 + (f[pos] - 48);
                    mul *= 10;
                }
                pos++;
            }
            ret = neg ? ret - (double) part / (double) mul : ret + (double) part / (double) mul;
        }

        // scientific part
        if (pos < end && (f[pos] == 101 || f[pos] == 69)) {
            pos++;
            if (pos < end) {
                neg = (f[pos] == 45);
                pos++;
                part = 0;
                while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
                    part = part * 10 + (f[pos++] - 48);
                }
                if (neg) {
                    ret = ret / Math.pow(10, part);
                } else {
                    ret = ret * Math.pow(10, part);
                }
            }
        }
        return (float) ret;
    }

    private float parseFloat(byte[] f, int start, int end, boolean takeLog) {
        double ret = 0f;         // return value
        int pos = start;          // read pointer position
        int part = 0;          // the current part (int, float and sci parts of the number)
        boolean neg = false;      // true if part is a negative number
        // the max long is 2147483647
        final int MAX_INT_BIT = 9;

        //ACSII
        //'0' : 48 
        //'9': 57
        //'-' 45
        //'.' 46
        //'e':101
        //'E':69
        while (f[pos] == ' ') {
            pos++;
        }
        // find start
        while (pos < end && (f[pos] < 48 || f[pos] > 57) && f[pos] != 45 && f[pos] != 46) {
            pos++;
        }

        // sign
        if (f[pos] == 45) {
            neg = true;
            pos++;
        }

        // integer part
        while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
            part = part * 10 + (f[pos++] - 48);
        }
        ret = neg ? (double) (part * -1) : (double) part;

        // float part
        if (pos < end && f[pos] == 46) {
            pos++;
            int mul = 1;
            part = 0;
            int num = 0;
            while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
                num++;
                if (num <= MAX_INT_BIT) {
                    part = part * 10 + (f[pos] - 48);
                    mul *= 10;
                }
                pos++;
            }
            ret = neg ? ret - (double) part / (double) mul : ret + (double) part / (double) mul;
        }

        // scientific part
        if (pos < end && (f[pos] == 101 || f[pos] == 69)) {
            pos++;
            if (pos < end) {
                neg = (f[pos] == 45);
                pos++;
                part = 0;
                while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
                    part = part * 10 + (f[pos++] - 48);
                }
                if (neg) {
                    ret = ret / Math.pow(10, part);
                } else {
                    ret = ret * Math.pow(10, part);
                }
            }
        }
        if (takeLog) {
            ret = Math.log10(ret);
        }
        return (float) ret;
    }

    private double parseDouble(byte[] f, int start, int end) {
        double ret = 0f;         // return value
        int pos = start;          // read pointer position
        int part = 0;          // the current part (int, float and sci parts of the number)
        boolean neg = false;      // true if part is a negative number
        // the max long is 2147483647
        final int MAX_INT_BIT = 9;

        //ACSII
        //'0' : 48 
        //'9': 57
        //'-' 45
        //'.' 46
        //'e':101
        //'E':69
        while (f[pos] == ' ') {
            pos++;
        }
        // find start
        while (pos < end && (f[pos] < 48 || f[pos] > 57) && f[pos] != 45 && f[pos] != 46) {
            pos++;
        }

        // sign
        if (f[pos] == 45) {
            neg = true;
            pos++;
        }

        // integer part
        while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
            part = part * 10 + (f[pos++] - 48);
        }
        ret = neg ? (double) (part * -1) : (double) part;

        // float part
        if (pos < end && f[pos] == 46) {
            pos++;
            int mul = 1;
            part = 0;
            int num = 0;
            while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
                num++;
                if (num <= MAX_INT_BIT) {
                    part = part * 10 + (f[pos] - 48);
                    mul *= 10;
                }
                pos++;
            }
            ret = neg ? ret - (double) part / (double) mul : ret + (double) part / (double) mul;
        }

        // scientific part
        if (pos < end && (f[pos] == 101 || f[pos] == 69)) {
            pos++;
            if (pos < end) {
                neg = (f[pos] == 45);
                pos++;
                part = 0;
                while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
                    part = part * 10 + (f[pos++] - 48);
                }
                if (neg) {
                    ret = ret / Math.pow(10, part);
                } else {
                    ret = ret * Math.pow(10, part);
                }
            }
        }
        return ret;
    }
    private int[] tempInt;

    private int[] tokenizeDelimiter(byte[] string, int startI, int endI, byte delimiter) {
        int tempLength = ((endI - startI) / 2) + 2;
        if (tempInt == null || tempInt.length < tempLength) {
            tempInt = new int[tempLength];
        }

        int i = startI;
        int j = 0;

        int wordCount = 0;
        int tt = startI;
        tempInt[wordCount] = startI;
        wordCount++;
        if (endI > string.length) {
            endI = string.length;
        }
        j = -1;
        for (; tt < endI; tt++) {
            if (string[tt] == delimiter) {
                j = tt;
                break;
            }
        }

        while (j >= 0 && j <= endI) {
            tempInt[wordCount] = j;
            wordCount++;
            i = j + 1;
            j = -1;
            for (tt = i; tt < endI; tt++) {
                if (string[tt] == delimiter) {
                    j = tt;
                    break;
                }
            }
        }

        if (i <= endI) {
            // tempInt[wordCount++] = i;
            tempInt[wordCount] = endI;
            wordCount++;
        }

        int[] result = new int[wordCount];
        System.arraycopy(tempInt, 0, result, 0, wordCount);
        return result;
    }

    private int tokenizeDelimiter(byte[] string, int startI, int endI, byte delimiter, int[] indexes) {
        int i = startI;
        int j = 0;
        int wordCount = 0;
        int tt = startI;
        j = -1;
        indexes[wordCount] = startI;
        wordCount++;
        if (endI > string.length) {
            endI = string.length;
        }
        for (; tt < endI; tt++) {
            if (string[tt] == delimiter) {
                j = tt;
                break;
            }
        }

        while (j >= 0 && j <= endI) {
            indexes[wordCount] = j;
            wordCount++;
            i = j + 1;
            j = -1;
            for (tt = i; tt < endI; tt++) {
                if (string[tt] == delimiter) {
                    j = tt;
                    break;
                }
            }
        }

        if (i <= endI) {
            // tempInt[wordCount++] = i;
            indexes[wordCount] = endI;
        }
        return wordCount + 1;
    }

}
