/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import cern.colt.list.IntArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;

import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFFont;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.cobi.kggseq.Constants;
import org.cobi.kggseq.Options;
import org.cobi.kggseq.entity.Chromosome;
import org.cobi.kggseq.entity.Genome;
import org.cobi.kggseq.entity.ProteinDomain;
import org.cobi.kggseq.entity.RefExon;
import org.cobi.kggseq.entity.RefGene;
import org.cobi.kggseq.entity.RefIntron;
import org.cobi.kggseq.entity.RefRefmRNA;
import org.cobi.kggseq.entity.ReferenceGenome;
import org.cobi.kggseq.entity.SeqSegment;
import org.cobi.kggseq.entity.RefmRNA;
import org.cobi.kggseq.entity.mRNA;
import org.cobi.util.download.stable.DownloadTaskEvent;
import org.cobi.util.download.stable.DownloadTaskListener;
import org.cobi.util.download.stable.HttpClient4API;
import org.cobi.util.download.stable.HttpClient4DownloadTask;
import org.cobi.util.file.LocalFileFunc;
import org.cobi.util.text.LocalExcelFile;
import org.cobi.util.text.Util;

/**
 *
 * @author mxli
 */
public class GeneRegionParser implements Constants {

    //  private static final Log LOG = Log.getInstance(GeneRegionParser.class);
    private static final Logger LOG = Logger.getLogger(GeneRegionParser.class);

    public void assignGeneSymb2Transcript(String vAFile, Genome genome) throws Exception {
        int indexmRNAName = 1;
        int indexName2 = 12;
        int maxColNum = indexmRNAName;
        maxColNum = Math.max(maxColNum, indexName2);
        String currentLine = null;
        int lineCounter = 0;

        File dataFile = new File(vAFile);
        BufferedReader br = null;

        String mRNAName;
        String geneSym;

        boolean incomplete;
        Map<String, String> mRNAGeneSymbMap = new HashMap<String, String>();

        StringBuilder tmpBuffer = new StringBuilder();
        try {
            br = LocalFileFunc.getBufferedReader(vAFile);

            while ((currentLine = br.readLine()) != null) {
                lineCounter++;
                StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles 
                mRNAName = null;
                geneSym = null;

                incomplete = false;
                for (int iCol = 0; iCol <= maxColNum; iCol++) {
                    if (st.hasMoreTokens()) {
                        tmpBuffer.delete(0, tmpBuffer.length());
                        tmpBuffer.append(st.nextToken().trim());
                        if (iCol == indexmRNAName) {
                            mRNAName = tmpBuffer.toString();
                        } else if (iCol == indexName2) {
                            geneSym = tmpBuffer.toString();
                        }
                    } else {
                        break;
                    }
                    if (iCol == maxColNum) {
                        incomplete = false;
                    }
                }

                if (incomplete) {
                    continue;
                }

                mRNAGeneSymbMap.put(mRNAName, geneSym);
                // System.out.println(currentLine);
            }

        } catch (NumberFormatException nex) {
            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
            // LOG.error(nex, info);
            throw new Exception(info);
        }
        br.close();

        Chromosome[] chroms = genome.getChromosomes();
        String geneSymb = null;

        for (int i = 0; i < chroms.length; i++) {
            if (chroms[i] == null) {
                continue;
            }

            for (mRNA mrna : chroms[i].mRNAList) {
                geneSymb = mRNAGeneSymbMap.get(mrna.refID);
                if (geneSymb != null) {
                    if (mrna.geneSymb != null) {
                        if (!geneSymb.equals(mrna.geneSymb)) {
                            LOG.info("Warning! Inconsisten gene symbols in the reference database ("
                                    + geneSymb + ") and input data (" + mrna.geneSymb + ") for " + mrna.refID + "! and the former will be used.");
                            mrna.geneSymb = geneSymb;
                        }
                    } else {
                        mrna.geneSymb = geneSymb;
                    }
                }
            }
        }
    }

    public Map<String, float[]> readRefGeneLength(String vAFile, int splicing) throws Exception {
        int indexmRNAName = 1;
        int indexChom = 2;
        int indexStrand = 3;
        int indexTxStart = 4;
        int indexTxEnd = 5;
        int indexCdsStart = 6;
        int indexCdsEnd = 7;
        int indexExonCount = 8;
        int indexExonStarts = 9;
        int indexExonEnds = 10;
        int indexName2 = 12;
        int indexSeq = 16;

        int maxColNum = indexmRNAName;
        maxColNum = Math.max(maxColNum, indexChom);
        maxColNum = Math.max(maxColNum, indexStrand);
        maxColNum = Math.max(maxColNum, indexTxStart);
        maxColNum = Math.max(maxColNum, indexTxEnd);
        maxColNum = Math.max(maxColNum, indexCdsStart);
        maxColNum = Math.max(maxColNum, indexCdsEnd);
        maxColNum = Math.max(maxColNum, indexExonCount);
        maxColNum = Math.max(maxColNum, indexExonStarts);
        maxColNum = Math.max(maxColNum, indexExonEnds);
        maxColNum = Math.max(maxColNum, indexName2);
        maxColNum = Math.max(maxColNum, indexSeq);

        String currentLine = null;
        String currChr = null;
        StringBuilder tmpBuffer = new StringBuilder();
        long lineCounter = 0;

        File dataFile = new File(vAFile);
        BufferedReader br = null;

        boolean incomplete = true;
        String mRNAName = null;
        char strand = '0';
        int cdsStart = -1;
        int cdsEnd = -1;
        int txStart = -1;
        int txEnd = -1;
        String exonStarts = null;
        String exonEnds = null;
        String geneSym = null;
        int geneNum = 0;
        int transcriptNum = 0;

        Map<String, float[]> geneSegAvgLen = new HashMap<String, float[]>();

        try {
            br = LocalFileFunc.getBufferedReader(vAFile);

            int UTR5Len = 0;
            int UTR3Len = 0;
            int exonLen = 0;
            int intronLen = 0;
            IntArrayList lens = null;

            Map<String, IntArrayList> mrnaSegLenAll = new HashMap<String, IntArrayList>();
            splicing = splicing * 2;
            int segNum = 0;
            float segNumf = 0.0f;
            float unit = 1000f;

            while ((currentLine = br.readLine()) != null) {
                lineCounter++;
                StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles
                incomplete = true;
                currChr = null;
                mRNAName = null;
                strand = '0';
                cdsStart = -1;
                cdsEnd = -1;
                txStart = -1;
                txEnd = -1;
                exonStarts = null;
                exonEnds = null;
                geneSym = null;

                for (int iCol = 0; iCol <= maxColNum; iCol++) {
                    if (st.hasMoreTokens()) {
                        tmpBuffer.delete(0, tmpBuffer.length());
                        tmpBuffer.append(st.nextToken().trim());
                        if (iCol == indexmRNAName) {
                            mRNAName = tmpBuffer.toString();
                        } else if (iCol == indexStrand) {
                            strand = tmpBuffer.charAt(0);
                        } else if (iCol == indexChom) {
                            currChr = tmpBuffer.toString();
                            currChr = currChr.substring(3);
                        } else if (iCol == indexTxStart) {
                            txStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexTxEnd) {
                            txEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsStart) {
                            cdsStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsEnd) {
                            cdsEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexExonStarts) {
                            exonStarts = tmpBuffer.toString();
                        } else if (iCol == indexExonEnds) {
                            exonEnds = tmpBuffer.toString();
                        } else if (iCol == indexName2) {
                            geneSym = tmpBuffer.toString();
                        }
                    } else {
                        break;
                    }
                    if (iCol == maxColNum) {
                        incomplete = false;
                    }
                }

                if (incomplete) {
                    continue;
                }

                if (cdsStart == cdsEnd) {
                    continue;
                }

                String[] bounderStartStr = exonStarts.split(",");
                String[] bounderEndStr = exonEnds.split(",");

                int[] bounderStarts = new int[bounderStartStr.length];
                int[] bounderEnds = new int[bounderEndStr.length];
                for (int i = 0; i < bounderStarts.length; i++) {
                    bounderStarts[i] = Util.parseInt(bounderStartStr[i]);
                    bounderEnds[i] = Util.parseInt(bounderEndStr[i]);
                }

                if (strand == '+') {
                    UTR5Len = 0;
                    UTR3Len = 0;
                    exonLen = 0;

                    //assume the boundary is not inclusive                     
                    for (int i = 0; i < bounderStarts.length; i++) {
                        exonLen += (bounderEnds[i] - bounderStarts[i]);

                        //sometimes the ccds do no start from the first exomes
                        if (cdsStart >= bounderEnds[i]) {
                            UTR5Len += (bounderEnds[i] - bounderStarts[i]);
                            exonLen -= (bounderEnds[i] - bounderStarts[i]);
                        } else if (cdsStart >= bounderStarts[i] && cdsStart <= bounderEnds[i]) {
                            UTR5Len += (cdsStart - bounderStarts[i]);
                            exonLen -= (cdsStart - bounderStarts[i]);
                        }

                        if (cdsEnd <= bounderStarts[i]) {
                            UTR3Len += (bounderEnds[i] - bounderStarts[i]);
                            exonLen -= (bounderEnds[i] - bounderStarts[i]);
                        } else if (cdsEnd >= bounderStarts[i] && cdsEnd <= bounderEnds[i]) {
                            UTR3Len += (bounderEnds[i] - cdsEnd);
                            exonLen -= (bounderEnds[i] - cdsEnd);
                        }
                    }
                    intronLen = 0;
                    for (int i = 1; i < bounderStarts.length; i++) {
                        intronLen += (bounderStarts[i] - bounderEnds[i - 1]);
                        //adjust for splicing lenght
                        intronLen -= splicing;
                        exonLen += splicing;
                    }
                    lens = mrnaSegLenAll.get(geneSym);
                    if (lens == null) {
                        lens = new IntArrayList();
                        mrnaSegLenAll.put(geneSym, lens);
                    }

                    lens.add(UTR5Len);
                    lens.add(UTR3Len);
                    lens.add(exonLen);
                    lens.add(intronLen);
                    mrnaSegLenAll.put(geneSym, lens);
                } else {
                    UTR5Len = 0;
                    UTR3Len = 0;
                    exonLen = 0;

                    //assume the boundary is not inclusive                     
                    for (int i = bounderStarts.length - 1; i >= 0; i--) {
                        exonLen += (bounderEnds[i] - bounderStarts[i]);

                        //sometimes the ccds do no start from the first exomes
                        if (cdsEnd <= bounderStarts[i]) {
                            UTR5Len += (bounderEnds[i] - bounderStarts[i]);
                            exonLen -= (bounderEnds[i] - bounderStarts[i]);
                        } else if (cdsEnd >= bounderStarts[i] && cdsEnd <= bounderEnds[i]) {
                            UTR5Len += (bounderEnds[i] - cdsEnd);
                            exonLen -= (bounderEnds[i] - cdsEnd);
                        }

                        if (cdsStart >= bounderEnds[i]) {
                            UTR3Len += (bounderEnds[i] - bounderStarts[i]);
                            exonLen -= (bounderEnds[i] - bounderStarts[i]);
                        } else if (cdsStart >= bounderStarts[i] && cdsStart <= bounderEnds[i]) {
                            UTR3Len += (cdsStart - bounderStarts[i]);
                            exonLen -= (cdsStart - bounderStarts[i]);
                        }
                    }

                    intronLen = 0;
                    for (int i = 1; i < bounderStarts.length; i++) {
                        intronLen += (bounderStarts[i] - bounderEnds[i - 1]);
                        //adjust for splicing lenght
                        intronLen -= splicing;
                        exonLen += splicing;
                    }
                    lens = mrnaSegLenAll.get(geneSym);
                    if (lens == null) {
                        lens = new IntArrayList();
                        mrnaSegLenAll.put(geneSym, lens);
                    }

                    lens.add(UTR5Len);
                    lens.add(UTR3Len);
                    lens.add(exonLen);
                    lens.add(intronLen);
                    mrnaSegLenAll.put(geneSym, lens);
                }

            }

            float UTR5LenT = 0;
            float UTR3LenT = 0;
            float exonLenT = 0;
            float intronLenT = 0;
            float UTR5Lenf = 0;
            float UTR3Lenf = 0;
            float exonLenf = 0;
            float intronLenf = 0;
            /*
             pp <-read.table("E:/home/mxli/MyJava/kggseq1/geneLen.txt", header=TRUE, na.strings = "NaN") ;            
             hist(pp$Exon, breaks=100);            
             * 
             */
            /*
             BufferedWriter bw = new BufferedWriter(new FileWriter("geneLen.txt"));
             bw.write("Gene\tExon\tUTR5\tUTR3\tintron\n");
             * 
             */

            //take the average length
            for (Map.Entry<String, IntArrayList> geneLens : mrnaSegLenAll.entrySet()) {
                segNum = geneLens.getValue().size() / 4;
                segNumf = segNum * unit;
                UTR5Lenf = 0;
                UTR3Lenf = 0;
                exonLenf = 0;
                intronLenf = 0;
                for (int i = 0; i < segNum; i++) {
                    UTR5Lenf += geneLens.getValue().getQuick(i * 4);
                    UTR3Lenf += geneLens.getValue().getQuick(i * 4 + 1);
                    exonLenf += geneLens.getValue().getQuick(i * 4 + 2);
                    intronLenf += geneLens.getValue().getQuick(i * 4 + 3);
                }

                UTR5Lenf /= segNumf;
                UTR3Lenf /= segNumf;
                exonLenf /= segNumf;
                intronLenf /= segNumf;

                if (UTR5Lenf < 0) {
                    UTR5Lenf = 0;
                }
                if (UTR3Lenf < 0) {
                    UTR3Lenf = 0;
                }
                if (exonLenf < 0) {
                    exonLenf = 0;
                }
                if (intronLenf < 0) {
                    intronLenf = 0;
                }

                float[] lenss = new float[]{exonLenf, UTR5Lenf, UTR3Lenf, intronLenf};
                geneSegAvgLen.put(geneLens.getKey(), lenss);
                // bw.write(geneLens.getKey() + "\t" + exonLenf + "\t" + UTR5Lenf + "\t" + UTR3Lenf + "\t" + intronLenf + "\n");
                UTR5LenT += UTR5Lenf;
                UTR3LenT += UTR3Lenf;
                exonLenT += exonLenf;
                intronLenT += intronLenf;
            }
            // bw.close();
            float[] lenss = new float[]{exonLenT, UTR5LenT, UTR3LenT, intronLenT};
            geneSegAvgLen.put("Total", lenss);
            //System.out.println("Total" + "\t" + exonLenT + "\t" + UTR5LenT + "\t" + UTR3LenT + "\t" + intronLenT);
        } catch (NumberFormatException nex) {
            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
            // LOG.error(nex, info);
            throw new Exception(info);
        }
        br.close();
        return geneSegAvgLen;
    }

    public Map<String, Double> readMergeRefGeneCodingLength(String[] vAFile, int splicing, boolean onlyCoding) throws Exception {
        Map<String, RefGene> mergedGeneMap = new HashMap<String, RefGene>();
        for (String fileStr : vAFile) {
            readDBGeneExons(fileStr, mergedGeneMap);
        }
        Map<String, Double> mergedGeneLen = new HashMap<String, Double>();

        for (Map.Entry<String, RefGene> geneLens : mergedGeneMap.entrySet()) {
            RefGene gene = geneLens.getValue();
            gene.combineExons(splicing, true);
            int totalMergedLen = 0;
            for (SeqSegment mgs : gene.getMergedSegments()) {
                totalMergedLen += (mgs.getEnd() - mgs.getStart() + 1);
            }
            gene.setTotalMergedLen(totalMergedLen);
            mergedGeneLen.put(geneLens.getKey(), geneLens.getValue().getTotalMergedLen() / 1000.0);
            //System.out.println(geneLens.getKey() + "\t" + geneLens.getValue().getTotalMergedLen() / 1000.0);
        }
        mergedGeneMap.clear();
        return mergedGeneLen;
    }

    public void readDBGeneExons(String vAFile, Map<String, RefGene> refGeneMap) throws Exception {
        int indexmRNAName = 1;
        int indexChom = 2;
        int indexStrand = 3;
        int indexTxStart = 4;
        int indexTxEnd = 5;
        int indexCdsStart = 6;
        int indexCdsEnd = 7;
        int indexExonCount = 8;
        int indexExonStarts = 9;
        int indexExonEnds = 10;
        int indexName2 = 12;
        //int indexSeq = 16;

        int maxColNum = indexmRNAName;
        maxColNum = Math.max(maxColNum, indexChom);
        maxColNum = Math.max(maxColNum, indexStrand);
        maxColNum = Math.max(maxColNum, indexTxStart);
        maxColNum = Math.max(maxColNum, indexTxEnd);
        maxColNum = Math.max(maxColNum, indexCdsStart);
        maxColNum = Math.max(maxColNum, indexCdsEnd);
        maxColNum = Math.max(maxColNum, indexExonCount);
        maxColNum = Math.max(maxColNum, indexExonStarts);
        maxColNum = Math.max(maxColNum, indexExonEnds);
        maxColNum = Math.max(maxColNum, indexName2);
        // maxColNum = Math.max(maxColNum, indexSeq);

        String currentLine = null;
        String currChr = null;
        StringBuilder tmpBuffer = new StringBuilder();
        long lineCounter = 0;

        File dataFile = new File(vAFile);
        BufferedReader br = null;

        boolean incomplete = true;
        String mRNAName = null;
        char strand = '0';
        int cdsStart = -1;
        int cdsEnd = -1;
        int txStart = -1;
        int txEnd = -1;
        String exonStarts = null;
        String exonEnds = null;
        String geneSym = null;

        try {
            br = LocalFileFunc.getBufferedReader(vAFile);

            RefGene gene = null;

            RefRefmRNA mran = null;
            RefExon exon = null;
            while ((currentLine = br.readLine()) != null) {
                lineCounter++;
                StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles
                incomplete = true;
                currChr = null;
                mRNAName = null;
                strand = '0';
                cdsStart = -1;
                cdsEnd = -1;
                txStart = -1;
                txEnd = -1;
                exonStarts = null;
                exonEnds = null;
                geneSym = null;

                for (int iCol = 0; iCol <= maxColNum; iCol++) {
                    if (st.hasMoreTokens()) {
                        tmpBuffer.delete(0, tmpBuffer.length());
                        tmpBuffer.append(st.nextToken().trim());
                        if (iCol == indexmRNAName) {
                            mRNAName = tmpBuffer.toString();
                        } else if (iCol == indexStrand) {
                            strand = tmpBuffer.charAt(0);
                        } else if (iCol == indexChom) {
                            if (tmpBuffer.charAt(0) == 'c') {
                                currChr = tmpBuffer.substring(3);
                            } else {
                                currChr = tmpBuffer.toString();
                            }
                        } else if (iCol == indexTxStart) {
                            txStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexTxEnd) {
                            txEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsStart) {
                            cdsStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsEnd) {
                            cdsEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexExonStarts) {
                            exonStarts = tmpBuffer.toString();
                        } else if (iCol == indexExonEnds) {
                            exonEnds = tmpBuffer.toString();
                        } else if (iCol == indexName2) {
                            geneSym = tmpBuffer.toString();
                        }
                    } else {
                        break;
                    }
                    if (iCol == maxColNum) {
                        incomplete = false;
                    }
                }

                if (incomplete) {
                    continue;
                }
//exclude non-coding genes
                if (cdsStart == cdsEnd) {
                    continue;
                }

                String[] bounderStartStr = exonStarts.split(",");
                String[] bounderEndStr = exonEnds.split(",");

                int[] bounderStarts = new int[bounderStartStr.length];
                int[] bounderEnds = new int[bounderEndStr.length];
                for (int i = 0; i < bounderStarts.length; i++) {
                    bounderStarts[i] = Util.parseInt(bounderStartStr[i]);
                    bounderEnds[i] = Util.parseInt(bounderEndStr[i]);
                }
                if (geneSym.startsWith("ENSG")) {
                    if (geneSym.endsWith("decay")) {
                        //System.out.println(geneSym);
                        continue;
                    }
                    String[] cells = geneSym.split(";");
                    geneSym = cells[1];
                }
                gene = refGeneMap.get(geneSym);
                if (gene == null) {
                    gene = new RefGene(geneSym);
                    gene.setChromosome(currChr);
                    gene.setStrand(strand);
                    refGeneMap.put(geneSym, gene);
                }

                mran = new RefRefmRNA(mRNAName, txStart - txEnd, null);
                gene.addmRNA(mran);

                if (strand == '+') {
                    //assume the boundary is not inclusive                     
                    for (int i = 0; i < bounderStarts.length; i++) {
                        exon = null;
                        //sometimes the ccds do no start from the first exomes
                        if (cdsStart >= bounderEnds[i]) {
                            exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i]);
                        } else if (cdsStart >= bounderStarts[i] && cdsStart <= bounderEnds[i]) {
                            if (cdsEnd >= bounderStarts[i] && cdsEnd <= bounderEnds[i]) {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i], cdsStart, cdsEnd, cdsEnd - cdsStart);
                            } else {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i], cdsStart, bounderEnds[i], bounderEnds[i] - cdsStart);
                            }
                        } else if (cdsStart <= bounderStarts[i]) {
                            if (cdsEnd <= bounderStarts[i]) {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i]);
                            } else if (cdsEnd >= bounderStarts[i] && cdsEnd <= bounderEnds[i]) {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i], bounderStarts[i], cdsEnd, cdsEnd - bounderStarts[i]);
                            } else {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i], bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i]);
                            }
                        }
                        mran.addExon(exon);

                        if (i < bounderStarts.length - 1) {
                            RefIntron intron = new RefIntron(bounderEnds[i], bounderStarts[i + 1], bounderStarts[i + 1] - bounderEnds[i]);
                            mran.addIntron(intron);
                        }
                    }

                } else {
                    //assume the boundary is not inclusive                     
                    for (int i = bounderStarts.length - 1; i >= 0; i--) {
                        //sometimes the ccds do no start from the first exomes
                        if (cdsEnd <= bounderStarts[i]) {
                            exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i]);
                        } else if (cdsEnd >= bounderStarts[i] && cdsEnd <= bounderEnds[i]) {
                            if (cdsStart >= bounderStarts[i] && cdsStart <= bounderEnds[i]) {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i], cdsStart, cdsEnd, cdsEnd - cdsStart);
                            } else {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i], bounderStarts[i], cdsEnd, cdsEnd - bounderStarts[i]);
                            }
                        } else {
                            if (cdsStart >= bounderEnds[i]) {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i]);
                            } else if (cdsStart >= bounderStarts[i] && cdsStart <= bounderEnds[i]) {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i], cdsStart, bounderEnds[i], bounderEnds[i] - cdsStart);
                            } else {
                                exon = new RefExon(bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i], bounderStarts[i], bounderEnds[i], bounderEnds[i] - bounderStarts[i]);
                            }
                        }
                        mran.addExon(exon);
                        if (i > 0) {
                            RefIntron intron = new RefIntron(bounderEnds[i - 1], bounderStarts[i], bounderStarts[i] - bounderEnds[i - 1]);
                            mran.addIntron(intron);
                        }
                    }
                }
            }
        } catch (NumberFormatException nex) {
            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
            // LOG.error(nex, info);
            throw new Exception(info);
        }
        br.close();

    }

    public void extractDBGenes(String vAFile, Set<String> geneSet) throws Exception {
        int indexmRNAName = 1;
        int indexChom = 2;
        int indexStrand = 3;
        int indexTxStart = 4;
        int indexTxEnd = 5;
        int indexCdsStart = 6;
        int indexCdsEnd = 7;
        int indexExonCount = 8;
        int indexExonStarts = 9;
        int indexExonEnds = 10;
        int indexName2 = 12;
        int indexSeq = 16;

        int maxColNum = indexmRNAName;
        maxColNum = Math.max(maxColNum, indexChom);
        maxColNum = Math.max(maxColNum, indexStrand);
        maxColNum = Math.max(maxColNum, indexTxStart);
        maxColNum = Math.max(maxColNum, indexTxEnd);
        maxColNum = Math.max(maxColNum, indexCdsStart);
        maxColNum = Math.max(maxColNum, indexCdsEnd);
        maxColNum = Math.max(maxColNum, indexExonCount);
        maxColNum = Math.max(maxColNum, indexExonStarts);
        maxColNum = Math.max(maxColNum, indexExonEnds);
        maxColNum = Math.max(maxColNum, indexName2);
        maxColNum = Math.max(maxColNum, indexSeq);

        String currentLine = null;
        String currChr = null;
        StringBuilder tmpBuffer = new StringBuilder();
        long lineCounter = 0;

        File dataFile = new File(vAFile);
        BufferedReader br = null;

        boolean incomplete = true;
        String mRNAName = null;
        char strand = '0';
        int cdsStart = -1;
        int cdsEnd = -1;
        int txStart = -1;
        int txEnd = -1;
        String exonStarts = null;
        String exonEnds = null;
        String geneSym = null;

        try {
            br = LocalFileFunc.getBufferedReader(vAFile);

            RefGene gene = null;

            RefRefmRNA mran = null;
            RefExon exon = null;
            while ((currentLine = br.readLine()) != null) {
                lineCounter++;
                StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles
                incomplete = true;
                currChr = null;
                mRNAName = null;
                strand = '0';
                cdsStart = -1;
                cdsEnd = -1;
                txStart = -1;
                txEnd = -1;
                exonStarts = null;
                exonEnds = null;
                geneSym = null;

                for (int iCol = 0; iCol <= maxColNum; iCol++) {
                    if (st.hasMoreTokens()) {
                        tmpBuffer.delete(0, tmpBuffer.length());
                        tmpBuffer.append(st.nextToken().trim());
                        if (iCol == indexmRNAName) {
                            mRNAName = tmpBuffer.toString();
                        } else if (iCol == indexStrand) {
                            strand = tmpBuffer.charAt(0);
                        } else if (iCol == indexChom) {
                            currChr = tmpBuffer.toString();
                            currChr = currChr.substring(3);
                        } else if (iCol == indexTxStart) {
                            txStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexTxEnd) {
                            txEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsStart) {
                            cdsStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsEnd) {
                            cdsEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexExonStarts) {
                            exonStarts = tmpBuffer.toString();
                        } else if (iCol == indexExonEnds) {
                            exonEnds = tmpBuffer.toString();
                        } else if (iCol == indexName2) {
                            geneSym = tmpBuffer.toString();
                        }
                    } else {
                        break;
                    }
                    if (iCol == maxColNum) {
                        incomplete = false;
                    }
                }

                if (incomplete) {
                    continue;
                }
//exclude non-coding genes
                if (cdsStart == cdsEnd) {
                    continue;
                }
                if (geneSet.contains(geneSym)) {
                    String[] cells = currentLine.split("\t");
                    System.out.print(cells[0]);
                    for (int t = 1; t < cells.length - 1; t++) {
                        System.out.print("\t");
                        System.out.print(cells[t]);
                    }
                    System.out.print("\n");
                }
                //klkl

            }
        } catch (NumberFormatException nex) {
            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
            // LOG.error(nex, info);
            throw new Exception(info);
        }
        br.close();

    }

    /**
     * @return @throws Exception
     * @pdOid f0621cff-9d97-421e-a77a-765bd0938dfb
     */
    public ReferenceGenome readRefGeneSeq(String vAFile, String dbLabel, int splicing, int nearGene, File domainFile, File idMapFile) throws Exception {
        int indexmRNAName = 1;
        int indexChom = 2;
        int indexStrand = 3;
        int indexTxStart = 4;
        int indexTxEnd = 5;
        int indexCdsStart = 6;
        int indexCdsEnd = 7;
        int indexExonCount = 8;
        int indexExonStarts = 9;
        int indexExonEnds = 10;
        int indexName2 = 12;
        int indexSeq = 16;
        int indexDelSite = 17;
        int indexInsSite = 18;
        int maxColNum = indexmRNAName;
        maxColNum = Math.max(maxColNum, indexChom);
        maxColNum = Math.max(maxColNum, indexStrand);
        maxColNum = Math.max(maxColNum, indexTxStart);
        maxColNum = Math.max(maxColNum, indexTxEnd);
        maxColNum = Math.max(maxColNum, indexCdsStart);
        maxColNum = Math.max(maxColNum, indexCdsEnd);
        maxColNum = Math.max(maxColNum, indexExonCount);
        maxColNum = Math.max(maxColNum, indexExonStarts);
        maxColNum = Math.max(maxColNum, indexExonEnds);
        maxColNum = Math.max(maxColNum, indexName2);
        maxColNum = Math.max(maxColNum, indexSeq);
        maxColNum = Math.max(maxColNum, indexDelSite);
        maxColNum = Math.max(maxColNum, indexInsSite);

        ReferenceGenome genome = new ReferenceGenome(splicing, nearGene, nearGene);
        genome.setName(dbLabel);
        String currentLine = null;
        String currChr = null;
        StringBuilder tmpBuffer = new StringBuilder();
        long lineCounter = 0;

        File dataFile = new File(vAFile);
        BufferedReader br = null;

        ProteinDomainRetriever pdR = new ProteinDomainRetriever();
        Map<String, String> refSeqUniportIDMap = null;
        Map<String, List<ProteinDomain>> uniportIDDomainMap = null;
        if (domainFile.exists()) {
            refSeqUniportIDMap = new HashMap<String, String>();
            uniportIDDomainMap = new HashMap<String, List<ProteinDomain>>();
            pdR.readProteinDomains(idMapFile.getCanonicalPath(), domainFile.getCanonicalPath(), refSeqUniportIDMap, uniportIDDomainMap);
        }

        boolean incomplete = true;
        String mRNAName = null;
        char strand = '0';
        int cdsStart = -1;
        int cdsEnd = -1;
        int txStart = -1;
        int txEnd = -1;
        String exonStarts = null;
        String exonEnds = null;
        String geneSym = null;
        int geneNum = 0;
        int transcriptNum = 0;
        int index = 0;
        Set<String> geneSymbSet = new HashSet<String>();
        String rnaSeq = null;
        Set<String> mRNAIDSet = new HashSet<String>();
        String delSite;
        String insSite;
        try {
            /*
             //skip to the head line 
             while ((currentLine = br.readLine()) != null) {
             lineCounter++;
             if (currentLine.startsWith("VAR")) {
             break;
             }
             }
             * 
             */

            br = LocalFileFunc.getBufferedReader(vAFile);

            while ((currentLine = br.readLine()) != null) {
                lineCounter++;
                StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles
                incomplete = true;
                currChr = null;
                mRNAName = null;
                strand = '0';
                cdsStart = -1;
                cdsEnd = -1;
                txStart = -1;
                txEnd = -1;
                exonStarts = null;
                exonEnds = null;
                geneSym = null;
                rnaSeq = null;
                delSite = null;
                insSite = null;

                for (int iCol = 0; iCol <= maxColNum; iCol++) {
                    if (st.hasMoreTokens()) {
                        tmpBuffer.delete(0, tmpBuffer.length());
                        tmpBuffer.append(st.nextToken().trim());
                        if (iCol == indexmRNAName) {
                            //remove the version id; it is usually not needed
                            index = tmpBuffer.indexOf(".");
                            if (index >= 0) {
                                mRNAName = tmpBuffer.substring(0, index);
                            } else {
                                mRNAName = tmpBuffer.toString();
                            }
                        } else if (iCol == indexStrand) {
                            strand = tmpBuffer.charAt(0);
                        } else if (iCol == indexChom) {
                            currChr = tmpBuffer.toString(); 
                            currChr = currChr.substring(3);
                        } else if (iCol == indexTxStart) {
                            txStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexTxEnd) {
                            txEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsStart) {
                            cdsStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsEnd) {
                            cdsEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexExonStarts) {
                            exonStarts = tmpBuffer.toString();
                        } else if (iCol == indexExonEnds) {
                            exonEnds = tmpBuffer.toString();
                        } else if (iCol == indexName2) {
                            geneSym = tmpBuffer.toString();
                        } else if (iCol == indexSeq) {
                            rnaSeq = tmpBuffer.toString();
                        } else if (iCol == indexDelSite) {
                            delSite = tmpBuffer.toString();
                        } else if (iCol == indexInsSite) {
                            insSite = tmpBuffer.toString();
                        }
                    } else {
                        break;
                    }

                }

                if (geneSym.endsWith("decay")) {
                    continue;
                }
                String[] infors = geneSym.split(";");
                if (infors.length > 1) {
                    geneSym = infors[1] + ":" + infors[0];
                } else {
                    geneSym = infors[0];
                }

                geneSymbSet.add(geneSym);
                RefmRNA mrna = new RefmRNA(mRNAName, txStart, txEnd, cdsStart, cdsEnd);
                mrna.setStrand(strand);
                if (refSeqUniportIDMap != null) {
                    String uniproID = refSeqUniportIDMap.get(mRNAName);
                    if (uniproID != null) {
                        mrna.setUniprotID(uniproID);
                        List<ProteinDomain> pdList = uniportIDDomainMap.get(uniproID);
                        if (pdList != null) {
                            mrna.setProteinDomainList(pdList);
                        }
                    }
                }

                String[] bounderStarts = exonStarts.split(",");
                String[] bounderEnds = exonEnds.split(",");
                // System.out.println(mRNAName);
                for (int i = 0; i < bounderStarts.length; i++) {
                    SeqSegment exon = new SeqSegment(Util.parseInt(bounderStarts[i]), Util.parseInt(bounderEnds[i]));
                    mrna.addExon(exon);
                }
                mrna.setmRnaSequenceStart(mrna.getExons().get(0).getStart());
                mrna.setmRnaSequence(rnaSeq.toUpperCase());
                if (delSite != null && !delSite.equals(".")) {
                    index = delSite.indexOf(":");
                    String[] delSites = delSite.substring(0, index).split(",");
                    int[] delSitesInt = new int[delSites.length];
                    for (int i = 0; i < delSitesInt.length; i++) {
                        delSitesInt[i] = Integer.parseInt(delSites[i]);
                    }
                    mrna.setDelSites(delSitesInt);
                    mrna.setDelSeq(delSite.substring(index + 1));
                }

                if (insSite != null && !insSite.equals(".")) {
                    index = insSite.indexOf(":");
                    String[] insSites = insSite.substring(0, index).split(",");
                    int[] insSitesInt = new int[insSites.length];
                    for (int i = 0; i < insSitesInt.length; i++) {
                        insSitesInt[i] = Integer.parseInt(insSites[i]);
                    }
                    mrna.setInsSites(insSitesInt);
                    mrna.setInsSeq(insSite.substring(index + 1));
                }

                mrna.makeAccuIntronLength();
                mrna.setGeneSymb(geneSym);
                int[] poss = genome.getmRNAPos(mrna.getRefID() + ":" + currChr + ":" + mrna.getCodingStart() + ":" + mrna.getCodingEnd());
                if (mRNAIDSet.contains(mrna.getRefID())) {
                    mrna.setMultipleMapping(true);
                } else {
                    mRNAIDSet.add(mrna.getRefID());
                }
                if (poss != null) {
                    mrna = genome.getmRNA(poss);
                    if (!currChr.equals(STAND_CHROM_NAMES[poss[0]])) {
                        //note a transcript can be mapped onto multiple locations
                        String info = "Duplicated refGene items: " + mRNAName;
                        LOG.info(info);
                        continue;
                    }
                } else {
                    transcriptNum++;
                    genome.addRefRNA(mrna, currChr);
                }
                // System.out.println(currentLine);
            }
            br.close();

            LOG.info(transcriptNum + " transcripts of " + geneSymbSet.size() + " genes have been read in the dataset " + dataFile.getName() + "!");

        } catch (NumberFormatException nex) {
            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
            // LOG.error(nex, info);
            throw new Exception(info);
        }

        genome.sortmRNAMakeIndexonChromosomes();
        return genome;
    }

    /**
     * @return @throws Exception
     * @pdOid f0621cff-9d97-421e-a77a-765bd0938dfb
     */
    public ReferenceGenome readRefGeneSeqUcsc(String refGeneFile, String refGeneFaFile, String dbLabel, int splicing, int nearGene, File domainFile, File idMapFile) throws Exception {
        int indexmRNAName = 1;
        int indexChom = 2;
        int indexStrand = 3;
        int indexTxStart = 4;
        int indexTxEnd = 5;
        int indexCdsStart = 6;
        int indexCdsEnd = 7;
        int indexExonCount = 8;
        int indexExonStarts = 9;
        int indexExonEnds = 10;
        int indexName2 = 12;

        int maxColNum = indexmRNAName;
        maxColNum = Math.max(maxColNum, indexChom);
        maxColNum = Math.max(maxColNum, indexStrand);
        maxColNum = Math.max(maxColNum, indexTxStart);
        maxColNum = Math.max(maxColNum, indexTxEnd);
        maxColNum = Math.max(maxColNum, indexCdsStart);
        maxColNum = Math.max(maxColNum, indexCdsEnd);
        maxColNum = Math.max(maxColNum, indexExonCount);
        maxColNum = Math.max(maxColNum, indexExonStarts);
        maxColNum = Math.max(maxColNum, indexExonEnds);
        maxColNum = Math.max(maxColNum, indexName2);

        ReferenceGenome genome = new ReferenceGenome(splicing, nearGene, nearGene);
        genome.setName(dbLabel);
        String currentLine = null;
        String currChr = null;
        StringBuilder tmpBuffer = new StringBuilder();
        long lineCounter = 0;

        BufferedReader br = LocalFileFunc.getBufferedReader(refGeneFaFile);
        ProteinDomainRetriever pdR = new ProteinDomainRetriever();
        Map<String, String> refSeqUniportIDMap = null;
        Map<String, List<ProteinDomain>> uniportIDDomainMap = null;
        if (domainFile.exists()) {
            refSeqUniportIDMap = new HashMap<String, String>();
            uniportIDDomainMap = new HashMap<String, List<ProteinDomain>>();
            pdR.readProteinDomains(idMapFile.getCanonicalPath(), domainFile.getCanonicalPath(), refSeqUniportIDMap, uniportIDDomainMap);
        }

        boolean incomplete = true;
        String mRNAName = null;
        char strand = '0';
        int cdsStart = -1;
        int cdsEnd = -1;
        int txStart = -1;
        int txEnd = -1;
        String exonStarts = null;
        String exonEnds = null;
        String geneSym = null;
        int geneNum = 0;
        int transcriptNum = 0;
        int index = 0;
        Set<String> geneSymbSet = new HashSet<String>();
        Set<String> mRNAIDSet = new HashSet<String>();

        try {
            /*
             //skip to the head line 
             while ((currentLine = br.readLine()) != null) {
             lineCounter++;
             if (currentLine.startsWith("VAR")) {
             break;
             }
             }
             * 
             */

            StringBuilder rnaSeq = new StringBuilder();

            // mrna.setmRnaSequence(rnaSeq.toUpperCase());
            String curMRNA = null;
            Map<String, String> rnaSeqMap = new HashMap<String, String>();
            while ((currentLine = br.readLine()) != null) {
                lineCounter++;
                currentLine = currentLine.trim();
                if (currentLine.isEmpty()) {
                    continue;
                }
                if (currentLine.startsWith(">")) {
                    if (curMRNA != null) {
                        rnaSeqMap.put(curMRNA, rnaSeq.toString().toUpperCase());
                    }
                    index = currentLine.indexOf(' ');
                    if (index > 0) {
                        curMRNA = currentLine.substring(1, index);
                    }
                    rnaSeq.delete(0, rnaSeq.length());
                } else {
                    rnaSeq.append(currentLine);
                }
                // System.out.println(currentLine);
            }
            if (curMRNA != null) {
                rnaSeqMap.put(curMRNA, rnaSeq.toString().toUpperCase());
                rnaSeq.delete(0, rnaSeq.length());
            }
            br.close();

            File dataFile = new File(refGeneFile);
            br = LocalFileFunc.getBufferedReader(refGeneFile);

            while ((currentLine = br.readLine()) != null) {
                lineCounter++;
                StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles
                incomplete = true;
                currChr = null;
                mRNAName = null;
                strand = '0';
                cdsStart = -1;
                cdsEnd = -1;
                txStart = -1;
                txEnd = -1;
                exonStarts = null;
                exonEnds = null;
                geneSym = null;

                for (int iCol = 0; iCol <= maxColNum; iCol++) {
                    if (st.hasMoreTokens()) {
                        tmpBuffer.delete(0, tmpBuffer.length());
                        tmpBuffer.append(st.nextToken().trim());
                        if (iCol == indexmRNAName) {
                            //remove the version id; it is usually not needed
                            index = tmpBuffer.indexOf(".");
                            if (index >= 0) {
                                mRNAName = tmpBuffer.substring(0, index);
                            } else {
                                mRNAName = tmpBuffer.toString();
                            }
                        } else if (iCol == indexStrand) {
                            strand = tmpBuffer.charAt(0);
                        } else if (iCol == indexChom) {
                            currChr = tmpBuffer.toString();
                            currChr = currChr.substring(3);
                        } else if (iCol == indexTxStart) {
                            txStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexTxEnd) {
                            txEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsStart) {
                            cdsStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsEnd) {
                            cdsEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexExonStarts) {
                            exonStarts = tmpBuffer.toString();
                        } else if (iCol == indexExonEnds) {
                            exonEnds = tmpBuffer.toString();
                        } else if (iCol == indexName2) {
                            geneSym = tmpBuffer.toString();
                        }
                    } else {
                        break;
                    }
                    if (iCol == maxColNum) {
                        incomplete = false;
                    }
                }

                if (incomplete) {
                    continue;
                }
                if (geneSym.endsWith("decay")) {
                    continue;
                }
                String[] infors = geneSym.split(";");
                if (infors.length > 1) {
                    geneSym = infors[1] + ":" + infors[0];
                } else {
                    geneSym = infors[0];
                }

                geneSymbSet.add(geneSym);
                RefmRNA mrna = new RefmRNA(mRNAName, txStart, txEnd, cdsStart, cdsEnd);
                mrna.setStrand(strand);
                if (refSeqUniportIDMap != null) {
                    String uniproID = refSeqUniportIDMap.get(mRNAName);
                    if (uniproID != null) {
                        mrna.setUniprotID(uniproID);
                        List<ProteinDomain> pdList = uniportIDDomainMap.get(uniproID);
                        if (pdList != null) {
                            mrna.setProteinDomainList(pdList);
                        }
                    }
                }

                String[] bounderStarts = exonStarts.split(",");
                String[] bounderEnds = exonEnds.split(",");
                // System.out.println(mRNAName);
                for (int i = 0; i < bounderStarts.length; i++) {
                    SeqSegment exon = new SeqSegment(Util.parseInt(bounderStarts[i]), Util.parseInt(bounderEnds[i]));
                    mrna.addExon(exon);
                }
                mrna.setmRnaSequenceStart(mrna.getExons().get(0).getStart());
                String seq = rnaSeqMap.get(mrna.getRefID());
                if (seq != null) {
                    mrna.setmRnaSequence(seq);
                }

                mrna.makeAccuIntronLength();
                mrna.setGeneSymb(geneSym);
                int[] poss = genome.getmRNAPos(mrna.getRefID() + ":" + currChr + ":" + mrna.getCodingStart() + ":" + mrna.getCodingEnd());
                if (mRNAIDSet.contains(mrna.getRefID())) {
                    mrna.setMultipleMapping(true);
                } else {
                    mRNAIDSet.add(mrna.getRefID());
                }
                if (poss != null) {
                    mrna = genome.getmRNA(poss);
                    if (!currChr.equals(STAND_CHROM_NAMES[poss[0]])) {
                        //note a transcript can be mapped onto multiple locations
                        String info = "Duplicated refGene items: " + mRNAName;
                        LOG.info(info);
                        continue;
                    }
                } else {
                    transcriptNum++;
                    //System.out.println(currChr);
                    genome.addRefRNA(mrna, currChr);

                }

                // System.out.println(currentLine);
            }
            br.close();

            LOG.info(transcriptNum + " transcripts of " + geneSymbSet.size() + " genes have been read in the dataset " + dataFile.getName() + "!");

        } catch (NumberFormatException nex) {
            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
            // LOG.error(nex, info);
            throw new Exception(info);
        }

        genome.sortmRNAMakeIndexonChromosomes();
        return genome;
    }

    /**
     * @return @throws Exception
     * @pdOid f0621cff-9d97-421e-a77a-765bd0938dfb
     */
    public Set<String> retrieveRefmRNA(String vAFile, String urlFolder, String resourceFolder) throws Exception {
        int indexmRNAName = 1;
        int indexChom = 2;
        int indexStrand = 3;
        int indexTxStart = 4;
        int indexTxEnd = 5;
        int indexCdsStart = 6;
        int indexCdsEnd = 7;
        int indexExonCount = 8;
        int indexExonStarts = 9;
        int indexExonEnds = 10;
        int indexName2 = 12;

        int maxColNum = indexmRNAName;
        maxColNum = Math.max(maxColNum, indexChom);
        maxColNum = Math.max(maxColNum, indexStrand);
        maxColNum = Math.max(maxColNum, indexTxStart);
        maxColNum = Math.max(maxColNum, indexTxEnd);
        maxColNum = Math.max(maxColNum, indexCdsStart);
        maxColNum = Math.max(maxColNum, indexCdsEnd);
        maxColNum = Math.max(maxColNum, indexExonCount);
        maxColNum = Math.max(maxColNum, indexExonStarts);
        maxColNum = Math.max(maxColNum, indexExonEnds);
        maxColNum = Math.max(maxColNum, indexName2);

        String currentLine = null;
        String currChr = null;
        StringBuilder tmpBuffer = new StringBuilder();
        long lineCounter = 0;

        File dataFile = new File(vAFile);
        BufferedReader br = LocalFileFunc.getBufferedReader(vAFile);

        boolean incomplete = true;
        System.out.print("Parse ");
        String mRNAName = null;
        char strand = '0';
        int cdsStart = -1;
        int cdsEnd = -1;
        int txStart = -1;
        int txEnd = -1;
        String exonStarts = null;
        String exonEnds = null;
        String geneSym = null;
        int geneNum = 0;
        int transcriptNum = 0;
        Set<String> chromNameSet = new HashSet<String>();

        try {
            while ((currentLine = br.readLine()) != null) {
                lineCounter++;
                StringTokenizer st = new StringTokenizer(currentLine.trim());
                //initialize varaibles
                incomplete = true;
                currChr = null;
                mRNAName = null;
                strand = '0';
                cdsStart = -1;
                cdsEnd = -1;
                txStart = -1;
                txEnd = -1;
                exonStarts = null;
                exonEnds = null;
                geneSym = null;
                for (int iCol = 0; iCol <= maxColNum; iCol++) {
                    if (st.hasMoreTokens()) {
                        tmpBuffer.delete(0, tmpBuffer.length());
                        tmpBuffer.append(st.nextToken().trim());
                        if (iCol == indexmRNAName) {
                            mRNAName = tmpBuffer.toString();
                        } else if (iCol == indexStrand) {
                            strand = tmpBuffer.charAt(0);
                        } else if (iCol == indexChom) {
                            currChr = tmpBuffer.toString();
                        } else if (iCol == indexTxStart) {
                            txStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexTxEnd) {
                            txEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsStart) {
                            cdsStart = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexCdsEnd) {
                            cdsEnd = Util.parseInt(tmpBuffer.toString());
                        } else if (iCol == indexExonStarts) {
                            exonStarts = tmpBuffer.toString();
                        } else if (iCol == indexExonEnds) {
                            exonEnds = tmpBuffer.toString();
                        } else if (iCol == indexName2) {
                            geneSym = tmpBuffer.toString();
                        }
                    } else {
                        break;
                    }
                    if (iCol == maxColNum) {
                        incomplete = false;
                    }
                }
                if (incomplete) {
                    continue;
                }
                chromNameSet.add(currChr);
            }

        } catch (NumberFormatException nex) {
            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
            // LOG.error(nex, info);
            throw new Exception(info);
        }
        br.close();

        int MAX_TASK = 2;
        boolean toDownload = false;
        ExecutorService exec = Executors.newFixedThreadPool(MAX_TASK);
        CompletionService serv = new ExecutorCompletionService(exec);
        int runningThread = 0;
        int i = 0;

        // GlobalManager.RESOURCE_PATH + "/" + options.refGenomeVersion + "/"
        long startTime = System.nanoTime();
        for (String chromName : chromNameSet) {
            File resourceFile = new File(resourceFolder + chromName + ".fa.gz");
            if (resourceFile.exists()) {
                long fileSize = resourceFile.length();
                long netFileLen = HttpClient4API.getContentLength(urlFolder + chromName + ".fa.gz");
                if (netFileLen <= 1024 || fileSize == netFileLen) {
                    continue;
                }

            } else {
                File parePath = resourceFile.getParentFile();
                if (!parePath.exists()) {
                    parePath.mkdirs();
                }
            }
            if (i == 0) {
                LOG.info("Downloading resource fasta from UCSC");
            }

            final HttpClient4DownloadTask task = new HttpClient4DownloadTask(urlFolder + chromName + ".fa.gz", 10);
            File filePath = new File(resourceFolder);
            if (!filePath.exists()) {
                filePath.mkdirs();
            }
            task.setLocalPath(resourceFile.getCanonicalPath());
            final String dbLabel = chromName;
            task.addTaskListener(new DownloadTaskListener() {

                @Override
                public void autoCallback(DownloadTaskEvent event) {
                    int progess = (int) (event.getTotalDownloadedCount() * 100.0 / event.getTotalCount());
                    String infor = progess + "%     Realtime Speed:" + event.getRealTimeSpeed() + " Global Speed:" + event.getGlobalSpeed();
                    System.out.print(infor);
                    char[] bs = new char[infor.length()];
                    Arrays.fill(bs, '\b');
                    System.out.print(bs);
                }

                @Override
                public void taskCompleted() throws Exception {
                    // File savedFile = new File(task.getLocalPath()); 

                    String msg1 = "Resource " + dbLabel + " has been downloaded!";
                    LOG.info(msg1);
                }
            });
            runningThread++;
            TimeUnit.MILLISECONDS.sleep(500);
            serv.submit(task);
            toDownload = true;
            i++;
        }

        for (int index = 0; index < runningThread; index++) {
            Future task = serv.take();
            String download = (String) task.get();
        }
        exec.shutdown();

        if (toDownload) {
            StringBuilder inforString = new StringBuilder();
            inforString.append("The lapsed time for downloading is : ");
            long endTime = System.nanoTime();
            inforString.append((endTime - startTime) / 1000000000.0);
            inforString.append(" Seconds.\n");
            LOG.info(inforString.toString());
        }
        return chromNameSet;
    }

    public void writeFullGeneFeatureSelectionReport(List<RefGene> geneList, String outFileName) throws Exception {
        XSSFWorkbook wb = new XSSFWorkbook();
        XSSFCellStyle style = wb.createCellStyle();
        //apply custom font to the text in the comment
        XSSFFont font = wb.createFont();

        List<String[]> tmpList = new ArrayList<String[]>();

        long smallerPos = 0;
        int listLen = geneList.size();

        XSSFSheet sheet3 = wb.createSheet("summary merged information ");
        // String[] titles3 = {"Gene_Symbol", "Gene_Description", "Pathway", "Chromosome", "Length",  "Segment_Number", "Total_Segment_Length", "Protein_Isoform_Number"};
        String[] titles3 = {"Gene_Symbol", "Gene_Description",
            "Segment_Number", "Total_Segment_Length", "Protein_Isoform_Number"};
        int colNum = titles3.length;
        tmpList.add(titles3);
        for (int geneCount = 0; geneCount < listLen; geneCount++) {
            RefGene gene = geneList.get(geneCount);
            String[] item = new String[colNum];
            item[0] = gene.getSymb();
            item[1] = gene.getDescription();
            // item[2] = gene.getPathwayInformation();
            // item[3] = gene.getChromosome();
            //  item[4] = String.valueOf(gene.getLength());
            List<SeqSegment> segList = gene.getMergedSegments();
            item[2] = String.valueOf(segList.size());
            int len = 0;

            for (int i = 0; i < segList.size(); i++) {
                SeqSegment seg = segList.get(i);
                len += (seg.getEnd() - seg.getStart()) + 1;
            }
            item[3] = String.valueOf(len);
            item[4] = String.valueOf(gene.getMRNAList().size());
            tmpList.add(item);
        }
        LocalExcelFile.writeArray2XLSXSheet(sheet3, wb, tmpList, true);
        tmpList.clear();

        XSSFSheet sheet2 = wb.createSheet("detailed mergered information ");
        String[] titles2 = {"Gene_Symbol", "Gene_Description", "Chromosome", "Length",
            "Segment_ID", "Segment_Feature", "Segment_Start", "Segment_End", "Segment_Length"};
        colNum = titles2.length;
        tmpList.add(titles2);
        for (int geneCount = 0; geneCount < listLen; geneCount++) {
            RefGene gene = geneList.get(geneCount);
            String[] item = new String[colNum];
            item[0] = gene.getSymb();
            item[1] = gene.getDescription();
            item[2] = gene.getChromosome();
            item[3] = String.valueOf(gene.getLength());
            List<SeqSegment> segList = gene.getMergedSegments();

            int segLen = 0;
            for (int i = 0; i < segList.size(); i++) {
                SeqSegment seg = segList.get(i);
                if (i > 0) {
                    item = new String[colNum];
                    item[0] = gene.getSymb();
                    item[2] = gene.getChromosome();
                }

                item[4] = String.valueOf(i + 1);
                item[5] = seg.getDescription();

                //assume the start is always smaller thatn the end
                item[6] = String.valueOf(seg.getStart());
                item[7] = String.valueOf(seg.getEnd());
                /*
                 if (gene.getStrand() == '-') {
                 item[6] = String.valueOf(gene.getEndPos() - seg.getStart() + 1);
                 item[7] = String.valueOf(gene.getEndPos() - seg.getEnd() + 1);
                 } else {
                 item[6] = String.valueOf(gene.getStartPos() + seg.getStart() - 1);
                 item[7] = String.valueOf(gene.getStartPos() + seg.getEnd() - 1);
                 }
                 */
                item[8] = String.valueOf(seg.getEnd() - seg.getStart() + 1);
                segLen += (seg.getEnd() - seg.getStart() + 1);
                tmpList.add(item);
            }
            // item = new String[colNum];
            //item[9] = String.valueOf(segLen);
        }
        LocalExcelFile.writeArray2XLSXSheet(sheet2, wb, tmpList, true);
        tmpList.clear();

        XSSFSheet sheet1 = wb.createSheet("original gene information");
        String[] titles1 = {"Gene_Symbol", "Chromosome", "Start_Position", "End_Position", "Length", "Transcript_ID",
            "Intron_Number", "Intron_Length",
            "Exon_ID", "Exon_Start_Position", "Exon_End_Position", "Exon_Length"
        };

        colNum = titles1.length;
        boolean firstGene = true;
        boolean firstMRNA = true;
        tmpList.add(titles1);
        smallerPos = 0;
        listLen = geneList.size();
        for (int geneCount = 0; geneCount < listLen; geneCount++) {
            RefGene gene = geneList.get(geneCount);
            firstGene = true;
            smallerPos = Math.min(gene.getStartPos(), gene.getEndPos());
            Map<String, RefRefmRNA> mRNAMap = gene.getMRNAList();
            for (Map.Entry<String, RefRefmRNA> m2 : mRNAMap.entrySet()) {
                firstMRNA = true;
                RefRefmRNA mrna = m2.getValue();
                List<RefExon> exonList = mrna.getExons();
                for (int i = 0; i < exonList.size(); i++) {
                    String[] item = new String[colNum];
                    if (firstGene) {
                        item[0] = gene.getSymb();
                        item[1] = gene.getChromosome();
                        item[2] = String.valueOf(gene.getStartPos());
                        item[3] = String.valueOf(gene.getEndPos());
                        item[4] = String.valueOf(gene.getLength());
                    }

                    if (firstMRNA) {
                        item[5] = mrna.getRefID();
                        item[6] = String.valueOf(mrna.getIntronsNum());
                        item[7] = String.valueOf(mrna.getIntronsLength());
                    }

                    RefExon exon = exonList.get(i);
                    item[8] = String.valueOf(i + 1);
                    item[9] = String.valueOf(exon.getStart() + smallerPos);
                    item[10] = String.valueOf(exon.getEnd() + smallerPos);
                    item[11] = String.valueOf(exon.getLength());
                    tmpList.add(item);
                    firstMRNA = false;
                    firstGene = false;
                }
            }
        }

        LocalExcelFile.writeArray2XLSXSheet(sheet1, wb, tmpList, true);
        tmpList.clear();

        // Write the output to a file
        if (!outFileName.endsWith(".xlsx")) {
            outFileName = outFileName + ".xlsx";
        }

        File outFile = new File(outFileName);
        FileOutputStream fileOut = new FileOutputStream(outFile);
        wb.write(fileOut);
        fileOut.close();
        LOG.info("The report is saved at " + outFile.getCanonicalPath());
    }

    public void readGeneSymbolFromFile(String fileName, Set<String> geneSet) throws Exception {
        File file = new File(fileName);
        if (!file.exists()) {
            return;
        }
        BufferedReader br1 = new BufferedReader(new FileReader(fileName));
        String line = null;
        int geneOder = 0;
        StringBuilder tmpStr = new StringBuilder();
        String geneSymb = null;
        while ((line = br1.readLine()) != null) {
            if (line.trim().length() == 0) {
                continue;
            }
            StringTokenizer tokenizer = new StringTokenizer(line);
            geneSymb = tmpStr.append(tokenizer.nextToken().trim()).toString();
            tmpStr.delete(0, tmpStr.length());
            geneSet.add(geneSymb);
        }

        LOG.info(geneSet.size() + " unique genes are read.");
        br1.close();
    }

    public static void main(String[] args) {
        GeneRegionParser gpp = new GeneRegionParser();
        try {

            GeneRegionParser parser = new GeneRegionParser();
            //parser.retrieveRefmRNA("resources/hg18/refGene.txt", "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/chromosomes/", "resources/hg18/fasta/");
            // parser.readRefGeneLength("resources/hg19/kggseq_hg19_refGene.txt.gz", 3);
            // Map<String, RefGene> refGeneMap = parser.readMergeRefGeneLength("resources/hg19/kggseq_hg19_refGene.txt.gz", 10,false);
            // Map<String, RefGene> refGeneMap = parser.readMergeRefGeneLength("resources/hg19/kggseq_hg19_GEncode.txt.gz", 10, false);
            // Map<String, RefGene> refGeneMap = parser.readMergeRefGeneLength("resources/hg19/kggseq_hg19_knownGene.txt.gz", 10, false);
            Map<String, RefGene> refGeneMap = new HashMap<String, RefGene>();
            //parser.readDBGeneExons("resources/hg19/kggseq_hg19_refGene.txt.gz", refGeneMap);
            // parser.readDBGeneExons("resources/hg19/kggseq_hg19_GEncode.txt.gz", refGeneMap);
            // parser.readDBGeneExons("resources/hg19/kggseq_hg19_knownGene.txt.gz", refGeneMap);
            //parser.readDBGeneExons("refseqs/ncbi_ref_GRCh37.p5_top_level.ucscformat", refGeneMap);
            //parser.readDBGeneExons("refseqs/refGene.txt", refGeneMap);
            parser.readDBGeneExons("refseqs/ensembl.ref.GRCh37.p3.toplevel.ucscformat1", refGeneMap);

            for (Map.Entry<String, RefGene> geneLens : refGeneMap.entrySet()) {
                RefGene gene = geneLens.getValue();

                gene.combineExons(0, true);
                int totalMergedLen = 0;
                for (SeqSegment mgs : gene.getMergedSegments()) {
                    totalMergedLen += (mgs.getEnd() - mgs.getStart() + 1);
                }
                gene.setTotalMergedLen(totalMergedLen);
            }

            //merge overlapped
            List<RefGene> reportGeneList = new ArrayList<RefGene>();
            Set<String> canidSet = new HashSet<String>();
            Set<String> validGeneSet = new HashSet<String>();
            //canidSet.add("STX16");
            //canidSet.add("JAK1");
            parser.readGeneSymbolFromFile("jk.txt", canidSet);
            // canidSet.add("TTN");

            // parser.extractDBGenes("resources/hg19/kggseq_hg19_refGene.txt.gz", canidSet);
            for (Map.Entry<String, RefGene> geneLens : refGeneMap.entrySet()) {
                if (canidSet != null && !canidSet.isEmpty() && !canidSet.contains(geneLens.getKey())) {
                    continue;
                }
                RefGene gene = geneLens.getValue();
                validGeneSet.add(geneLens.getKey());
                /// gene.combineExons(3, true);
                //gene.selectCombinedIntron(100);
                // bw.write(geneLens.getKey() + "\t" + exonLenf + "\t" + UTR5Lenf + "\t" + UTR3Lenf + "\t" + intronLenf + "\n");
                gene.getMergedSegments().get(0).setStart(gene.getMergedSegments().get(0).getStart() - 0);
                gene.getMergedSegments().get(gene.getMergedSegments().size() - 1).setEnd(gene.getMergedSegments().get(gene.getMergedSegments().size() - 1).getEnd() + 0);
                reportGeneList.add(gene);
            }
            canidSet.removeAll(validGeneSet);

            LOG.info("Ignored gene: " + canidSet.toString());
            // geneInforRetriever.writeBriefSegmentalReport(reportGeneList,option.outFileName);
            // geneInforRetriever.writeDetailedSegmentalReport(reportGeneList,option.outFileName);
            gpp.writeFullGeneFeatureSelectionReport(reportGeneList, "ensembl");

            // bw.close();
            //float[] lenss = new float[]{exonLenT, UTR5LenT, UTR3LenT, intronLenT};
            //mergedGeneLen.put("Total", lenss);
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }
}
