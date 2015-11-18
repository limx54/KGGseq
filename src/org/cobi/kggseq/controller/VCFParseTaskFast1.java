/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Output;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.zip.GZIPOutputStream;
import org.apache.log4j.Logger;
import org.cobi.kggseq.Constants;
import static org.cobi.kggseq.Constants.STAND_CHROM_NAMES;
import org.cobi.kggseq.GlobalManager;
import org.cobi.kggseq.entity.Genome;
import org.cobi.kggseq.entity.Individual;
import org.cobi.kggseq.entity.Variant;
import org.cobi.util.text.BGZFInputStream.BZSpider;
import org.cobi.util.text.BZIP2InputStream.Spider;

import org.cobi.util.text.Util;
import org.cobi.util.thread.Task;
import org.objenesis.strategy.StdInstantiatorStrategy;

/**
 *
 * @author mxli
 */
public class VCFParseTaskFast1 extends Task implements Callable<String>, Constants {

    private final Logger LOG = Logger.getLogger(VCFParseTaskFast1.class);
    private Map<String, Integer> chromNameIndexMap = new HashMap<String, Integer>();
    final private String UNKNOWN_CHROM_NAME0 = "Un";
    final private String UNKNOWN_CHROM_NAME1 = "GL";
    private List<Variant>[] varChroms = null;

    private String storagePath;
    private Kryo kryo = new Kryo();

    //As array is much faster than list; I try to use array when it does not was
    private int threadID = -1;
    private int[] caeSetID;
    private int[] controlSetID;

    private int indexCHROM;
    private int indexPOS;
    private int indexID;
    private int indexREF;
    private int indexALT;
    private int indexQUAL;
    private int indexFILTER;
    private int indexFORMAT;
    private int indexINFO;

    private double avgSeqQualityThrehsold;
    private double minMappingQual;
    private double maxStrandBias;
    private double maxFisherStrandBias;
    private double gtyQualityThrehsold;
    private int minGtySeqDepth;
    private int minSeqDepth;
    private double altAlleleFracRefHomThrehsold;
    private double altAlleleFractHetThrehsold;
    private double altAlleleFractAltHomThrehsold;
    private int minSecondPL;
    private double minBestGP;
    private int minOBS;
    private double sampleMafOver;
    private double sampleMafLess;
    private int maxGtyAlleleNum;
    private int[] tempInt;
    private Set<String> vcfLabelSet;

    private boolean noGtyVCF;
    private boolean considerSNP;
    private boolean considerIndel;
    private boolean needGty;
    private boolean needReadsInfor;
    private boolean needGtyQual;

    //result variables
    private int ignoredLowQualGtyNum = 0;
    private int ignoredLowDepthGtyNum = 0;
    private int ignoredBadAltFracGtyNum = 0;
    private int ignoredLowPLGtyNum = 0;
    private int ignoredLowGPGtyNum = 0;
    private int ignoreStrandBiasSBNum = 0;
    private int missingGtyNum = 0;
    private final int formatProbVarNum = 0;
    private int filterOutLowQualNum = 0;
    private int vcfFilterOutNum = 0;
    private int ignoredLineNumMinOBS = 0;
    private int ignoredLineNumMinMAF = 0;
    private int ignoredLineNumMaxMAF = 0;
    private int nonRSVariantNum = 0;
    private int ignoreMappingQualNum = 0;
    private int ignoreStrandBiasFSNum = 0;
    private int indelNum = 0, snvNum = 0;
    private int ignoredInproperChromNum = 0;
    private int ignoredVarBymaxGtyAlleleNum = 0;
    private int ignoredVarByRegionsInNum = 0;
    private int ignoredVarByRegionsOutNum = 0;

    private int totalVarNum = 0;
    private int totalAcceptVarNum = 0;
    //temp variables to save time
    private boolean checkVCFfilter = false;
    private double sampleMafOverC = 0;
    private double sampleMafLessC = 0;
    private int[] gtyQuality = null;
    private int[] gtys = null;
    private int[] gtyDepth = null;
    private int[] readCounts = null;

    private float[] readFractions = null;
    private int[] secondMostGtyPL = null;
    private int[] bestGtyGP = null;
    private boolean needAccoundAffect = false;
    private boolean needAccoundUnaffect = false;
    private boolean needMAFQCOver = false;
    private boolean needMAFQCLess = false;
    private double maf = 0;
    private boolean hasOrginalGenome = false;
    private int effectiveIndivNum = 0;
    private int totalPedSubjectNum = 0;
    private int maxVcfIndivNum = 0;
    private int maxEffectiveColVCF = -1;

    private int controlSize;
    private int caseSize;

    BZSpider bzSpider;
    int maxVarNum;

    //temp variables
    private int iGty = 0;
    private int index1 = 0;
    private int index2 = 0;

    private int t = 0;

    private int p = 0;

    private int index = 0;
    private int g11 = 0, g12 = 0, g22 = 0;
    private int indexA, indexB;
    private boolean isLowQualBreak = false;

    private boolean isIndel = false;
    private boolean isInvalid = false;
    private int alleleNum = 0;

    private int maxIndex2 = -1;
    private int ii = 0;
    private String tmpStr = null;

    private boolean isRegionInModel = false;
    private boolean isRegionOutModel = false;

    private StringBuilder tmpSB = new StringBuilder();
    private BufferedWriter[] vcfWriters;
    private int[][] regionsIn;
    private int[][] regionsOut;

    public void setRegionsIn(int[][] regions) {
        if (regions != null) {
            isRegionInModel = true;
            regionsIn = regions;
        }
    }

    public void setRegionsOut(int[][] regions) {
        if (regions != null) {
            isRegionOutModel = true;
            regionsOut = regions;
        }
    }

    public int getIgnoredVarByRegionsInNum() {
        return ignoredVarByRegionsInNum;
    }

    public int getIgnoredVarByRegionsOutNum() {
        return ignoredVarByRegionsOutNum;
    }

    public String getStoragePath() {
        return storagePath;
    }

    public void setStoragePath(String storagePath) {
        this.storagePath = storagePath;
    }

    public int getTotalVarNum() {
        return totalVarNum;
    }

    public int getTotalAcceptVarNum() {
        return totalAcceptVarNum;
    }

    public void setBr(BZSpider br) {
        this.bzSpider = br;
    }

    public void setMaxVarNum(int maxVarNum) {
        this.maxVarNum = maxVarNum;
    }

    public int getFormatProbVarNum() {
        return formatProbVarNum;
    }

    public int getFilterOutLowQualNum() {
        return filterOutLowQualNum;
    }

    public int getVcfFilterOutNum() {
        return vcfFilterOutNum;
    }

    public int getMissingGtyNum() {
        return missingGtyNum;
    }

    public int getIgnoredLowQualGtyNum() {
        return ignoredLowQualGtyNum;
    }

    public int getIgnoredLowDepthGtyNum() {
        return ignoredLowDepthGtyNum;
    }

    public int getIgnoredBadAltFracGtyNum() {
        return ignoredBadAltFracGtyNum;
    }

    public int getIgnoredLowPLGtyNum() {
        return ignoredLowPLGtyNum;
    }

    public int getIgnoredLowGPGtyNum() {
        return ignoredLowGPGtyNum;
    }

    public int getIgnoreStrandBiasSBNum() {
        return ignoreStrandBiasSBNum;
    }

    public int getIgnoredLineNumMinOBS() {
        return ignoredLineNumMinOBS;
    }

    public int getIgnoredLineNumMinMAF() {
        return ignoredLineNumMinMAF;
    }

    public int getIgnoredLineNumMaxMAF() {
        return ignoredLineNumMaxMAF;
    }

    public void setIgnoredLineNumMaxMAF(int ignoredLineNumMaxMAF) {
        this.ignoredLineNumMaxMAF = ignoredLineNumMaxMAF;
    }

    public int getIgnoredVarBymaxGtyAlleleNum() {
        return ignoredVarBymaxGtyAlleleNum;
    }

    public int getIgnoreMappingQualNum() {
        return ignoreMappingQualNum;
    }

    public int getIgnoreStrandBiasFSNum() {
        return ignoreStrandBiasFSNum;
    }

    public int getIndelNum() {
        return indelNum;
    }

    public int getSnvNum() {
        return snvNum;
    }

    public int getIgnoredInproperChromNum() {
        return ignoredInproperChromNum;
    }

    public void prepareTempVariables() {
        checkVCFfilter = false;
        if (vcfLabelSet != null) {
            checkVCFfilter = true;
        }
        sampleMafOverC = 1 - sampleMafOver;
        sampleMafLessC = 1 - sampleMafLess;
        needAccoundAffect = false;
        if (caeSetID != null) {
            needAccoundAffect = true;
            caseSize = caeSetID.length;
        }
        needAccoundUnaffect = false;
        if (controlSetID != null) {
            needAccoundUnaffect = true;
            controlSize = controlSetID.length;
        }
        if (orgGenome != null) {
            hasOrginalGenome = true;
        }

        if (sampleMafOver >= 0) {
            needMAFQCOver = true;
        }
        if (sampleMafLess < 0.5) {
            needMAFQCLess = true;
        }

        effectiveIndivNum = effectIndivIDInVCF.size();
        totalPedSubjectNum = pedVCFIDMap.length;
        maxEffectiveColVCF = indexCHROM;
        maxEffectiveColVCF = Math.max(maxEffectiveColVCF, indexPOS);
        maxEffectiveColVCF = Math.max(maxEffectiveColVCF, indexID);
        maxEffectiveColVCF = Math.max(maxEffectiveColVCF, indexREF);
        maxEffectiveColVCF = Math.max(maxEffectiveColVCF, indexALT);
        maxEffectiveColVCF = Math.max(maxEffectiveColVCF, indexQUAL);
        maxEffectiveColVCF = Math.max(maxEffectiveColVCF, indexFILTER);
        maxEffectiveColVCF = Math.max(maxEffectiveColVCF, indexFORMAT);
        maxEffectiveColVCF = Math.max(maxEffectiveColVCF, indexINFO);

        if (!noGtyVCF) {
            for (int i = 0; i < effectIndivIDInVCF.size(); i++) {
                if (maxEffectiveColVCF < effectIndivIDInVCF.getQuick(i) + indexFORMAT + 1) {
                    maxEffectiveColVCF = effectIndivIDInVCF.getQuick(i) + indexFORMAT + 1;
                    maxVcfIndivNum = effectIndivIDInVCF.getQuick(i);
                }
            }

            maxVcfIndivNum += 1;

            gtyQuality = new int[maxVcfIndivNum];
            gtys = new int[maxVcfIndivNum + maxVcfIndivNum];
            gtyDepth = new int[maxVcfIndivNum];
            //read counts 0, 1,2,3
            readCounts = new int[maxVcfIndivNum + maxVcfIndivNum];
            secondMostGtyPL = new int[maxVcfIndivNum];
            readFractions = new float[maxVcfIndivNum];
            bestGtyGP = new int[maxVcfIndivNum];

            Arrays.fill(readFractions, Float.NaN);

        }

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

    private int indexof(byte[] src, int start, int end, byte[] tar) {
        if (end - start < tar.length) {
            return -1;
        }
        boolean allMatched = false;
        end = end - tar.length;
        for (int i = start; i < end; i++) {
            allMatched = true;
            for (int j = 0; j < tar.length; j++) {
                if (tar[j] != src[i + j]) {
                    allMatched = false;
                    break;
                }
            }
            if (allMatched) {
                return i;
            }
        }
        return -1;
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
    /*
     public int[] tokenizeDelimiter(final byte[] string, final byte delimiter) {
     int strLen = string.length;
     int tempLength = (strLen / 2) + 2;
     if (tempInt == null || tempInt.length < tempLength) {
     tempInt = new int[tempLength];
     }

     int i = 0;
     int j = 0;
     int wordCount = 1;
     tempInt[wordCount] = 0;
     int tt = 0;

     j = -1;
     for (tt = 0; tt < strLen; tt++) {
     if (string[tt] == delimiter) {
     j = tt;
     break;
     }
     }

     while (j >= 0) {
     tempInt[++wordCount] = j;
     i = j + 1;
     j = -1;
     for (tt = i; tt < strLen; tt++) {
     if (string[tt] == delimiter) {
     j = tt;
     break;
     }
     }
     }

     if (i <= strLen) {
     tempInt[++wordCount] = i;
     }
     tempInt[wordCount] = strLen;
     int[] result = new int[wordCount];
     tempInt[0] = ++wordCount;
     System.arraycopy(tempInt, 1, result, 0, wordCount - 1);
     return result;
     }

   
     public void tokenizeDelimiter(final byte[] string, final byte delimiter, int[] fixedInt) {
     int strLen = string.length;
     int i = 0;
     int j = 0;
     int wordCount = 1;
     fixedInt[wordCount] = 0;
     int tt = 0;
     j = -1;
     for (tt = 0; tt < strLen; tt++) {
     if (string[tt] == delimiter) {
     j = tt;
     break;
     }
     }

     while (j >= 0) {
     fixedInt[++wordCount] = j;
     i = j + 1;
     j = -1;
     for (tt = i; tt < strLen; tt++) {
     if (string[tt] == delimiter) {
     j = tt;
     break;
     }
     }
     }

     if (i <= strLen) {
     fixedInt[++wordCount] = i;
     }
     fixedInt[wordCount] = strLen;
     fixedInt[0] = ++wordCount;
     }

     public int[] tokenizeDelimiter(byte[] string, int startI, int endI, byte delimiter, int maxIndex) {
     int tempLength = ((endI - startI) / 2) + 2;
     if (tempInt == null || tempInt.length < tempLength) {
     tempInt = new int[tempLength];

     }

     int i = startI;
     int j = 0;

     int wordCount = 0;

     int tt = startI;
     j = -1;
     for (; tt < endI; tt++) {
     if (string[tt] == delimiter) {
     j = tt;
     break;
     }
     }

     while (j >= 0 && j <= endI) {
     tempInt[wordCount] = j;
     if (wordCount >= maxIndex) {
     wordCount++;
     break;
     }
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

     if (wordCount <= maxIndex + 1) {
     if (i <= endI) {
     // tempInt[wordCount++] = i;
     tempInt[wordCount++] = i;
     }
     }

     tempInt[wordCount] = endI;
     wordCount++;
     int[] result = new int[wordCount];
     System.arraycopy(tempInt, 0, result, 0, wordCount);
     return result;
     }
 
     public static void tokenizeDelimiter(byte[] string, int startI, int endI, byte delimiter, int maxIndex, int[] indexes) {
     int i = startI;
     int j = 0;
     int wordCount = 0;

     int tt = startI;
     j = -1;
     for (; tt < endI; tt++) {
     if (string[tt] == delimiter) {
     j = tt;
     break;
     }
     }

     while (j >= 0 && j <= endI) {
     indexes[wordCount] = j;
     if (wordCount >= maxIndex) {
     wordCount++;
     break;
     }
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

     if (wordCount <= maxIndex + 1) {
     if (i <= endI) {
     // tempInt[wordCount++] = i;
     indexes[wordCount++] = i;
     }
     }
     indexes[wordCount] = endI;
     }
     */

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

    private float parseFloat(byte[] f, int start, int end) {
        float ret = 0f;         // return value
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
        ret = neg ? (float) (part * -1) : (float) part;

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
            ret = neg ? ret - (float) part / (float) mul : ret + (float) part / (float) mul;
        }

        // scientific part
        if (pos < end && (f[pos] == 101 || f[pos] == 69)) {
            pos++;
            neg = (f[pos] == 45);
            pos++;
            part = 0;
            while (pos < end && !(f[pos] > 57 || f[pos] < 48)) {
                part = part * 10 + (f[pos++] - 48);
            }
            if (neg) {
                ret = ret / (float) Math.pow(10, part);
            } else {
                ret = ret * (float) Math.pow(10, part);
            }
        }
        return ret;
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

    int gtyIndexInInfor = -1;
    int gtyQualIndexInInfor = -1;
    int gtyDepthIndexInInfor = -1;
    int gtyAlleleDepthIndexInInfor = -1;
    int gtyAltAlleleFracIndexInInfor = -1;
    int gtyPLIndexInInfor = -1;
    int gtyGPIndexInInfor = -1;

    boolean hasIndexGT = false;
    boolean hasIndexGQ = false;
    boolean hasIndexDP = false;
    boolean hasIndexAD = false;
    boolean hasIndexFA = false;
    boolean hasIndexPL = false;
    boolean hasIndexGP = false;

    byte[] gtBytes = "GT".getBytes();
    byte[] gqBytes = "GQ".getBytes();
    byte[] dpBytes = "DP".getBytes();
    byte[] adBytes = "AD".getBytes();
    byte[] faBytes = "FA".getBytes();
    byte[] plBytes = "PL".getBytes();
    byte[] gpBytes = "GP".getBytes();

    public int parseVariantsInFileOnlyFastToken() {
        byte[] currentLine = null;
        int[] cellDelimts = new int[maxEffectiveColVCF + 2];
        int acceptVarNum = 0;

        byte[] sByte;
        byte[] tmpByte = new byte[100];
        int lastTmpByteIndex = 0;

        int obsS;
        String ref = null;
        String alt = null;

        String currChr = null;
        int makerPostion = 0;
        String varLabel = null;
        String vcfLabel = null;
        boolean incomplete = true;

        double avgSeqQuality = 0;
        double mappingQual;
        double strandBias;
        try {

            //at most use 5 bits represent a genotype
        /*
             2 bits for an unphased -genotype of a bi-allelic sequence variants
             3 bits for a phased -genotype of a bi-allelic sequence variants
             3 bits for an unphased -genotype of a tri-allelic sequence variants
             4 bits for a phased -genotype of a tri-allelic sequence variants
             4 bits for an unphased -genotype of a quad-allelic sequence variants
             5 bits for a phased -genotype of a quad-allelic sequence variants        
             */
            /*
             if (caseIDSet.isEmpty()) {
             throw new Exception("It seems that you have no specified patients (labeled with \'2\') in your sample!"
             + " You need patient samples to proceed on KGGSeq!");
             }
             * 
             */
            //unfornatuely, mine splitter is faster than guava
            //Splitter niceCommaSplitter = Splitter.on('\t').limit(maxEffectiveColVCF + 1);
            //String[] cellDelimts = new String[maxEffectiveColVCF + 1];
            String inforS = null;

            int[] cellsBT = new int[maxGtyAlleleNum + 1];

            int[] cellsPL = new int[maxGtyAlleleNum * maxGtyAlleleNum + 1];
            int len;
            //warning if all of the genotypes are missing, it will have a problem.
            //decide the whether genotypes are phased or not //at most consider 3 alternative alleles
            //note in the current version of GATK, the PGT are given as alternative
            int[] formatCellIndexes;

            /*
             if (indexof(currentLine, 0, currentLine.length, ("0|0").getBytes()) >= 0 || indexof(currentLine, 0, currentLine.length, ("0|1").getBytes()) >= 0 || indexof(currentLine, 0, currentLine.length, ("1|0").getBytes()) >= 0 || indexof(currentLine, 0, currentLine.length, ("1|1").getBytes()) >= 0) {
             isPhased = true;
             }*/
            //isPhased = false; 
            int s = 0;
            int sc = 0;
            int subjectNum;
            boolean hasNotCheckedPhase = true;
            Integer chromID;
            boolean isWithinRegion;
            int vcfID;
            int regionNum;
            while ((currentLine = bzSpider.readLine(cellDelimts)) != null) {
                totalVarNum++;
                // System.out.println(currentLine); 
                if (cellDelimts[0] == 0) {
                    continue;
                }
                if (tmpByte.length <= currentLine.length) {
                    //the currentLine does not have line ending symbol
                    tmpByte = new byte[currentLine.length + 1];
                }
                lastTmpByteIndex = 0;
                //indexCHROM==0;
                sByte = Arrays.copyOfRange(currentLine, 0, cellDelimts[indexCHROM + 1]);
                currChr = new String(sByte);
                if (currChr.startsWith(UNKNOWN_CHROM_NAME0) || currChr.startsWith(UNKNOWN_CHROM_NAME1)) {
                    ignoredInproperChromNum++;
                    break;
                }

                if (currChr.charAt(0) == 'c' || currChr.charAt(0) == 'C') {
                    currChr = currChr.substring(3);
                }
                //Mitochondrion
                if (currChr.charAt(0) == 'M' || currChr.charAt(0) == 'm') {
                    currChr = "M";
                }

                //indexPOS=1;
                makerPostion = parseInt(currentLine, cellDelimts[indexPOS] + 1, cellDelimts[indexPOS + 1]);
                chromID = chromNameIndexMap.get(currChr);
                if (chromID == null) {
                    //System.err.println("Unrecognized chromosome name: " + currChr);
                    continue;
                }

                if (isRegionInModel) {
                    int[] pos = regionsIn[chromID];
                    if (pos == null) {
                        ignoredVarByRegionsInNum++;
                        continue;
                    }
                    isWithinRegion = false;
                    regionNum = pos.length / 2;
                    for (int k = 0; k < regionNum; k++) {
                        if (makerPostion >= pos[k * 2] && makerPostion <= pos[k * 2 + 1]) {
                            isWithinRegion = true;
                            break;
                        }
                    }

                    if (!isWithinRegion) {
                        ignoredVarByRegionsInNum++;
                        continue;
                    }
                }

                if (isRegionOutModel) {
                    int[] pos = regionsOut[chromID];
                    if (pos != null) {
                        isWithinRegion = false;
                        regionNum = pos.length / 2;
                        for (int k = 0; k < regionNum; k++) {
                            if (makerPostion >= pos[k * 2] && makerPostion <= pos[k * 2 + 1]) {
                                isWithinRegion = true;
                                break;
                            }
                        }
                        if (isWithinRegion) {
                            ignoredVarByRegionsOutNum++;
                            continue;
                        }
                    }
                }

                /*
                 if (makerPostion == 56102470) {
                 int sss = 0;
                 }*/
                //indexID=2
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexID] + 1, cellDelimts[indexID + 1]);
                varLabel = new String(sByte);
                if (varLabel.charAt(0) != 'r') {
                    nonRSVariantNum++;
                }
                //indexREF=3
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexREF] + 1, cellDelimts[indexREF + 1]);
                ref = new String(sByte);
                //alt = currentLine.substring(cellDelimts[3] + 1, cellDelimts[4]);
                //indexALT=4
                //for reference data, sometimes we do not have alternative alleles
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexALT] + 1, cellDelimts[indexALT + 1]);
                alt = new String(sByte);

                if (alt.charAt(0) == '.') {
                    continue;
                }
                isInvalid = false;
                String[] altAllelesIndexes = Util.tokenize(alt, ',');
                alleleNum = altAllelesIndexes.length + 1;
                if (alleleNum > maxGtyAlleleNum) {
                    ignoredVarBymaxGtyAlleleNum++;
                    continue;
                }
                isIndel = false;
                for (int ss = 0; ss < altAllelesIndexes.length; ss++) {
                    alt = altAllelesIndexes[ss];
                    //only one alternative alleles; the most common  scenario
                    if (ref.length() == alt.length()) {
                        //substitution
                        //now it can sonsider double polymorphsom
                        altAllelesIndexes[ss] = alt;
                    } else if (ref.length() < alt.length()) {
                        //insertion
                                /*examples 
                         insertion1
                         chr1 1900106 . TCT TCTCCT 217 . INDEL;DP=62;AF1=0.5;CI95=0.5,0.5;DP4=17,9,18,12;MQ=60;FQ=217;PV4=0.78,1,1,0.1 GT:PL:DP:SP:GQ 0/1:255,0,255:56:-991149567:99
                        
                         insertion2
                         chr1 109883576 . C CAT 214 . INDEL;DP=15;AF1=1;CI95=1,1;DP4=0,0,1,11;MQ=60;FQ=-70.5 GT:PL:DP:SP:GQ 1/1:255,36,0:12:-991149568:69
                         * 
                         */
                        //for Indel TTCC TT--
                        //for Insertion T +TTTT
                        tmpSB.delete(0, tmpSB.length());
                        if (alt.startsWith(ref)) {
                            for (t = ref.length(); t > 0; t--) {
                                tmpSB.append('+');
                            }
                            tmpSB.append(alt.substring(ref.length()));
                            altAllelesIndexes[ss] = tmpSB.toString();
                        } else if (alt.endsWith(ref)) {
                            tmpSB.append(alt.substring(0, alt.length() - ref.length()));
                            for (t = ref.length(); t > 0; t--) {
                                tmpSB.append('+');
                            }
                            altAllelesIndexes[ss] = tmpSB.toString();
                        }
                        isIndel = true;
                    } else if (ref.length() > alt.length()) {
                        //deletion     
                                /*examples
                         deletion1
                         chr1 113659065 . ACTCT ACT 214 . INDEL;DP=61;AF1=1;CI95=1,1;DP4=0,0,22,34;MQ=60;FQ=-204 GT:PL:DP:SP:GQ 1/1:255,169,0:56:-991149568:99
                         deletion2
                         chr1 1289367 . CTG C 101 . INDEL;DP=14;AF1=0.5;CI95=0.5,0.5;DP4=5,2,5,1;MQ=60;FQ=104;PV4=1,0.4,1,1 GT:PL:DP:SP:GQ 0/1:139,0,168:13:-991149568:99
                         */
                        //Note it can work for multiple deletion alleles like:chr1	158164305	.	TAA	TA,T

                        //for Indel TTCC TT--
                        //for Insertion T +TTTT
                        tmpSB.delete(0, tmpSB.length());
                        if (ref.startsWith(alt)) {
                            tmpSB.append(alt);
                            for (t = ref.length() - alt.length(); t > 0; t--) {
                                tmpSB.append('-');
                            }
                            altAllelesIndexes[ss] = tmpSB.toString();
                        } else if (ref.endsWith(alt)) {
                            for (t = ref.length() - alt.length(); t > 0; t--) {
                                tmpSB.append('-');
                            }
                            tmpSB.append(alt);
                            altAllelesIndexes[ss] = tmpSB.toString();
                        }
                        isIndel = true;
                    } else {
                        StringBuilder info = new StringBuilder("Unexpected (REF	ALT) format when parsing line :" + currentLine);
                        LOG.warn(info);
                        isInvalid = true;
                        // throw new Exception(info.toString());
                    }
                }
                if (isInvalid) {
                    continue;
                }

                if (isIndel) {
                    indelNum++;
                } else {
                    snvNum++;
                }

                if (!considerSNP || !considerIndel) {
                    //a lazy point 
                    incomplete = true;

                    //only consider Indel
                    if (!considerSNP && isIndel) {
                        incomplete = false;
                    } else if (!considerIndel && !isIndel) {
                        incomplete = false;
                    }

                    if (incomplete) {
                        continue;
                    }
                }

                /*
                 if (hasOrginalGenome) {
                 Variant[] vars = orgGenome.lookupVariants(currChr, makerPostion, isIndel, ref, altAllelesIndexes);
                 if (vars == null) {
                 continue;
                 }
                 }*/
                //initialize varaibles
                incomplete = true;
                obsS = 0;
                hasIndexGT = false;
                hasIndexGQ = false;
                hasIndexDP = false;
                hasIndexAD = false;
                hasIndexFA = false;
                isLowQualBreak = false;

                hasIndexPL = false;
                hasIndexGP = false;
                gtyPLIndexInInfor = -1;
                gtyGPIndexInInfor = -1;
                gtyIndexInInfor = -1;
                gtyQualIndexInInfor = -1;
                gtyDepthIndexInInfor = -1;
                gtyAlleleDepthIndexInInfor = -1;
                gtyAltAlleleFracIndexInInfor = -1;
                mappingQual = Integer.MAX_VALUE;
                strandBias = Integer.MIN_VALUE;

                inforS = null;

//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  
//chr1    109     .       A       T       237.97  PASS    AC=21;AF=0.328;AN=64;DP=47;Dels=0.02;HRun=0;HaplotypeScore=1.9147;MQ=44.81;MQ0=48;QD=5.53;SB=-28.76;sumGLbyD=9.00       GT:AD:DP:GQ:PL  0/1:6,1:3:15.67:16,0,64 0/0:3,0:1:3.01:0,3,33 
//chr1	53598	.	CCTA	C	447.88	PASS	AC=2;AF=1.00;AN=2;DP=0;Dels=0.50;HRun=0;HaplotypeScore=0.0000;MQ=20.50;MQ0=5;QD=40.72;SB=-138.61;sumGLbyD=46.54	GT:AD:DP:GQ:PL	./.	1/1:2,1:0:3.01:66,3,0
                sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexFILTER] + 1, cellDelimts[indexFILTER + 1]);
                vcfLabel = new String(sByte);
                if (checkVCFfilter && cellDelimts[indexFILTER + 1] - cellDelimts[indexFILTER] > 1 && !vcfLabelSet.contains(vcfLabel)) {
                    vcfFilterOutNum++;
                    isLowQualBreak = true;
                    continue;
                }

                if (cellDelimts[0] > indexQUAL && avgSeqQualityThrehsold > 0) {
                    // seqQualS = cellDelimts[indexQUAL];
                    if (currentLine[cellDelimts[indexQUAL] + 1] != 46) {
                        avgSeqQuality = parseFloat(currentLine, cellDelimts[indexQUAL] + 1, cellDelimts[indexQUAL + 1]);
                    } else {
                        //sometimes . denotes for ignored sequence qaulity information
                        avgSeqQuality = Integer.MAX_VALUE;
                    }

                    if (avgSeqQuality < avgSeqQualityThrehsold) {
                        filterOutLowQualNum++;
                        isLowQualBreak = true;

                        continue;
                    }
                }

                //I do notknow the minimun threshold for this 
                if (cellDelimts[0] > indexINFO) {
                    sByte = Arrays.copyOfRange(currentLine, cellDelimts[indexINFO] + 1, cellDelimts[indexINFO + 1]);
                    inforS = new String(sByte);
                    int inforSLen = inforS.length();
                    if (minMappingQual > 0) {
                        index1 = inforS.indexOf("MQ=");
                        if (index1 >= 0) {
                            index1 += 3;
                            index2 = index1 + 1;
                            while (index2 < inforSLen && inforS.charAt(index2) != ';') {
                                index2++;
                            }
                            if (Util.isNumeric(inforS.substring(index1, index2))) {
                                mappingQual = Util.parseInt(inforS.substring(index1, index2));
                            }
                            if (mappingQual < minMappingQual) {
                                ignoreMappingQualNum++;
                                isLowQualBreak = true;
                                continue;
                            }
                        }
                    }
                    if (maxStrandBias > 0) {
                        index1 = inforS.indexOf("SB=");
                        if (index1 >= 0) {
                            index1 += 3;
                            index2 = index1 + 1;
                            while (index2 < inforSLen && inforS.charAt(index2) != ';') {
                                index2++;
                            }
                            if (Util.isNumeric(inforS.substring(index1, index2))) {
                                strandBias = Util.parseFloat(inforS.substring(index1, index2));
                            }

                            if (strandBias > maxStrandBias) {
                                ignoreStrandBiasSBNum++;
                                isLowQualBreak = true;
                                continue;
                            }
                        }
                    }

                    if (maxFisherStrandBias > 0) {
                        index1 = inforS.indexOf("FS=");
                        if (index1 >= 0) {
                            index1 += 3;
                            index2 = index1 + 1;
                            while (index2 < inforSLen && inforS.charAt(index2) != ';') {
                                index2++;
                            }
                            if (Util.isNumeric(inforS.substring(index1, index2))) {
                                strandBias = Util.parseFloat(inforS.substring(index1, index2));
                            }

                            if (strandBias > maxFisherStrandBias) {
                                ignoreStrandBiasFSNum++;
                                isLowQualBreak = true;
                                continue;
                            }
                        }
                    }
                }

                //it is no gty vcf
                if (indexFORMAT < 0 || noGtyVCF) {
                    Variant var = new Variant(makerPostion, ref, altAllelesIndexes);
                    var.setIsIndel(isIndel);
                    var.setLabel(varLabel);
                    varChroms[chromID].add(var);

                    acceptVarNum++;
                    if (acceptVarNum >= maxVarNum) {
                        totalAcceptVarNum += acceptVarNum;
                        writeChromosomeToDiskClean();
                        acceptVarNum = 0;
                    }
                    continue;
                }

                if (maxVcfIndivNum > 0) {
                    Arrays.fill(gtys, -1);
                }

                //System.out.println(currentLine);
                //StringTokenizer st1 = new StringTokenizer(tmpBuffer.toString(), ":"); 
                formatCellIndexes = tokenizeDelimiter(currentLine, cellDelimts[indexFORMAT] + 1, cellDelimts[indexFORMAT + 1], (byte) ':');

                // String[] cells1 = Util.tokenize(cellDelimts[indexFORMAT], ':');
                ii = 0;
                maxIndex2 = 0;
                len = formatCellIndexes.length - 1;
                for (ii = 0; ii < len; ii++) {
                    if (equal(currentLine, ii == 0 ? formatCellIndexes[ii] : formatCellIndexes[ii] + 1, formatCellIndexes[ii + 1], gtBytes)) {
                        gtyIndexInInfor = ii;
                        hasIndexGT = true;
                        maxIndex2 = ii;
                    } else if (equal(currentLine, ii == 0 ? formatCellIndexes[ii] : formatCellIndexes[ii] + 1, formatCellIndexes[ii + 1], gqBytes)) {
                        if (gtyQualityThrehsold > 0) {
                            gtyQualIndexInInfor = ii;
                            hasIndexGQ = true;
                            maxIndex2 = ii;
                            Arrays.fill(gtyQuality, 0);
                        }
                    } else if (equal(currentLine, ii == 0 ? formatCellIndexes[ii] : formatCellIndexes[ii] + 1, formatCellIndexes[ii + 1], dpBytes)) {
                        if (minGtySeqDepth > 0) {
                            gtyDepthIndexInInfor = ii;
                            hasIndexDP = true;
                            maxIndex2 = ii;
                            Arrays.fill(gtyDepth, 0);
                        }
                    } else if (equal(currentLine, ii == 0 ? formatCellIndexes[ii] : formatCellIndexes[ii] + 1, formatCellIndexes[ii + 1], adBytes)) {
                        if (altAlleleFracRefHomThrehsold < 1 || altAlleleFractHetThrehsold > 0 || altAlleleFractAltHomThrehsold > 0 || needReadsInfor) {
                            gtyAlleleDepthIndexInInfor = ii;
                            hasIndexAD = true;
                            maxIndex2 = ii;
                            Arrays.fill(readCounts, -1);
                            Arrays.fill(readFractions, Float.NaN);
                        }
                    } else if (equal(currentLine, ii == 0 ? formatCellIndexes[ii] : formatCellIndexes[ii] + 1, formatCellIndexes[ii + 1], faBytes)) {
                        if (altAlleleFracRefHomThrehsold < 1 || altAlleleFractHetThrehsold > 0 || altAlleleFractAltHomThrehsold > 0) {
                            gtyAltAlleleFracIndexInInfor = ii;
                            hasIndexFA = true;
                            maxIndex2 = ii;
                            Arrays.fill(readFractions, Float.NaN);
                        }
                    } else if (equal(currentLine, ii == 0 ? formatCellIndexes[ii] : formatCellIndexes[ii] + 1, formatCellIndexes[ii + 1], plBytes)) {
                        if (minSecondPL > 0) {
                            gtyPLIndexInInfor = ii;
                            hasIndexPL = true;
                            maxIndex2 = ii;
                            Arrays.fill(secondMostGtyPL, 0);
                        }
                    } else if (equal(currentLine, ii == 0 ? formatCellIndexes[ii] : formatCellIndexes[ii] + 1, formatCellIndexes[ii + 1], gpBytes)) {
                        if (minBestGP > 0) {
                            gtyGPIndexInInfor = ii;
                            hasIndexGP = true;
                            maxIndex2 = ii;
                            Arrays.fill(bestGtyGP, 0);
                        }
                    }
                }

                //1/1:0,2:2:6.02:70,6,0	./.
                indexA = -1;
                for (int k = 0; k < effectiveIndivNum; k++) {
                    iGty = effectIndivIDInVCF.getQuick(k);
                    s = iGty + indexFORMAT + 1;
                    /*
                     //strange format "1/1:0,2:2:6.02:70,6,0"	"./."
                     if (currentLine.charAt(cellDelimts[s-1]+1) == '"') {
                     cellDelimts[s] = cellDelimts[s-1]+1; //does not work for the end index
                     }
                     */
//ASCII .
                    if (currentLine[cellDelimts[s] + 1] == 46) {
                        gtys[((iGty << 1))] = -1;
                        gtys[((iGty << 1)) + 1] = -1;
                        continue;
                    }
//                    System.out.println(k);
                    //':':58
                    tokenizeDelimiter(currentLine, cellDelimts[s] + 1, cellDelimts[s + 1], (byte) 58, formatCellIndexes);
                    if (hasNotCheckedPhase) {
                        //'|' 124
                        if (indexof(currentLine, gtyIndexInInfor == 0 ? formatCellIndexes[gtyIndexInInfor] : formatCellIndexes[gtyIndexInInfor] + 1, formatCellIndexes[gtyIndexInInfor + 1], (byte) 124) > 0
                                && indexof(currentLine, gtyIndexInInfor == 0 ? formatCellIndexes[gtyIndexInInfor] : formatCellIndexes[gtyIndexInInfor] + 1, formatCellIndexes[gtyIndexInInfor + 1], (byte) 47) < 0) {
                            isPhased = true;
                            hasNotCheckedPhase = false;
                        } else if (indexof(currentLine, gtyIndexInInfor == 0 ? formatCellIndexes[gtyIndexInInfor] : formatCellIndexes[gtyIndexInInfor] + 1, formatCellIndexes[gtyIndexInInfor + 1], (byte) 124) < 0
                                && indexof(currentLine, gtyIndexInInfor == 0 ? formatCellIndexes[gtyIndexInInfor] : formatCellIndexes[gtyIndexInInfor] + 1, formatCellIndexes[gtyIndexInInfor + 1], (byte) 47) > 0) {
                            isPhased = false;
                            hasNotCheckedPhase = false;
                        }
                    }
//Note this function is slower
//sByte = Arrays.copyOfRange(currentLine, cellDelimts[s] + 1, cellDelimts[s + 1]);
                    //tokenizeDelimiter(sByte, 0, sByte.length, (byte) ':', formatCellIndexes); 
                    if (gtyIndexInInfor >= 0) {
                        Arrays.fill(cellsBT, -1);
                        if (isPhased) {
                            //'|' 124
                            tokenizeDelimiter(currentLine, gtyIndexInInfor == 0 ? formatCellIndexes[gtyIndexInInfor] : formatCellIndexes[gtyIndexInInfor] + 1, formatCellIndexes[gtyIndexInInfor + 1], (byte) 124, cellsBT);
                            //tokenizeDelimiter(sByte, gtyIndexInInfor == 0 ? formatCellIndexes[gtyIndexInInfor] : formatCellIndexes[gtyIndexInInfor] + 1, formatCellIndexes[gtyIndexInInfor + 1], (byte) '|', cellsBT);
                        } else {
                            //'/' 47
                            tokenizeDelimiter(currentLine, gtyIndexInInfor == 0 ? formatCellIndexes[gtyIndexInInfor] : formatCellIndexes[gtyIndexInInfor] + 1, formatCellIndexes[gtyIndexInInfor + 1], (byte) 47, cellsBT);
                        }
                        gtys[(iGty << 1)] = parseInt(currentLine, cellsBT[0], cellsBT[1]);
                        //gtys[iGty][0] = parseInt(sByte, cellsBT[0], cellsBT[1]);
                        if (cellsBT[2] > 0) {
                            gtys[(iGty << 1) + 1] = parseInt(currentLine, cellsBT[1] + 1, cellsBT[2]);
                            // gtys[iGty][1] = parseInt(sByte, cellsBT[1] + 1, cellsBT[2]);
                        } else {
                            gtys[(iGty << 1) + 1] = gtys[(iGty << 1)];
                        }
                    }

                    if (gtyQualIndexInInfor >= 0) {
                        if (currentLine[gtyQualIndexInInfor == 0 ? formatCellIndexes[gtyQualIndexInInfor] : formatCellIndexes[gtyQualIndexInInfor] + 1] == 46) {
                            gtyQuality[iGty] = 0;
                        } else {
                            gtyQuality[iGty] = parseInt(currentLine, gtyQualIndexInInfor == 0 ? formatCellIndexes[gtyQualIndexInInfor] : formatCellIndexes[gtyQualIndexInInfor] + 1, formatCellIndexes[gtyQualIndexInInfor + 1]);
                        }

                    }
                    if (gtyDepthIndexInInfor >= 0) {
                        if (currentLine[gtyDepthIndexInInfor == 0 ? formatCellIndexes[gtyDepthIndexInInfor] : formatCellIndexes[gtyDepthIndexInInfor] + 1] == 46) {
                            gtyDepth[iGty] = 0;
                        } else {
                            gtyDepth[iGty] = parseInt(currentLine, gtyDepthIndexInInfor == 0 ? formatCellIndexes[gtyDepthIndexInInfor] : formatCellIndexes[gtyDepthIndexInInfor] + 1, formatCellIndexes[gtyDepthIndexInInfor + 1]);
                            if (indexA == -1 && readCounts[(iGty << 1)] > 0 && gtyAlleleDepthIndexInInfor != -1 && gtyDepthIndexInInfor > gtyAlleleDepthIndexInInfor) {
                                readCounts[(iGty << 1) + 1] = gtyDepth[iGty] - readCounts[(iGty << 1)];
                            }
                        }
                    }

                    // the format is a little bit compex. e.g., GT:AD:DP:GQ:PL	3/3:0,0,0,40:60:56:1200,1307,1629,1257,1440,1527,56,93,139,0
                    if (gtyAlleleDepthIndexInInfor >= 0) {
                        if (currentLine[gtyAlleleDepthIndexInInfor == 0 ? formatCellIndexes[gtyAlleleDepthIndexInInfor] : formatCellIndexes[gtyAlleleDepthIndexInInfor] + 1] != 46) {
                            Arrays.fill(cellsBT, -1);
                            len = tokenizeDelimiter(currentLine, gtyAlleleDepthIndexInInfor == 0 ? formatCellIndexes[gtyAlleleDepthIndexInInfor] : formatCellIndexes[gtyAlleleDepthIndexInInfor] + 1, formatCellIndexes[gtyAlleleDepthIndexInInfor + 1], (byte) 44, cellsBT);
                            if (len < 3) {
                                indexA = -1;
                                if (gtyDepthIndexInInfor < gtyAlleleDepthIndexInInfor) {
                                    readCounts[(iGty << 1) + 1] = parseInt(currentLine, cellsBT[0], cellsBT[1] + 1);
                                    readCounts[(iGty << 1)] = gtyDepth[iGty] - readCounts[(iGty << 1) + 1];
                                }
                            } else {
                                readCounts[(iGty << 1)] = parseInt(currentLine, cellsBT[0], cellsBT[1] + 1);
                                readCounts[(iGty << 1) + 1] = 0;
                                len = len - 1;
                                for (t = 1; t < len; t++) {
                                    //.
                                    if (currentLine[cellsBT[t] + 1] == 46) {
                                        continue;
                                    }
                                    //'0' : 48
                                    if (currentLine[cellsBT[t] + 1] == 48) {
                                        continue;
                                    }
                                    sc = parseInt(currentLine, cellsBT[t] + 1, cellsBT[t + 1]);
                                    readCounts[(iGty << 1) + 1] += sc;
                                }
                            }
                        }
                    }
                    if (gtyAltAlleleFracIndexInInfor >= 0) {
                        if (currentLine[gtyAltAlleleFracIndexInInfor == 0 ? formatCellIndexes[gtyAltAlleleFracIndexInInfor] : formatCellIndexes[gtyAltAlleleFracIndexInInfor] + 1] == 46) {
                            readFractions[iGty] = Float.NaN;
                        } else {
                            readFractions[iGty] = parseFloat(currentLine, gtyAltAlleleFracIndexInInfor == 0 ? formatCellIndexes[gtyAltAlleleFracIndexInInfor] : formatCellIndexes[gtyAltAlleleFracIndexInInfor] + 1, formatCellIndexes[gtyAltAlleleFracIndexInInfor + 1]);
                        }
                    }

                    if (gtyPLIndexInInfor >= 0) {
                        if (currentLine[gtyPLIndexInInfor == 0 ? formatCellIndexes[gtyPLIndexInInfor] : formatCellIndexes[gtyPLIndexInInfor] + 1] != 46) {
                            len = tokenizeDelimiter(currentLine, gtyPLIndexInInfor == 0 ? formatCellIndexes[gtyPLIndexInInfor] : formatCellIndexes[gtyPLIndexInInfor] + 1, formatCellIndexes[gtyPLIndexInInfor + 1], (byte) 44, cellsPL);

                            if (currentLine[gtyPLIndexInInfor == 0 ? formatCellIndexes[gtyPLIndexInInfor] : formatCellIndexes[gtyPLIndexInInfor] + 1] != 46) {
                                secondMostGtyPL[iGty] = parseInt(currentLine, cellsPL[0], cellsPL[1]);
                                if (secondMostGtyPL[iGty] == 0) {
                                    secondMostGtyPL[iGty] = Integer.MAX_VALUE;
                                }
                            } else {
                                secondMostGtyPL[iGty] = Integer.MAX_VALUE;
                            }
                            len = len - 1;
                            for (t = 1; t < len; t++) {
                                //.
                                if (currentLine[cellsPL[t] + 1] == 46) {
                                    continue;
                                }
                                //'0' : 48
                                if (currentLine[cellsPL[t] + 1] == 48) {
                                    continue;
                                }
                                sc = parseInt(currentLine, cellsPL[t] + 1, cellsPL[t + 1]);

                                if (sc < secondMostGtyPL[iGty]) {
                                    secondMostGtyPL[iGty] = sc;
                                }
                            }
                        }
                    }
                    if (gtyGPIndexInInfor >= 0) {
                        if (currentLine[gtyGPIndexInInfor == 0 ? formatCellIndexes[gtyGPIndexInInfor] : formatCellIndexes[gtyGPIndexInInfor] + 1] != 46) {
                            len = tokenizeDelimiter(currentLine, gtyGPIndexInInfor == 0 ? formatCellIndexes[gtyGPIndexInInfor] : formatCellIndexes[gtyGPIndexInInfor] + 1, formatCellIndexes[gtyGPIndexInInfor], (byte) 44, cellsPL);

                            if (currentLine[formatCellIndexes[gtyGPIndexInInfor - 1] + 1] != 46) {
                                bestGtyGP[iGty] = parseInt(currentLine, cellsPL[0], cellsPL[1]);
                                if (bestGtyGP[iGty] == 0) {
                                    bestGtyGP[iGty] = Integer.MAX_VALUE;
                                }
                            } else {
                                bestGtyGP[iGty] = Integer.MAX_VALUE;
                            }
                            len = len - 1;
                            for (t = 1; t < len; t++) {
                                if (currentLine[cellsPL[t] + 1] == 46) {
                                    continue;
                                }
                                sc = parseInt(currentLine, cellsPL[t] + 1, cellsPL[t + 1]);
                                if (sc == 0) {
                                    continue;
                                }
                                if (sc < bestGtyGP[iGty]) {
                                    bestGtyGP[iGty] = sc;
                                }
                            }
                        }
                    }
                }

                /*
                 if (!isLowQualBreak && currChr.indexOf(UNKNOWN_CHROM_NAME0) < 0 && currChr.indexOf(UNKNOWN_CHROM_NAME1) < 0 && hasAlt) {
                 formatProbVarNum++;
                 LOG.error("Format error at line " + fileLineCounter + ": " + currentLine);
                 continue;
                 }
                 * 
                 */
                //QC
                obsS = genotypeQC();
                /*
                 if (obsS == 0) {
                 ignoredLineNumNoVar++;
                 continue;
                 }
                 * 
                 */
                if (obsS < minOBS) {
                    ignoredLineNumMinOBS++;
                    continue;
                }

                if ((!needAccoundAffect && !needAccoundUnaffect) || needMAFQCOver || needMAFQCLess) {
                    g11 = 0;
                    g12 = 0;
                    g22 = 0;

                    for (int i = 0; i < totalPedSubjectNum; i++) {
                        int idLabel = pedVCFIDMap[i];
                        if (idLabel < 0) {
                            continue;
                        }
                        if (gtys[(idLabel << 1)] == -1) {
                            continue;
                        }
                        if (gtys[(idLabel << 1) + 1] == -1) {
                            if (gtys[(idLabel << 1)] == 0) {
                                g11++;
                            } else {
                                g22++;
                            }
                        } else {
                            if (gtys[(idLabel << 1)] != gtys[(idLabel << 1) + 1]) {
                                g12++;
                            } else if (gtys[(idLabel << 1)] == 0) {
                                g11++;
                            } else {
                                g22++;
                            }
                        }
                    }
                }
                if (needMAFQCOver || needMAFQCLess) {
                    maf = (g12 * 0.5 + g22) / (g11 + g12 + g22);
                    if (needMAFQCOver) {
                        if (maf <= sampleMafOver || maf >= sampleMafOverC) {
                            ignoredLineNumMinMAF++;
                            continue;
                        }
                    }
                    if (needMAFQCLess) {
                        if (maf >= sampleMafLess && maf <= sampleMafLessC) {
                            ignoredLineNumMaxMAF++;
                            continue;
                        }
                    }
                }

                Variant var = new Variant(makerPostion, ref, altAllelesIndexes);
                var.setIsIndel(isIndel);
                var.setLabel(varLabel);
                varChroms[chromID].add(var);

                if (needGtyQual) {
                    System.arraycopy(currentLine, 0, tmpByte, 0, cellDelimts[9]);
                    lastTmpByteIndex = cellDelimts[9];
                }

                if (needGty) {
                    encodeGenotype(effectiveIndivNum, var);
                }

                if (needReadsInfor) {
                    var.readInfor = new char[effectiveIndivNum * 2];
                    //two chars for a gty 
                    for (index = 0; index < totalPedSubjectNum; index++) {
                        vcfID = pedVCFIDMap[index];
                        if (vcfID < 0) {
                            continue;
                        }

                        var.readInfor[(pedEncodeGytID[index] << 1)] = (char) readCounts[(vcfID << 1)];
                        var.readInfor[(pedEncodeGytID[index] << 1) + 1] = (char) readCounts[(vcfID << 1) + 1];
                    }
                }

                if (needGtyQual) {
                    for (index = 0; index < totalPedSubjectNum; index++) {
                        int idLabel = pedVCFIDMap[index];
                        //if the subject is not in the original vcf, just ignore it
                        if (idLabel < 0) {
                            /*
                             tmpByte[lastTmpByteIndex] = 46;
                             lastTmpByteIndex++;
                             tmpByte[lastTmpByteIndex] = 47;
                             lastTmpByteIndex++;
                             tmpByte[lastTmpByteIndex] = 46;
                             lastTmpByteIndex++;
                             */
                            continue;
                        }

                        //\t:9
                        tmpByte[lastTmpByteIndex] = 9;
                        lastTmpByteIndex++;

                        if (gtys[(idLabel << 1)] == -1) {
                            tmpByte[lastTmpByteIndex] = 46;
                            lastTmpByteIndex++;
                            /*
                             tmpByte[lastTmpByteIndex] = 47;
                             lastTmpByteIndex++;
                             tmpByte[lastTmpByteIndex] = 46;
                             lastTmpByteIndex++;
                             */
                            continue;
                        }
                        s = idLabel + indexFORMAT + 1;
                        t = cellDelimts[s + 1] - cellDelimts[s] - 1;
                        System.arraycopy(currentLine, cellDelimts[s] + 1, tmpByte, lastTmpByteIndex, t);
                        lastTmpByteIndex += t;
                    }

                    tmpByte[lastTmpByteIndex] = 10;
                    lastTmpByteIndex++;
                    String newString = new String(Arrays.copyOfRange(tmpByte, 0, lastTmpByteIndex));
                    vcfWriters[chromID].write(newString);
                }

                if (needAccoundAffect) {
                    g11 = 0;
                    g12 = 0;
                    g22 = 0;

                    for (int i = 0; i < caseSize; i++) {
                        int idLabel = pedVCFIDMap[caeSetID[i]];
                        if (idLabel < 0) {
                            continue;
                        }
                        if (gtys[(idLabel << 1)] == -1) {
                            continue;
                        }

                        if (gtys[(idLabel << 1)] != gtys[(idLabel << 1) + 1]) {
                            g12++;
                        } else if (gtys[(idLabel << 1)] == 0) {
                            g11++;
                        } else {
                            g22++;
                        }
                    }
                    var.setAffectedRefHomGtyNum(g11);
                    var.setAffectedHetGtyNum(g12);
                    var.setAffectedAltHomGtyNum(g22);
                }

                if (needAccoundUnaffect) {
                    g11 = 0;
                    g12 = 0;
                    g22 = 0;

                    for (int i = 0; i < controlSize; i++) {
                        int idLabel = pedVCFIDMap[controlSetID[i]];
                        if (idLabel < 0) {
                            continue;
                        }
                        if (gtys[(idLabel << 1)] == -1) {
                            continue;
                        }

                        if (gtys[(idLabel << 1)] != gtys[(idLabel << 1) + 1]) {
                            g12++;
                        } else if (gtys[(idLabel << 1)] == 0) {
                            g11++;
                        } else {
                            g22++;
                        }
                    }
                    var.setUnaffectedRefHomGtyNum(g11);
                    var.setUnaffectedHetGtyNum(g12);
                    var.setUnaffectedAltHomGtyNum(g22);
                }

                if (!needAccoundAffect && !needAccoundUnaffect) {
                    var.setAffectedRefHomGtyNum(g11);
                    var.setAffectedHetGtyNum(g12);
                    var.setAffectedAltHomGtyNum(g22);
                }
                acceptVarNum++;
                if (acceptVarNum >= maxVarNum) {
                    // System.out.println(acceptVarNum);
                    totalAcceptVarNum += acceptVarNum;
                    writeChromosomeToDiskClean();
                    acceptVarNum = 0;
                }
        
            }

            if (acceptVarNum > 0) {
                writeChromosomeToDiskClean();
                totalAcceptVarNum += acceptVarNum;
                acceptVarNum = 0;
            }
        } catch (Exception nex) {
            nex.printStackTrace();
            String inforS = new String(currentLine);
            String info = nex.toString() + " when parsing at line : " + inforS;
            LOG.error(info);
            isInvalid = true;
            // throw new Exception(info);
            // LOG.error(nex, info);
        }

        //change the order to be consisten with the pedigree file
        // effectIndivIDInVCF.quickSort();
        return acceptVarNum;
    }

    private int genotypeQC() {
        int obsS = 0;
        for (int k = totalPedSubjectNum - 1; k >= 0; k--) {
            index = pedVCFIDMap[k];
            if (index < 0) {
                continue;
            }

            //ignore variants with missing genotypes
            if (gtys[(index << 1)] == -1) {
                missingGtyNum++;
                continue;
            }

            if (hasIndexGQ && gtyQuality[index] < gtyQualityThrehsold) {
                gtys[(index << 1)] = -1;
                gtys[(index << 1) + 1] = -1;
                ignoredLowQualGtyNum++;
                continue;
            }
            if (hasIndexDP && gtyDepth[index] < minGtySeqDepth) {
                gtys[(index << 1)] = -1;
                gtys[(index << 1) + 1] = -1;
                ignoredLowDepthGtyNum++;
                continue;
            }
            if (hasIndexPL && secondMostGtyPL[index] < minSecondPL) {
                ignoredLowPLGtyNum++;
                gtys[(index << 1)] = -1;
                gtys[(index << 1) + 1] = -1;
                continue;
            }
            if (hasIndexGP && bestGtyGP[index] < minBestGP) {
                ignoredLowGPGtyNum++;
                gtys[(index << 1)] = -1;
                gtys[(index << 1) + 1] = -1;
                continue;
            }

            if (hasIndexAD || hasIndexFA) {
                if (gtys[(index << 1)] == 0 && gtys[(index << 1) + 1] == 0) {
                    if (altAlleleFracRefHomThrehsold < 1) {
                        //the AD infor may be missing
                        if (readCounts[(index << 1)] == -1 && Float.isNaN(readFractions[index])) {
                            continue;
                        }
                        if (Float.isNaN(readFractions[index])) {
                            readFractions[index] = (float) (((double) readCounts[(index << 1) + 1]) / (readCounts[(index << 1)] + readCounts[(index << 1) + 1]));
                        }

                        if (Float.isNaN(readFractions[index])) {
                            gtys[(index << 1)] = -1;
                            gtys[(index << 1) + 1] = -1;
                            ignoredBadAltFracGtyNum++;
                            continue;
                        } else if (readFractions[index] > altAlleleFracRefHomThrehsold) {
                            gtys[(index << 1)] = -1;
                            gtys[(index << 1) + 1] = -1;
                            ignoredBadAltFracGtyNum++;
                            continue;
                        }

                    }
                } else if (gtys[(index << 1)] != gtys[(index << 1) + 1]) {
                    if (altAlleleFractHetThrehsold > 0) {
                        //the AD infor may be missing
                        if (readCounts[(index << 1)] == -1 && Float.isNaN(readFractions[index])) {
                            continue;
                        }
                        if (Float.isNaN(readFractions[index])) {
                            readFractions[index] = (float) (((double) readCounts[(index << 1) + 1]) / (readCounts[(index << 1)] + readCounts[(index << 1) + 1]));
                        }
                        if (Float.isNaN(readFractions[index])) {
                            gtys[(index << 1)] = -1;
                            gtys[(index << 1) + 1] = -1;
                            ignoredBadAltFracGtyNum++;
                            continue;
                        } else if (readFractions[index] < altAlleleFractHetThrehsold) {
                            gtys[(index << 1)] = -1;
                            gtys[(index << 1) + 1] = -1;
                            ignoredBadAltFracGtyNum++;
                            continue;
                        }
                    }

                } else {
                    if (altAlleleFractAltHomThrehsold > 0) {
                        //the AD infor may be missing
                        if (readCounts[(index << 1) + 1] == -1 && Float.isNaN(readFractions[index])) {
                            continue;
                        }
                        if (Float.isNaN(readFractions[index])) {
                            readFractions[index] = (float) (((double) readCounts[(index << 1) + 1]) / (readCounts[(index << 1)] + readCounts[(index << 1) + 1]));
                        }
                        if (Float.isNaN(readFractions[index])) {
                            gtys[(index << 1)] = -1;
                            gtys[(index << 1) + 1] = -1;
                            ignoredBadAltFracGtyNum++;
                            continue;
                        } else if (readFractions[index] < altAlleleFractAltHomThrehsold) {
                            gtys[(index << 1)] = -1;
                            gtys[(index << 1) + 1] = -1;
                            ignoredBadAltFracGtyNum++;
                            continue;
                        }
                    }
                }
            }
            obsS++;
        }
        return obsS;
    }

    private void encodeGenotype(int effSubjectNum, Variant var) {
        int base, bitNum, byteNum, alleleNumShift, byteIndex1, byteIndex2, s1;
        //as the genotypes may be use for other purpose so we need record it before filtering 
        if (isPhased) {
            base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
            bitNum = base * effSubjectNum;

            if (bitNum % 32 != 0) {
                byteNum = bitNum / 32 + 1;
            } else {
                byteNum = bitNum / 32;
            }

            var.encodedGty = new int[byteNum];
            Arrays.fill(var.encodedGty, 0);
            alleleNumShift = alleleNum << 16;
            for (index = 0; index < totalPedSubjectNum; index++) {
                int vcfID = pedVCFIDMap[index];
                if (vcfID < 0) {
                    continue;
                }

                bitNum = base * pedEncodeGytID[index];

                if (gtys[(vcfID << 1)] == -1) {
                    //missing value         
                    continue;
                } else {
                    switch (alleleNum) {
                        case 2:
                            /*       
                             missing	Reference homozygous	Heterozygous 	Heterozygous 	Alternative homozygous
                             VCF genotype	.|.	0|0	0|1	1|0	1|1
                             Bits	        000  	001	010	011	100
                             Order	0	1	2	3	4                
               
                             II.II Tri-allelic sequence variant (4 bits)
                             missing 	Reference homozygous 	Heterozygous 	Heterozygous 	Heterozygous 	Heterozygous 	Alternative homozygous
                             VCF genotype 	.|. 	0|0 	0|1 	0|2 	1|0 	1|1 	1|2
                             Bits      	000 	0001 	0010 	0011 	0100 	0101 	0110
                             Decimal 	0 	1 	2 	3 	4 	5 	6
                             Heterozygous 	Heterozygous 	Alternative homozygous
                             VCF genotype 	2|0 	2|1 	2|2
                             Bits     	0111 	1000 	1001
                             Decimal 	7 	8 	9     
                             */

                            //to speedup the analysis
                            if (gtys[(vcfID << 1)] == 0 && gtys[(vcfID << 1) + 1] == 0) {
                                byteIndex1 = (bitNum + 2) / 32;
                                byteIndex2 = (bitNum + 2) % 32;
                                // System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                                var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                                //System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                            } else if (gtys[(vcfID << 1)] == 0 && gtys[(vcfID << 1) + 1] == 1) {
                                byteIndex1 = (bitNum + 1) / 32;
                                byteIndex2 = (bitNum + 1) % 32;
                                // System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                                var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                                // System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                            } else if (gtys[(vcfID << 1)] == 1 && gtys[(vcfID << 1) + 1] == 0) {
                                byteIndex1 = (bitNum + 1) / 32;
                                byteIndex2 = (bitNum + 1) % 32;
                                // System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                                var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                                byteIndex1 = (bitNum + 2) / 32;
                                byteIndex2 = (bitNum + 2) % 32;
                                // System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                                var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                            } else if (gtys[(vcfID << 1)] == 1 && gtys[(vcfID << 1) + 1] == 1) {
                                byteIndex1 = (bitNum) / 32;
                                byteIndex2 = (bitNum) % 32;
                                var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                            }
                            break;
                        default:
                            boolean[] bits = GlobalManager.phasedGtyCodingMap.get(gtys[(vcfID << 1)] | (gtys[(vcfID << 1) + 1] << 8) | (alleleNumShift));
                            for (int i = 0; i < base; i++) {
                                if (bits[i]) {
                                    byteIndex1 = (bitNum + i) / 32;
                                    byteIndex2 = (bitNum + i) % 32;
                                    var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                                }
                            }
                    }
                }
            }
        } else {
            base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
            bitNum = base * effSubjectNum;

            if (bitNum % 32 != 0) {
                byteNum = bitNum / 32 + 1;
            } else {
                byteNum = bitNum / 32;
            }

            var.encodedGty = new int[byteNum];
            Arrays.fill(var.encodedGty, 0);
            alleleNumShift = alleleNum << 16;
            for (index = 0; index < totalPedSubjectNum; index++) {
                int vcfID = pedVCFIDMap[index];
                if (vcfID < 0) {
                    continue;
                }

                /*
                 missing	Reference homozygous	Heterozygous 	Alternative homozygous
                 VCF genotype	./.	0/0	0/1	1/1
                 Bits	00  	01	10	11
                 Order	0	1	2	3        
                 */

                /*
                 missing	Reference homozygous	Heterozygous 	Heterozygous	Alternative homozygous	Heterozygous	Alternative homozygous
                 VCF genotype	./.	0/0	0/1	0/2	1/1	1/2	2/2
                 Bits	        000	001	010	011	100	101	110
                 Order	0	1	2	3	4	5	6 
                 */
                /*
                 I.III Quad-allelic sequence variant (4 bits)
                 missing 	Reference homozygous 	Heterozygous 	Heterozygous 	Heterozygous 	Alternative homozygous 	Heterozygous
                 VCF genotype 	./. 	0/0 	0/1 	0/2 	0/3 	1/1 	1/2
                 Bits 	      000 	0001 	0010 	0011 	0100 	0101 	0110
                 Decimal 	0 	1 	2 	3 	4 	5 	6
                 Heterozygous 	Alternative homozygous 	Heterozygous 	Alternative homozygous
                 VCF genotype 	1/3 	2/2 	2/3 	3/3
                 Bits 	     0111 	1000 	1001 	1010
                 Decimal 	7 	8 	9 	10                               
                 */
                bitNum = base * pedEncodeGytID[index];
                if (gtys[(vcfID << 1)] == -1) {
                    //missing value         
                    continue;
                } else {
                    switch (alleleNum) {
                        case 2:
                            /*
                             missing	Reference homozygous	Heterozygous 	Alternative homozygous
                             VCF genotype	./.	0/0	0/1	1/1
                             Bits	00  	01	10	11
                             Order	0	1	2	3        
                             */

                            //to speedup the analysis
                            if (gtys[(vcfID << 1)] == 0 && gtys[(vcfID << 1) + 1] == 0) {
                                byteIndex1 = (bitNum + 1) / 32;
                                byteIndex2 = (bitNum + 1) % 32;
                                // System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                                var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                                //System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                            } else if (gtys[(vcfID << 1)] == 0 && gtys[(vcfID << 1) + 1] == 1 || gtys[(vcfID << 1)] == 1 && gtys[(vcfID << 1) + 1] == 0) {
                                byteIndex1 = (bitNum) / 32;
                                byteIndex2 = (bitNum) % 32;
                                // System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                                var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                                // System.out.println(Integer.toBinaryString(var.encodedGty[byteIndex1]));
                            } else if (gtys[(vcfID << 1)] == 1 && gtys[(vcfID << 1) + 1] == 1) {
                                byteIndex1 = (bitNum) / 32;
                                byteIndex2 = (bitNum) % 32;
                                var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                                byteIndex1 = (bitNum + 1) / 32;
                                byteIndex2 = (bitNum + 1) % 32;
                                var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                            }
                            break;
                        default:
                            if (gtys[(vcfID << 1)] > gtys[(vcfID << 1) + 1]) {
                                s1 = gtys[(vcfID << 1)];
                                gtys[(vcfID << 1)] = gtys[(vcfID << 1) + 1];
                                gtys[(vcfID << 1) + 1] = s1;
                            }
                            boolean[] bits = GlobalManager.unphasedGtyCodingMap.get(gtys[(vcfID << 1)] | (gtys[(vcfID << 1) + 1] << 8) | (alleleNumShift));

                            for (int i = 0; i < base; i++) {
                                if (bits[i]) {
                                    byteIndex1 = (bitNum + i) / 32;
                                    byteIndex2 = (bitNum + i) % 32;
                                    var.encodedGty[byteIndex1] = var.encodedGty[byteIndex1] | GlobalManager.intOpers[byteIndex2];
                                }
                            }
                    }
                }
            }
        }
    }

    private void writeChromosomeToDiskClean() {
        int chromID = -1;
        String chromeName;
        try {
            for (List<Variant> varList : varChroms) {
                chromID++;
                chromeName = STAND_CHROM_NAMES[chromID];
                if (varList.isEmpty()) {
                    continue;
                }
                int fileIndex = -1;
                String chrNameP = "Chromosome." + chromeName;
                chrNameP = threadID + "." + chrNameP;
                File fileName = null;
                File folder = new File(storagePath);
                if (folder.exists()) {
                    do {
                        fileIndex++;
                        fileName = new File(storagePath + File.separator + chrNameP + ".var.obj." + fileIndex);
                    } while (fileName.exists());
                } else {
                    fileIndex++;
                    folder.mkdirs();
                }

                //comments: both Kryo and FSTObjectOutput are excellent tools for Serializationl. However, the former produced slightly smaller file and was slightly faster. So I used Kryo
                fileName = new File(storagePath + File.separator + chrNameP + ".var.obj." + fileIndex);
                //slower
                //  kryo.setInstantiatorStrategy(new SerializingInstantiatorStrategy());            
                Output output = new Output(new FileOutputStream(fileName), 1024 * 1024);
                //  output.setBuffer(buffer);

                for (Variant var : varList) {
                    kryo.writeObject(output, var);
                }
                output.flush();
                output.close();
                varList.clear();
            }
            System.gc();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void setBooleanFilter(boolean considerSNP, boolean considerIndel, boolean needGty, boolean needReadsInfor, boolean needGtyQual, boolean noGtyVCF) {
        this.considerSNP = considerSNP;
        this.considerIndel = considerIndel;
        this.needGty = needGty;
        this.needReadsInfor = needReadsInfor;
        this.needGtyQual = needGtyQual;
        this.noGtyVCF = noGtyVCF;
    }

    private IntArrayList effectIndivIDInVCF;
    private int[] pedVCFIDMap;
    private int[] pedEncodeGytID;
    private List<Individual> subjectList;
    private boolean isPhased = false;
    Genome orgGenome;

    public boolean isIsPhased() {
        return isPhased;
    }

    public List<Variant>[] getVarChroms() {
        return varChroms;
    }

    public List<Individual> getSubjectList() {
        return subjectList;
    }

    public void setOrgGenome(Genome orgGenome) {
        this.orgGenome = orgGenome;
    }

    public void setGenotypesAndSubjects(IntArrayList effectIndivID, List<Individual> subjectList1, int[] pedVCFIDMap, int[] pedEncodeGytIDMap, int[] caeSetID, int[] controlSetID) {
        this.effectIndivIDInVCF = effectIndivID;
        this.subjectList = new ArrayList<Individual>();
        if (subjectList1 != null && !subjectList1.isEmpty()) {
            subjectList.addAll(subjectList1);
        }
        this.caeSetID = caeSetID;
        this.controlSetID = controlSetID;
        this.pedVCFIDMap = pedVCFIDMap;
        this.pedEncodeGytID = pedEncodeGytIDMap;
    }

    public VCFParseTaskFast1(int threadID) {
        this.threadID = threadID;
        varChroms = new ArrayList[STAND_CHROM_NAMES.length];
        for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
            varChroms[i] = new ArrayList<>();
        }
        for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
            chromNameIndexMap.put(STAND_CHROM_NAMES[i], i);
        }
        kryo.setReferences(false);
        kryo.setRegistrationRequired(false);
        kryo.setInstantiatorStrategy(new StdInstantiatorStrategy());
        kryo.register(char[].class);
        kryo.register(long[].class);
        kryo.register(float[].class);
        kryo.register(String[].class);

        kryo.register(StringBuilder.class);
        kryo.register(Variant.class);
        kryo.register(OpenIntIntHashMap.class);
    }

    public void setVcfLabelSet(Set<String> vcfLabelSet) {
        this.vcfLabelSet = vcfLabelSet;
    }

    public void setQuantitativeQCParams(double avgSeqQualityThrehsold, double minMappingQual, double maxStrandBias, double maxFisherStrandBias, int maxGtyAlleleNum, double gtyQualityThrehsold, int minGtySeqDepth, int minSeqDepth, double altAlleleFracRefHomThrehsold, double altAlleleFractHetThrehsold, double altAlleleFractAltHomThrehsold, int minSecondPL, double minBestGP, int minOBS, double sampleMafOver, double sampleMafLess) {
        this.avgSeqQualityThrehsold = avgSeqQualityThrehsold;
        this.minMappingQual = minMappingQual;
        this.maxStrandBias = maxStrandBias;
        this.maxFisherStrandBias = maxFisherStrandBias;
        this.gtyQualityThrehsold = gtyQualityThrehsold;
        this.minGtySeqDepth = minGtySeqDepth;
        this.minSeqDepth = minSeqDepth;
        this.altAlleleFracRefHomThrehsold = altAlleleFracRefHomThrehsold;
        this.altAlleleFractHetThrehsold = altAlleleFractHetThrehsold;
        this.altAlleleFractAltHomThrehsold = altAlleleFractAltHomThrehsold;
        this.minSecondPL = minSecondPL;
        this.minBestGP = minBestGP;
        this.minOBS = minOBS;
        this.sampleMafOver = sampleMafOver;
        this.sampleMafLess = sampleMafLess;

        this.maxGtyAlleleNum = maxGtyAlleleNum;
    }

    public void setColIndex(int indexCHROM, int indexPOS, int indexID, int indexREF, int indexALT, int indexQUAL, int indexFILTER, int indexINFO, int indexFORMAT) {
        this.indexCHROM = indexCHROM;
        this.indexPOS = indexPOS;
        this.indexID = indexID;
        this.indexREF = indexREF;
        this.indexALT = indexALT;
        this.indexQUAL = indexQUAL;
        this.indexFILTER = indexFILTER;
        this.indexFORMAT = indexFORMAT;
        this.indexINFO = indexINFO;
    }

    @Override
    public String call() throws Exception {
        long startTime = System.currentTimeMillis();

        try {
            if (needGtyQual) {
                vcfWriters = new BufferedWriter[STAND_CHROM_NAMES.length];
                for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
                    String chrNameP = "Chromosome." + STAND_CHROM_NAMES[i];
                    chrNameP = threadID + "." + chrNameP;
                    int fileIndex = -1;

                    File fileName = null;
                    File folder = new File(storagePath);
                    if (folder.exists()) {
                        do {
                            fileIndex++;
                            fileName = new File(storagePath + File.separator + chrNameP + ".vcf.gz." + fileIndex);
                            //the initial size of the compressed file is 
                            if (fileName.length() <= 20) {
                                break;
                            }
                        } while (fileName.exists());
                    } else {
                        fileIndex++;
                        folder.mkdirs();
                    }

                    //comments: both Kryo and FSTObjectOutput are excellent tools for Serializationl. However, the former produced slightly smaller file and was slightly faster. So I used Kryo
                    fileName = new File(storagePath + File.separator + chrNameP + ".vcf.gz." + fileIndex);

                    GZIPOutputStream gzOut = new GZIPOutputStream(new FileOutputStream(fileName));
                    vcfWriters[i] = new BufferedWriter(new OutputStreamWriter(gzOut));
                }
            }
        } catch (Exception ex) {
            LOG.error(ex);
        }
        parseVariantsInFileOnlyFastToken();

        fireTaskComplete();
        try {
            if (needGtyQual) {
                for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
                    vcfWriters[i].close();
                }
            }
        } catch (Exception ex) {
            LOG.error(ex);
        }
        String info = "Elapsed time: " + (System.currentTimeMillis() - startTime) / 1000 + " seconds.";
        //  System.out.println(info);
        //return info;
        return info;
        // 
    }

}
