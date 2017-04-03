/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;
import org.apache.log4j.Logger;
import org.cobi.kggseq.Constants;
import org.cobi.kggseq.GlobalManager;
import org.cobi.kggseq.entity.Chromosome;
import org.cobi.kggseq.entity.Genome;
import org.cobi.kggseq.entity.Individual;
import org.cobi.kggseq.entity.Variant;
import org.cobi.util.file.LocalFileFunc;

/**
 *
 * @author mxli
 */
public class BinaryGtyProcessor implements Constants {

    private static final Logger LOG = Logger.getLogger(BinaryGtyProcessor.class);
    protected String pedigreeFileName;
    protected String mapFileName;
    //if kggseqBinaryFileName is not null, it must be a kggseq binary file
    protected String kggseqBinaryFileName;

//tmp variants
    static int byteIndex1;
    static int byteIndex2;
    static int bitNum;
    static boolean[] bits = new boolean[32];
    static StringBuilder stringBuilder = new StringBuilder();
    static boolean needGzExtension = false;

    public BinaryGtyProcessor(String prefixName) {
        this.pedigreeFileName = prefixName + ".kam";
        this.mapFileName = prefixName + ".kim";
        this.kggseqBinaryFileName = prefixName + ".ked";
    }

    public static final int[] getUnphasedGtyAt(byte[] gtys, int alleleNum, int base, int indivID, int blockSize) {
        bitNum = indivID;
        for (int i = 0; i < base; i++) {
            byteIndex1 = (bitNum) / 8;
            byteIndex2 = (bitNum) % 8;
            bits[i] = (gtys[byteIndex1] & GlobalManager.byteOpers[byteIndex2]) == GlobalManager.byteOpers[byteIndex2];
            bitNum += blockSize;
        }
        switch (base) {
            case 2:
                /*
                 missing	  Reference homozygous	 Heterozygous 	Alternative homozygous
                 VCF genotype	  0/0	       0/1	         1/1  ./.	           
                 Bits           00          01	       10	         11
                 Order	         0	              1	                2	         3        
                 */
                if (!bits[0] && !bits[1]) {
                    return new int[]{0, 0};
                } else if (!bits[0] && bits[1]) {
                    return new int[]{0, 1};
                } else if (bits[0] && !bits[1]) {
                    return new int[]{1, 1};
                } else if (bits[0] && bits[1]) {
                    return null;
                }
                break;
            default:
                stringBuilder.delete(0, stringBuilder.length());
                for (int i = 0; i < base; i++) {
                    if (bits[i]) {
                        stringBuilder.append(1);
                    } else {
                        stringBuilder.append(0);
                    }
                }
                stringBuilder.append(':').append(alleleNum);
                int[] alleles = GlobalManager.codingUnphasedGtyCodingMap.get(stringBuilder.toString());
                // String infor = "Sorry!!! squence variants with over 4 alleles are not supported and will be ignored!";
                // System.out.println(infor);
                return alleles;
        }
        return null;
    }

    public static final int[] getUnphasedGtyBool(boolean[] gtyBits, int alleleNum, int base, int indivID) {
        switch (base) {
            case 2:
                /*
                 missing	  Reference homozygous	 Heterozygous 	Alternative homozygous
                 VCF genotype	  0/0	       0/1	         1/1  ./.	           
                 Bits           00          01	       10	         11
                 Order	         0	              1	                2	         3        
                 */
                if (!gtyBits[0] && !gtyBits[1]) {
                    return new int[]{0, 0};
                } else if (!gtyBits[0] && gtyBits[1]) {
                    return new int[]{0, 1};
                } else if (gtyBits[0] && !gtyBits[1]) {
                    return new int[]{1, 1};
                } else if (gtyBits[0] && gtyBits[1]) {
                    return null;
                }
                break;
            default:
                stringBuilder.delete(0, stringBuilder.length());
                for (int i = 0; i < base; i++) {
                    if (gtyBits[i]) {
                        stringBuilder.append(1);
                    } else {
                        stringBuilder.append(0);
                    }
                }
                stringBuilder.append(':').append(alleleNum);
                int[] alleles = GlobalManager.codingUnphasedGtyCodingMap.get(stringBuilder.toString());
                // String infor = "Sorry!!! squence variants with over 4 alleles are not supported and will be ignored!";
                // System.out.println(infor);
                return alleles;
        }
        return null;
    }

    public static final void getPhasedGtyAt1(byte[] gtys, int base, int indivID, boolean[] bits1, int blockSize) {
        bitNum = indivID;
        for (int i = 0; i < base; i++) {
            byteIndex1 = (bitNum) / 8;
            byteIndex2 = (bitNum) % 8;
            bits1[i] = (gtys[byteIndex1] & GlobalManager.byteOpers[byteIndex2]) == GlobalManager.byteOpers[byteIndex2];
            bitNum += blockSize;
        }
    }

    public static final int[] getPhasedGtyAt(byte[] gtys, int alleleNum, int base, int indivID, int blockSize) {
        bitNum = indivID;
        for (int i = 0; i < base; i++) {
            byteIndex1 = (bitNum) / 8;
            byteIndex2 = (bitNum) % 8;
            bits[i] = (gtys[byteIndex1] & GlobalManager.byteOpers[byteIndex2]) == GlobalManager.byteOpers[byteIndex2];
            bitNum += blockSize;
        }

        switch (base) {
            case 2:
                /*       
                 missing	Reference homozygous	Heterozygous 	Heterozygous 	Alternative homozygous
                 VCF genotype		0|0	0|1	1|0	1|1 .|.
                 Bits	        000  	001	010	011	100
                 Order	0	1	2	3	4                
               
                 II.II Tri-allelic sequence variant (4 bits)
                 missing 	Reference homozygous 	Heterozygous 	Heterozygous 	Heterozygous 	Heterozygous 	Alternative homozygous
                 VCF genotype 	0|0 	0|1 	0|2 	1|0 	1|1 	1|2
                 Bits      	000 	0001 	0010 	0011 	0100 	0101 	0110
                 Decimal 	0 	1 	2 	3 	4 	5 	6
                 Heterozygous 	Heterozygous 	Alternative homozygous
                 VCF genotype 	2|0 	2|1 	2|2 	.|. 
                 Bits     	0111 	1000 	1001
                 Decimal 	7 	8 	9     
                 */
                if (!bits[0] && !bits[1] && !bits[2]) {
                    return new int[]{0, 0};
                } else if (!bits[0] && !bits[1] && bits[2]) {
                    return new int[]{0, 1};
                } else if (!bits[0] && bits[1] && !bits[2]) {
                    return new int[]{1, 0};
                } else if (!bits[0] && bits[1] && bits[2]) {
                    return new int[]{1, 1};
                } else if (bits[0] && !bits[1] && !bits[2]) {
                    return null;
                }
                break;
            default:
                stringBuilder.delete(0, stringBuilder.length());
                for (int i = 0; i < base; i++) {
                    if (bits[i]) {
                        stringBuilder.append(1);
                    } else {
                        stringBuilder.append(0);
                    }
                }
                stringBuilder.append(':').append(alleleNum);
                int[] alleles = GlobalManager.codingPhasedGtyCodingMap.get(stringBuilder.toString());
                // String infor = "Sorry!!! squence variants with over 4 alleles are not supported and will be ignored!";
                // System.out.println(infor);
                return alleles;
        }
        return null;
    }

    public static final int[] getPhasedGtyBool(boolean[] gtyBits, int alleleNum, int base, int indivID) {
        switch (base) {
            case 2:
                /*       
                 missing	Reference homozygous	Heterozygous 	Heterozygous 	Alternative homozygous
                 VCF genotype		0|0	0|1	1|0	1|1 .|.
                 Bits	        000  	001	010	011	100
                 Order	0	1	2	3	4                
               
                 II.II Tri-allelic sequence variant (4 bits)
                 missing 	Reference homozygous 	Heterozygous 	Heterozygous 	Heterozygous 	Heterozygous 	Alternative homozygous
                 VCF genotype 	0|0 	0|1 	0|2 	1|0 	1|1 	1|2
                 Bits      	000 	0001 	0010 	0011 	0100 	0101 	0110
                 Decimal 	0 	1 	2 	3 	4 	5 	6
                 Heterozygous 	Heterozygous 	Alternative homozygous
                 VCF genotype 	2|0 	2|1 	2|2 	.|. 
                 Bits     	0111 	1000 	1001
                 Decimal 	7 	8 	9     
                 */
                if (!gtyBits[0] && !gtyBits[1] && !gtyBits[2]) {
                    return new int[]{0, 0};
                } else if (!gtyBits[0] && !gtyBits[1] && gtyBits[2]) {
                    return new int[]{0, 1};
                } else if (!gtyBits[0] && gtyBits[1] && !gtyBits[2]) {
                    return new int[]{1, 0};
                } else if (!gtyBits[0] && gtyBits[1] && gtyBits[2]) {
                    return new int[]{1, 1};
                } else if (gtyBits[0] && !gtyBits[1] && !gtyBits[2]) {
                    return null;
                }
                break;
            default:
                stringBuilder.delete(0, stringBuilder.length());
                for (int i = 0; i < base; i++) {
                    if (gtyBits[i]) {
                        stringBuilder.append(1);
                    } else {
                        stringBuilder.append(0);
                    }
                }
                stringBuilder.append(':').append(alleleNum);
                int[] alleles = GlobalManager.codingPhasedGtyCodingMap.get(stringBuilder.toString());
                // String infor = "Sorry!!! squence variants with over 4 alleles are not supported and will be ignored!";
                // System.out.println(infor);
                return alleles;
        }
        return null;
    }

    public static final void getUnphasedGtyAt1(byte[] gtys, int base, int indivID, boolean[] bits1, int blockSize) {
        bitNum = indivID;
        for (int i = 0; i < base; i++) {
            byteIndex1 = (bitNum) / 8;
            byteIndex2 = (bitNum) % 8;
            bits1[i] = (gtys[byteIndex1] & GlobalManager.byteOpers[byteIndex2]) == GlobalManager.byteOpers[byteIndex2];
            bitNum += blockSize;
        }
    }

    public BinaryGtyProcessor(String pedigreeFileName, String mapFileName, String kggseqBinaryFileName) {
        this.pedigreeFileName = pedigreeFileName;
        this.mapFileName = mapFileName;
        this.kggseqBinaryFileName = kggseqBinaryFileName;
    }

    public boolean avaibleFiles() {
        File file = new File(pedigreeFileName);
        if (!file.exists()) {
            return false;
        }
        file = new File(mapFileName);
        if (!file.exists()) {
            return false;
        }
        file = new File(kggseqBinaryFileName);
        if (!file.exists()) {
            return false;
        }
        return true;
    }

    public void readPedigreeFile(List<Individual> indivList) throws Exception {
        StringBuilder tmpBuffer = new StringBuilder();
        String line = null;
        String delimiter = "\t\" \",/";
        BufferedReader br = null;

        File file = new File(pedigreeFileName + ".gz");
        if (file.exists()) {
            br = LocalFileFunc.getBufferedReader(file.getCanonicalPath());
            needGzExtension = true;
        } else {
            file = new File(pedigreeFileName);
            if (file.exists()) {
                br = LocalFileFunc.getBufferedReader(file.getCanonicalPath());
            } else {
                LOG.error(file.getCanonicalPath() + " does not exist!");
                return;
            }
        }

        while ((line = br.readLine()) != null) {
            line = line.toUpperCase();
            StringTokenizer tokenizer = new StringTokenizer(line, delimiter);
            Individual indiv = new Individual();
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setFamilyID(tmpBuffer.toString());
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setIndividualID(tmpBuffer.toString());
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setDadID(tmpBuffer.toString());
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setMomID(tmpBuffer.toString());
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setGender(Integer.valueOf(tmpBuffer.toString()));
            // indiv.setLabelInChip(indiv.getFamilyID() + "@*@" + indiv.getIndividualID());
            indiv.setLabelInChip(indiv.getIndividualID());
            if (tokenizer.hasMoreTokens()) {
                indiv.setAffectedStatus(Integer.parseInt(tokenizer.nextToken().trim()));
            }

            //System.out.println(indiv.getLabelInChip());
            indivList.add(indiv);
        }
        br.close();
    }

    // public void calcualte
    public StringBuilder readBinaryGenotype(List<Individual> subjectList, Genome genome, int minOBS, int maxGtyAlleleNum, double sampleMafOver,
            double sampleMafLess, boolean considerSNP, boolean considerIndel) throws Exception {
        int indiSize = subjectList.size();
        IntArrayList caseSetID = new IntArrayList();
        IntArrayList controlSetID = new IntArrayList();

        for (int i = 0; i < indiSize; i++) {
            if (subjectList.get(i).getAffectedStatus() == 2) {
                caseSetID.add(i);
            } else if (subjectList.get(i).getAffectedStatus() == 1) {
                controlSetID.add(i);
            }
        }
        int subID = 0;
        int caseNum = caseSetID.size();
        int controlNum = controlSetID.size();

        String info = ("Reading genotype bit-file from [ " + kggseqBinaryFileName + " ] \n");
        //openkggseqPedFormat
        DataInputStream in = null;
        if (needGzExtension) {
            in = new DataInputStream(new BufferedInputStream(new GZIPInputStream(new FileInputStream(kggseqBinaryFileName + ".gz"))));
        } else {
            in = new DataInputStream(new BufferedInputStream(new FileInputStream(kggseqBinaryFileName)));
        }

        byte bt1 = in.readByte();
        byte bt2 = in.readByte();
        // System.out.println(bytesToHexString(new byte[]{bt1}));
        // System.out.println(bytesToHexString(new byte[]{bt2}));
        // printBitSet(b); 
        // Magic numbers for .ked file: 10011110 10000010  =9E 82
        if (bt1 != -98 || bt2 != -126) {
            throw new Exception("The " + kggseqBinaryFileName + " is not a valid kggseq binary file!!!");
        }
        boolean isPhased = false;

        bt1 = in.readByte();
        if (bt1 == 1) {
            isPhased = true;
            genome.setIsPhasedGty(true);
        }

        int indiviNum = subjectList.size();
        int alleleNum = 0;
        Chromosome[] chroms = genome.getChromosomes();
        int base = 0;
        int bitNum = 0;
        int intNum = 0;
        int[] gtys = null;
        boolean needAccoundAffect = false;
        boolean needAccoundUnaffect = false;
        int g11 = 0;
        int g12 = 0;
        int g22 = 0;
        int missingG = 0;
        if (caseNum > 0) {
            needAccoundAffect = true;
        }
        if (controlNum > 0) {
            needAccoundUnaffect = true;
        }

        List<Variant> tmpList = new ArrayList<Variant>();
        int ignoredLineNumMinOBS = 0, ignoredLineNumMinMAF = 0, ignoredLineNumMaxMAF = 0, ignoredVarBymaxGtyAlleleNum = 0;
        int obsS = 0;

        double sampleMafOverC = 1 - sampleMafOver;
        double sampleMafLessC = 1 - sampleMafLess;
        boolean needMAFQCOver = false;
        if (sampleMafOver >= 0) {
            needMAFQCOver = true;
        }
        boolean needMAFQCLess = false;
        if (sampleMafLess < 0.5) {
            needMAFQCLess = true;
        }
        boolean incomplete;
        double maf;
        final int gtyLen = 8;
        int codeI;
        byte[] tmpBytes;

        for (int chromID = 0; chromID < chroms.length; chromID++) {
            if (chroms[chromID] == null) {
                continue;
            }
            tmpList.clear();

            for (Variant var : chroms[chromID].variantList) {
                alleleNum = var.getAltAlleles().length + 1;

                if (isPhased) {
                    if (var.compressedGtyLabel < 0) {
                        base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
                        bitNum = base * indiviNum;
                        intNum = bitNum / gtyLen;
                        if (bitNum % gtyLen != 0) {
                            intNum++;
                        }
                    } else {
                        intNum = var.compressedGtyLabel;
                    }

                    var.encodedGty = new byte[intNum];
                    intNum = in.read(var.encodedGty);
                    if (intNum != var.encodedGty.length) {
                        String infor = "Error occurred when reading binary genotypes. Your binary genotypes may be broken!";
                        LOG.error(infor);
                        System.exit(1);
                    }
                } else {
                    if (var.compressedGtyLabel < 0) {
                        base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
                        bitNum = base * indiviNum;
                        intNum = bitNum / gtyLen;
                        if (bitNum % gtyLen != 0) {
                            intNum++;
                        }
                    } else {
                        intNum = var.compressedGtyLabel;
                    }
                    var.encodedGty = new byte[intNum];
                    intNum = in.read(var.encodedGty);
                    if (intNum != var.encodedGty.length) {
                        String infor = "Error occurred when reading binary genotypes. Your binary genotypes may be broken!";
                        LOG.error(infor);
                        System.exit(1);
                    }
                }
                //In any case, it must read the binary genotype data
                if (alleleNum > maxGtyAlleleNum) {
                    ignoredVarBymaxGtyAlleleNum++;
                    continue;
                }
                if (!considerSNP || !considerIndel) {
                    //a lazy point 
                    incomplete = true;

                    //only consider Indel
                    if (!considerSNP && var.isIndel) {
                        incomplete = false;
                    } else if (!considerIndel && !var.isIndel) {
                        incomplete = false;
                    }

                    if (incomplete) {
                        continue;
                    }
                }

                if (var.compressedGtyLabel >= 0) {
                    if (isPhased) {
                        base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
                    } else {
                        base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
                    }
                    bitNum = base * indiviNum;
                    intNum = bitNum / gtyLen;
                    if (bitNum % gtyLen != 0) {
                        intNum++;
                    }

                    tmpBytes = new byte[intNum];
                    Arrays.fill(tmpBytes, (byte) 0);
                    for (int j = 0; j < var.compressedGtyLabel; j += 3) {
                        codeI = (var.encodedGty[j + 2] & 0xFF) | ((var.encodedGty[j + 1] & 0xFF) << 8) | ((var.encodedGty[j] & 0x0F) << 16);
                        byteIndex1 = (codeI) / 8;
                        byteIndex2 = (codeI) % 8;
                        tmpBytes[byteIndex1] = (byte) (tmpBytes[byteIndex1] | GlobalManager.byteOpers[byteIndex2]);
                    }
                    var.encodedGty = tmpBytes;
                }

                obsS = 0;
                g11 = 0;
                g12 = 0;
                g22 = 0;
                missingG = 0;
                for (int j = 0; j < indiSize; j++) {
                    subID = j;
                    if (isPhased) {
                        gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, indiSize);
                    } else {
                        gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, indiSize);
                    }
                    if (gtys == null) {
                        missingG++;
                        continue;
                    }
                    if (gtys[0] != gtys[1]) {
                        g12++;
                    } else if (gtys[0] == 0) {
                        g11++;
                    } else {
                        g22++;
                    }
                    obsS++;
                }
                var.setAffectedRefHomGtyNum(g11);
                var.setAffectedHetGtyNum(g12);
                var.setAffectedAltHomGtyNum(g22);
                var.setMissingtyNum(missingG);
                if (obsS < minOBS) {
                    ignoredLineNumMinOBS++;
                    continue;
                }
                maf = (g12 * 0.5 + g22) / (g11 + g12 + g22);
                if (needMAFQCOver || needMAFQCLess) {
                    if (needMAFQCOver) {
                        if (Double.isNaN(maf) || maf <= sampleMafOver || maf >= sampleMafOverC) {
                            ignoredLineNumMinMAF++;
                            continue;
                        }
                    }
                    if (needMAFQCLess) {
                        if (Double.isNaN(maf) || maf >= sampleMafLess && maf <= sampleMafLessC) {
                            ignoredLineNumMaxMAF++;
                            continue;
                        }
                    }
                }

                var.localAltAF = (float) maf;
                missingG = 0;
                if (needAccoundAffect) {
                    g11 = 0;
                    g12 = 0;
                    g22 = 0;
                    for (int j = 0; j < caseNum; j++) {
                        subID = caseSetID.getQuick(j);
                        if (isPhased) {
                            gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, indiSize);
                        } else {
                            gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, indiSize);
                        }
                        if (gtys == null) {
                            missingG++;
                            continue;
                        }
                        if (gtys[0] != gtys[1]) {
                            g12++;
                        } else if (gtys[0] == 0) {
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
                    for (int i = 0; i < controlNum; i++) {
                        subID = controlSetID.getQuick(i);
                        if (isPhased) {
                            gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, indiSize);
                        } else {
                            gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, indiSize);
                        }
                        if (gtys == null) {
                            missingG++;
                            continue;
                        }
                        if (gtys[0] != gtys[1]) {
                            g12++;
                        } else if (gtys[0] == 0) {
                            g11++;
                        } else {
                            g22++;
                        }
                    }
                    var.setUnaffectedRefHomGtyNum(g11);
                    var.setUnaffectedHetGtyNum(g12);
                    var.setUnaffectedAltHomGtyNum(g22);
                    var.setMissingtyNum(missingG);
                }
                tmpList.add(var);
            }
            chroms[chromID].variantList.clear();
            chroms[chromID].variantList.addAll(tmpList);
            tmpList.clear();
        }
        in.close();
        StringBuilder message = new StringBuilder();
        message.append("Quality control summaries:\n");
        if (ignoredLineNumMinOBS > 0) {
            message.append(" ").append(ignoredLineNumMinOBS).append(" variants are ignored due to the number of non-null genotypes in sample <").append(minOBS).append('\n');
        }

        if (ignoredLineNumMinMAF > 0) {
            message.append(" ").append(ignoredLineNumMinMAF).append(" variants are ignored due to their minor allele frequency (MAF) in sample <=").append(sampleMafOver).append('\n');
        }
        if (ignoredLineNumMaxMAF > 0) {
            message.append(" ").append(ignoredLineNumMaxMAF).append(" variants are ignored due to their minor allele frequency (MAF) in sample >=").append(sampleMafLess).append('\n');
        }
        if (ignoredVarBymaxGtyAlleleNum > 0) {
            message.append(" ").append(ignoredVarBymaxGtyAlleleNum).append(" variants are ignored because the number of alleles is > ").append(maxGtyAlleleNum).append(";\n");
        }

        return message;
    }

    public static String bytesToHexString(byte[] src) {
        StringBuilder stringBuilder = new StringBuilder("");
        if (src == null || src.length <= 0) {
            return null;
        }
        for (int i = 0; i < src.length; i++) {
            int v = src[i] & 0xFF;
            String hv = Integer.toHexString(v);
            if (hv.length() < 2) {
                stringBuilder.append(0);
            }
            stringBuilder.append(hv);
        }
        return stringBuilder.toString();
    }

    public Genome readVariantsMapFile(int[] counts) throws Exception {
        BufferedReader br = null;
        File file = new File(mapFileName + ".gz");
        if (file.exists()) {
            br = LocalFileFunc.getBufferedReader(file.getCanonicalPath());
        } else {
            file = new File(mapFileName);
            if (file.exists()) {
                br = LocalFileFunc.getBufferedReader(file.getCanonicalPath());
            } else {
                LOG.error(file.getCanonicalPath() + " does not exist!");
                return null;
            }
        }

        String line = null;

        int chromIndex = 0;
        int labelIndex = 1;
        int positionIndex = 3;
        int refIndex = 4;
        int altIndex = 5;
        int codeIndex = 6;
        int maxIndex = codeIndex;

        Genome genome = new Genome("KGGSeqGenome", "ked");
        genome.removeTempFileFromDisk();
        String chrom;
        int position = -1;
        String label = null;
        String ref = null;
        String alt = null;
        String code = null;
        StringBuilder tmpBuffer = new StringBuilder();
        int index;
        int effectiveSNPNum = 0;
        int lineCounter = 0;
        int indelNum = 0;
        try {
            while ((line = br.readLine()) != null) {
                if (line.trim().length() == 0) {
                    continue;
                }
                lineCounter++;
                StringTokenizer tokenizer = new StringTokenizer(line);
                index = 0;

                chrom = null;
                position = -1;
                label = null;
                ref = null;
                alt = null;
                code = null;
                while (tokenizer.hasMoreTokens()) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    tmpBuffer.append(tokenizer.nextToken().trim());

                    if (index == chromIndex) {
                        chrom = tmpBuffer.toString();
                    } else if (index == positionIndex) {
                        position = Integer.parseInt(tmpBuffer.toString());
                    } else if (index == labelIndex) {
                        label = tmpBuffer.toString();
                    } else if (index == refIndex) {
                        ref = tmpBuffer.toString();
                    } else if (index == altIndex) {
                        alt = tmpBuffer.toString();
                    } else if (index == codeIndex) {
                        code = tmpBuffer.toString();
                    }

                    if (index == maxIndex) {
                        break;
                    }
                    index++;
                }

                effectiveSNPNum++;
                Variant var = new Variant(position, ref, alt.split(","));
                if (alt.indexOf('+') >= 0 || alt.indexOf('-') >= 0) {
                    indelNum++;
                    var.isIndel = true;
                }

                var.compressedGtyLabel = Integer.parseInt(code);
                var.setLabel(label);
                genome.addVariant(chrom, var);
            }

            counts[0] = lineCounter;
            counts[1] = effectiveSNPNum;
            counts[2] = indelNum;

            /*
             StringBuilder runningInfo = new StringBuilder();
             runningInfo.append("The number of SNPs  in map file ");
             runningInfo.append(mapFile.getName());
             runningInfo.append(" is ");
             runningInfo.append(effectiveSNPNum);
             runningInfo.append(".");
             */
        } finally {
            br.close();
        }
        genome.setVarNum(effectiveSNPNum);
        return genome;
    }

    public static void main(String[] args) {
        try {
            //  byte byteInfo=-128;
            //  System.out.println(BitByteUtil.byteToBinaryString((byte) (byteInfo )));
            //                   System.out.println(BitByteUtil.byteToBinaryString((byte) (byteInfo>>> 2&0X3F)));                             
            List<Individual> indivList = new ArrayList<Individual>();
            BinaryGtyProcessor bgp = new BinaryGtyProcessor("./kggseq");
            bgp.readPedigreeFile(indivList);
            Genome genome = bgp.readVariantsMapFile(null);
            // bgp.readBinaryGenotype(indivList, genome);
            //   genome.export2FlatTextPlink(indivList, "./test");
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }
}
