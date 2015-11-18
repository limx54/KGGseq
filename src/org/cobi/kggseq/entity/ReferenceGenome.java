/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.entity;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.cobi.kggseq.Constants;
import org.cobi.kggseq.GlobalManager;
import org.cobi.util.file.LocalFileFunc;

/**
 *
 * @author mxli
 */
public class ReferenceGenome implements Constants {

    private static final Logger LOG = Logger.getLogger(ReferenceGenome.class);
    Map<String, Byte> fullchromNameIndexMap = new HashMap<String, Byte>();
    Map<String, Byte> chromNameIndexMap = new HashMap<String, Byte>();
    private List[] chromosomes = null;
    private List[] chromosomemRNABoundaryIndexList = null;
    private Map<String, int[]> geneChromPosMap = new HashMap<String, int[]>();
    boolean hasNotSorted = true;
    int upstreamDis = 1000;
    int donwstreamDis = 1000;
    int splicingDis = 2;
    private String name;
//    static final RNABoundaryIndex searchedmRNABoundaryIndex = new RNABoundaryIndex(0);
    static final RNABoundaryIndexComparator rNAGeneBoundaryIndexComparator = new RNABoundaryIndexComparator();

    public ReferenceGenome() {
        chromosomes = new ArrayList[STAND_CHROM_NAMES.length];
        chromosomemRNABoundaryIndexList = new ArrayList[STAND_CHROM_NAMES.length];
        for (byte i = 0; i < STAND_CHROM_NAMES.length; i++) {
            fullchromNameIndexMap.put("chr" + STAND_CHROM_NAMES[i], i);
            chromNameIndexMap.put(STAND_CHROM_NAMES[i], i);
        }
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public ReferenceGenome(int spd, int usd, int dsd) {
        splicingDis = spd;
        upstreamDis = usd;
        donwstreamDis = dsd;
        chromosomes = new ArrayList[STAND_CHROM_NAMES.length];
        chromosomemRNABoundaryIndexList = new ArrayList[STAND_CHROM_NAMES.length];
        for (byte i = 0; i < STAND_CHROM_NAMES.length; i++) {
            fullchromNameIndexMap.put("chr" + STAND_CHROM_NAMES[i], i);
            chromNameIndexMap.put(STAND_CHROM_NAMES[i], i);
        }
    }

    //a function to read the dataset refGeneMrna.fa.gz prepared by UCSC
    public void readRefGeneSequence(String sequenceFile) throws Exception {
        //Note: a transcipt ID can be mapped onto multiple regions
        Map<String, String> transciptSequence = new HashMap<String, String>();
        File dataFile = new File(sequenceFile);
        BufferedReader br = LocalFileFunc.getBufferedReader(sequenceFile);
        String currentLine;
        int lineCounter = 0;
        String lastTranscripID = null;
        String chromIDPostion;
        int index = 0;
        int index1 = 0;
        int strLen = 0;
        StringBuilder sb = new StringBuilder();

        while ((currentLine = br.readLine()) != null) {
            lineCounter++;

            if (currentLine.startsWith(">")) {

                //System.out.println(lastTranscripID);
                // >NM_001101330 1
//atggtttccgctagtgggacatcattttttaagggtatgttgcttgggag
//catttcctgggttttgataactatgtttggccaaattcacattcgacaca
//gaggtcagactcaagaccacgagcaccatcaccttcgtccacctaacagg
                if (lastTranscripID != null) {
                    String mrnaSeq = transciptSequence.get(lastTranscripID);
                    if (mrnaSeq == null) {
                        transciptSequence.put(lastTranscripID, sb.toString().toUpperCase());
                    } else {
                        System.out.println(lastTranscripID + " has more than 1 sequence");
                    }
                    sb.delete(0, sb.length());
                }
                index = currentLine.indexOf(' ');
                lastTranscripID = currentLine.substring(1, index);

            } else {
                sb.append(currentLine);
            }
        }
        //the last transcript
        if (sb.length() > 0) {
            String mrnaSeq = transciptSequence.get(lastTranscripID);
            sb.delete(0, sb.length());
            if (mrnaSeq == null) {
                transciptSequence.put(lastTranscripID, sb.toString().toUpperCase());
            } else {
                System.out.println(lastTranscripID + " has more than 1 sequence");
            }
        }
        br.close();
        String info = "Read sequence of " + transciptSequence.size() + " transcripts.";
        System.out.println(info);
        int noSequenceRNA = 0;
        StringBuilder noSequenceRNAIDs = new StringBuilder();

        //assign the sequence to genes
        for (int iChrom = 0; iChrom < chromosomes.length; iChrom++) {
            List<RefmRNA> chrom = chromosomes[iChrom];
            if (chrom == null) {
                continue;
            }
            for (int imRNA = 0; imRNA < chrom.size(); imRNA++) {
                RefmRNA mrna = chrom.get(imRNA);
                if (mrna.codingStart == mrna.codingEnd) {
                    String mixID = mrna.getRefID() + ":@hr" + STAND_CHROM_NAMES[iChrom] + ":" + mrna.getStart();
                    noSequenceRNA++;
                    noSequenceRNAIDs.append(mixID);
                    noSequenceRNAIDs.append(" ");
                    continue;
                }

                //keep all mRNA whether it has sequence or not
                String msqSeq = transciptSequence.get(mrna.getRefID());
                if (msqSeq == null) {
                    String mixID = mrna.getRefID() + "@chr" + STAND_CHROM_NAMES[iChrom] + ":" + mrna.getStart();
                    noSequenceRNA++;
                    noSequenceRNAIDs.append(mixID);
                    noSequenceRNAIDs.append(" ");
                    continue;
                }
                mrna.setmRnaSequenceStart(mrna.getStart());
                mrna.setmRnaSequence(msqSeq);
            }
        }
        transciptSequence.clear();
        if (noSequenceRNA > 0) {
            info = noSequenceRNA + " transcripts have no sequence data.";
            System.out.println(info);
            //System.out.println(noSequenceRNAIDs);
        }
    }

    public RefmRNA getmRNA(int[] poss) {
        return (RefmRNA) chromosomes[poss[0]].get(poss[1]);
    }

    public RefCNV getCNV(int[] poss) {
        return (RefCNV) chromosomes[poss[0]].get(poss[1]);
    }

    public RefDup getDup(int[] poss) {
        return (RefDup) chromosomes[poss[0]].get(poss[1]);
    }

    public RefmRNA getmRNA(String syb) {
        int[] poss = geneChromPosMap.get(syb);
        if (poss == null) {
            return null;
        }
        return (RefmRNA) chromosomes[poss[0]].get(poss[1]);
    }

    public int[] getmRNAPos(String syb) {
        return geneChromPosMap.get(syb);
    }

    public void exportGeneRegions(String outPath) throws Exception {
        int upstreamExtendLen = 5000;
        int downstreamExtendLen = 5000;
        BufferedWriter bw = new BufferedWriter(new FileWriter(outPath));
        bw.write("--regions  ");
        for (int iChrom = 0; iChrom < STAND_CHROM_NAMES.length; iChrom++) {
            List<RefmRNA> chrom = chromosomes[iChrom];
            if (chrom == null) {
                continue;
            }
            for (int iGene = 0; iGene < chrom.size(); iGene++) {
                RefmRNA mrna = chrom.get(iGene);

                if (mrna.codingStart == mrna.codingEnd) {
                    //   continue;
                }
                int start = (mrna.getStart() - upstreamExtendLen);
                if (start < 0) {
                    start = 0;
                }
                String mixID = "chr" + STAND_CHROM_NAMES[iChrom] + ":"
                        + start + "-" + (mrna.getEnd() + downstreamExtendLen) + ",";
                bw.write(mixID);
                // System.out.println(mrna.getRefID() + "\t" + mrna.getStart() + "\t" + mrna.getEnd());

            }
        }
        bw.close();

    }

    //as RefmRNA can be mappled onto multiple chrosomes so it would be very difficult to use a mrna to orgnaize the RefmRNA
    public void addRefRNA(RefmRNA mrna, String chrom) {
        Byte chromID = chromNameIndexMap.get(chrom);
        if (chromID == null) {
            return;
        }
        if (chromosomes[chromID] == null) {
            chromosomes[chromID] = new ArrayList<RefmRNA>();
            chromosomes[chromID].add(mrna);
        } else {
            chromosomes[chromID].add(mrna);
        }

        geneChromPosMap.put(mrna.getRefID() + ":" + chrom + ":" + mrna.codingStart + ":" + mrna.codingEnd, new int[]{chromID, chromosomes[chromID].size() - 1});
        hasNotSorted = true;
    }

    //as RefmRNA can be mappled onto multiple chrosomes so it would be very difficult to use a mrna to orgnaize the RefmRNA
    public void addRefCNV(RefCNV cnv, String chrom) {
        Byte chromID = chromNameIndexMap.get(chrom);
        if (chromID == null) {
            return;
        }
        if (chromosomes[chromID] == null) {
            chromosomes[chromID] = new ArrayList<RefCNV>();
            chromosomes[chromID].add(cnv);
        } else {
            chromosomes[chromID].add(cnv);
        }

        geneChromPosMap.put(cnv.getDescription() + ":" + chrom + ":" + cnv.getStart() + ":" + cnv.getEnd(), new int[]{chromID, chromosomes[chromID].size() - 1});
        hasNotSorted = true;
    }

    //as RefmRNA can be mappled onto multiple chrosomes so it would be very difficult to use a mrna to orgnaize the RefmRNA
    public void addRefDup(RefDup cnv, String chrom) {
        Byte chromID = chromNameIndexMap.get(chrom);
        if (chromID == null) {
            return;
        }
        if (chromosomes[chromID] == null) {
            chromosomes[chromID] = new ArrayList<RefCNV>();
            chromosomes[chromID].add(cnv);
        } else {
            chromosomes[chromID].add(cnv);
        }

        geneChromPosMap.put(cnv.getDescription() + ":" + chrom + ":" + cnv.getStart() + ":" + cnv.getEnd(), new int[]{chromID, chromosomes[chromID].size() - 1});
        hasNotSorted = true;
    }
    /*
     * A reference from Annovar
     Feature	Value 	Explanation
     nonsynonymous	1	Variants result in a codon coding for a different amino acid (missense) and an amino acid codon to a stop codon (stopgain) and a stop codon to an amino acid codon (stoplos)
     synonymous	2	
     splicing	3	variant is within 2-bp of a splicing junction (use -splicing_threshold to change this) 
     ncRNA	4	variant overlaps a transcript without coding annotation in the region definition (see Notes below for more explanation) 
     5UTR	5	variant overlaps a 5' untranslated region 
     3UTR	6	variant overlaps a 3' untranslated region 
     intronic	7	variant overlaps an intron 
     upstream	8	variant overlaps 1-kb region upstream of transcription start site 
     downstream	9	variant overlaps 1-kb region downtream of transcription end site (use -neargene to change this) 
     intergenic	10	variant is in intergenic region     
     */

    public GeneFeature getVarFeature(String chrom, Variant var, boolean isForwardStrandInput, Set<Byte> priorFeatureSet, RNABoundaryIndex searchedmRNABoundaryIndex) throws Exception {
        Byte chromID = chromNameIndexMap.get(chrom);
        if (chromID == null) {
            throw new Exception("Unknown chrosomsome ID " + chrom);
        }

        if (hasNotSorted) {
            throw new Exception("The genome has to been sorted before search the gene-feature of a variant.");
        }
        List<GeneFeature> featureList = new ArrayList<GeneFeature>();
        // System.out.println(chrom + " : " + var.physicalPosition);
        byte intergID = GlobalManager.VarFeatureIDMap.get("intergenic");

        if (chromosomemRNABoundaryIndexList[chromID] == null) {
            // String info = "Variant at chr" + chrom + ":" + var.refStartPosition + " cannot be annotated by the specified gene annotation database " + name + "!";
            //System.out.println(info);
            //  LOG.warn(info);
            if (var.smallestFeatureID > GlobalManager.VarFeatureIDMap.get("unknown")) {
                var.smallestFeatureID = GlobalManager.VarFeatureIDMap.get("unknown");
            }
            return new GeneFeature(GlobalManager.VarFeatureIDMap.get("unknown"), null);
        }

        /*
         if (chrom.equals("6")){
         var.refStartPosition=	31322303;
         var.setRefAllele('C');
         var.setAltAlleles(new String[]{"G"});
         }
         * 
         */
        //note for forward and reverse strand, the upstream and downsteam could be different
        int pos = var.refStartPosition;
        searchedmRNABoundaryIndex.position = pos;

        // System.out.println(pos);
        int headIndex = 0;
        int tailIndex = Collections.binarySearch(chromosomemRNABoundaryIndexList[chromID], searchedmRNABoundaryIndex, rNAGeneBoundaryIndexComparator);
        if (tailIndex < 0) {
            tailIndex = -tailIndex - 1;
            if (tailIndex == chromosomemRNABoundaryIndexList[chromID].size()) {
                if (var.smallestFeatureID > intergID) {
                    var.smallestFeatureID = intergID;
                }
                return new GeneFeature(intergID, null);
            }
            headIndex = tailIndex - 1;
            if (headIndex < 0) {
                if (var.smallestFeatureID > intergID) {
                    var.smallestFeatureID = intergID;
                }
                return new GeneFeature(intergID, null);
            }

            RNABoundaryIndex rbi1 = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[chromID].get(headIndex);
            //impossible to be empty; it at least inlude the RefmRNA itselft
            /*
             if (rbi1.mRNAIndexList.isEmpty()) {
             var.geneFeatureID = intergID;
             return geneSymbs;
             }
             * 
             */

            RNABoundaryIndex rbi2 = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[chromID].get(tailIndex);
            //impossible to be empty; it at least inlude the RefmRNA itselft
            /*
             if (rbi2.mRNAIndexList.isEmpty()) {
             var.geneFeatureID = intergID;
             return geneSymbs;
             }
             * 
             */

            int size1 = rbi1.mRNAIndexList.size();
            int size2 = rbi2.mRNAIndexList.size();
            if (size1 < size2) {
                for (int i = 0; i < size1; i++) {
                    int searchIndex = rbi1.mRNAIndexList.getQuick(i);
                    if (rbi2.mRNAIndexList.contains(searchIndex)) {
                        RefmRNA mRNA = (RefmRNA) chromosomes[chromID].get(searchIndex);
                        GeneFeature gf = mRNA.findFeature(chrom, var, isForwardStrandInput, upstreamDis, donwstreamDis, splicingDis);
                        if (gf != null) {
                            //gf.setName(mRNA.geneSymb + ":" + gf.getName());
                            featureList.add(gf);
                        }
                    }
                }
            } else {
                for (int i = 0; i < size2; i++) {
                    int searchIndex = rbi2.mRNAIndexList.getQuick(i);
                    if (rbi1.mRNAIndexList.contains(searchIndex)) {
                        RefmRNA mRNA = (RefmRNA) chromosomes[chromID].get(searchIndex);
                        GeneFeature gf = mRNA.findFeature(chrom, var, isForwardStrandInput, upstreamDis, donwstreamDis, splicingDis);
                        if (gf != null) {
                            //gf.setName(mRNA.geneSymb + ":" + gf.getName());
                            featureList.add(gf);
                        }
                    }
                }
            }

        } else {
            RNABoundaryIndex rbi = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[chromID].get(tailIndex);
            if (rbi.mRNAIndexList.isEmpty()) {
                if (var.smallestFeatureID > intergID) {
                    var.smallestFeatureID = intergID;
                }
                return new GeneFeature(intergID, null);
            }
            for (int i = 0; i < rbi.mRNAIndexList.size(); i++) {
                int searchIndex = rbi.mRNAIndexList.getQuick(i);
                RefmRNA mRNA = (RefmRNA) chromosomes[chromID].get(searchIndex);
                GeneFeature gf = mRNA.findFeature(chrom, var, isForwardStrandInput, upstreamDis, donwstreamDis, splicingDis);
                if (gf != null) {
                    //gf.setName(mRNA.geneSymb + ":" + gf.getName());
                    featureList.add(gf);
                }
            }
        }

        if (featureList.isEmpty()) {
            if (var.smallestFeatureID > intergID) {
                var.smallestFeatureID = intergID;
            }
            return new GeneFeature(intergID, null);
        } else {
            int gfSize = featureList.size();
            if (gfSize == 1) {
                if (var.smallestFeatureID > featureList.get(0).id) {
                    var.smallestFeatureID = featureList.get(0).id;
                    if (featureList.get(0).id < intergID) {
                        String ftName = featureList.get(0).name;
                        int index = ftName.indexOf(':');
                        if (index >= 0) {
                            var.geneSymb = ftName.substring(0, index);
                        }
                    }
                }
                return featureList.get(0);
            } else {
                Collections.sort(featureList, new GeneFeatureComparator());
                //by default use the first description 
                boolean hasAssignSeledtedID = false;
                GeneFeature selectGf = featureList.get(0);
                if (priorFeatureSet.contains(selectGf.id)) {
                    hasAssignSeledtedID = true;
                }

                if (var.smallestFeatureID > selectGf.id) {
                    var.smallestFeatureID = selectGf.id;
                    if (featureList.get(0).id < intergID) {
                        String ftName = featureList.get(0).name;
                        int index = ftName.indexOf(':');
                        if (index >= 0) {
                            var.geneSymb = ftName.substring(0, index);
                        }
                    }
                }

                for (int i = 1; i < gfSize; i++) {
                    GeneFeature gf = featureList.get(i);
                    if (!hasAssignSeledtedID && priorFeatureSet.contains(gf.id)) {
                        selectGf.id = gf.id;
                    }
                    if (!gf.getName().contains("intergenic")) {
                        selectGf.setName(gf.getName() + ";" + selectGf.getName());
                    }
                }
                return selectGf;
            }
        }
    }

    public List<RefCNV> getCNVFeature(String chrom, Variant var, RNABoundaryIndex searchedmRNABoundaryIndex) throws Exception {
        Byte chromID = chromNameIndexMap.get(chrom);
        if (chromID == null) {
            throw new Exception("Unknown chrosomsome ID " + chrom);
        }

        if (hasNotSorted) {
            throw new Exception("The genome has to been sorted before search the gene-feature of a variant.");
        }
        List<RefCNV> cnvList = new ArrayList<RefCNV>();
        // System.out.println(chrom + " : " + var.physicalPosition);

        if (chromosomemRNABoundaryIndexList[chromID] == null) {
            //String info = "A variant at " + var.refStartPosition + " of chromosome " + chrom + " cannot be annotated!";
            //System.out.println(info);
            //LOG.warn(info);           
            return cnvList;
        }

        //note for forward and reverse strand, the upstream and downsteam could be different
        int pos = var.refStartPosition;
        searchedmRNABoundaryIndex.position = pos;

        // System.out.println(pos);
        int headIndex = 0;
        int tailIndex = Collections.binarySearch(chromosomemRNABoundaryIndexList[chromID], searchedmRNABoundaryIndex, rNAGeneBoundaryIndexComparator);
        if (tailIndex < 0) {
            tailIndex = -tailIndex - 1;
            if (tailIndex == chromosomemRNABoundaryIndexList[chromID].size()) {
                return cnvList;
            }
            headIndex = tailIndex - 1;
            if (headIndex < 0) {
                return cnvList;
            }

            RNABoundaryIndex rbi1 = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[chromID].get(headIndex);
            //impossible to be empty; it at least inlude the RefmRNA itselft
            /*
             if (rbi1.mRNAIndexList.isEmpty()) {
             var.geneFeatureID = intergID;
             return geneSymbs;
             }
             * 
             */

            RNABoundaryIndex rbi2 = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[chromID].get(tailIndex);
            //impossible to be empty; it at least inlude the RefmRNA itselft
            /*
             if (rbi2.mRNAIndexList.isEmpty()) {
             var.geneFeatureID = intergID;
             return geneSymbs;
             }
             * 
             */

            int size1 = rbi1.mRNAIndexList.size();
            int size2 = rbi2.mRNAIndexList.size();
            if (size1 < size2) {
                for (int i = 0; i < size1; i++) {
                    int searchIndex = rbi1.mRNAIndexList.getQuick(i);
                    if (rbi2.mRNAIndexList.contains(searchIndex)) {
                        RefCNV cnv = (RefCNV) chromosomes[chromID].get(searchIndex);
                        //gf.setName(mRNA.geneSymb + ":" + gf.getName());
                        cnvList.add(cnv);
                    }
                }
            } else {
                for (int i = 0; i < size2; i++) {
                    int searchIndex = rbi2.mRNAIndexList.getQuick(i);
                    if (rbi1.mRNAIndexList.contains(searchIndex)) {
                        RefCNV cnv = (RefCNV) chromosomes[chromID].get(searchIndex);
                        cnvList.add(cnv);
                    }
                }
            }

        } else {
            RNABoundaryIndex rbi = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[chromID].get(tailIndex);
            if (rbi.mRNAIndexList.isEmpty()) {
                return cnvList;
            }
            for (int i = 0; i < rbi.mRNAIndexList.size(); i++) {
                int searchIndex = rbi.mRNAIndexList.getQuick(i);
                RefCNV cnv = (RefCNV) chromosomes[chromID].get(searchIndex);
                cnvList.add(cnv);
            }
        }
        return cnvList;
    }

    public List<RefDup> getDupFeature(String chrom, Variant var, RNABoundaryIndex searchedmRNABoundaryIndex) throws Exception {
        Byte chromID = chromNameIndexMap.get(chrom);
        if (chromID == null) {
            throw new Exception("Unknown chrosomsome ID " + chrom);
        }

        if (hasNotSorted) {
            throw new Exception("The genome has to been sorted before search the gene-feature of a variant.");
        }
        List<RefDup> cnvList = new ArrayList<RefDup>();
        // System.out.println(chrom + " : " + var.physicalPosition);

        if (chromosomemRNABoundaryIndexList[chromID] == null) {
            //String info = "A variant at " + var.refStartPosition + " of chromosome " + chrom + " cannot be annotated!";
            //System.out.println(info);
            //LOG.warn(info);           
            return cnvList;
        }

        //note for forward and reverse strand, the upstream and downsteam could be different
        int pos = var.refStartPosition;
        searchedmRNABoundaryIndex.position = pos;

        // System.out.println(pos);
        int headIndex = 0;
        int tailIndex = Collections.binarySearch(chromosomemRNABoundaryIndexList[chromID], searchedmRNABoundaryIndex, rNAGeneBoundaryIndexComparator);
        if (tailIndex < 0) {
            tailIndex = -tailIndex - 1;
            if (tailIndex == chromosomemRNABoundaryIndexList[chromID].size()) {
                return cnvList;
            }
            headIndex = tailIndex - 1;
            if (headIndex < 0) {
                return cnvList;
            }

            RNABoundaryIndex rbi1 = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[chromID].get(headIndex);
            //impossible to be empty; it at least inlude the RefmRNA itselft
            /*
             if (rbi1.mRNAIndexList.isEmpty()) {
             var.geneFeatureID = intergID;
             return geneSymbs;
             }
             * 
             */

            RNABoundaryIndex rbi2 = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[chromID].get(tailIndex);
            //impossible to be empty; it at least inlude the RefmRNA itselft
            /*
             if (rbi2.mRNAIndexList.isEmpty()) {
             var.geneFeatureID = intergID;
             return geneSymbs;
             }
             * 
             */

            int size1 = rbi1.mRNAIndexList.size();
            int size2 = rbi2.mRNAIndexList.size();
            if (size1 < size2) {
                for (int i = 0; i < size1; i++) {
                    int searchIndex = rbi1.mRNAIndexList.getQuick(i);
                    if (rbi2.mRNAIndexList.contains(searchIndex)) {
                        RefDup cnv = (RefDup) chromosomes[chromID].get(searchIndex);
                        //gf.setName(mRNA.geneSymb + ":" + gf.getName());
                        cnvList.add(cnv);
                    }
                }
            } else {
                for (int i = 0; i < size2; i++) {
                    int searchIndex = rbi2.mRNAIndexList.getQuick(i);
                    if (rbi1.mRNAIndexList.contains(searchIndex)) {
                        RefDup cnv = (RefDup) chromosomes[chromID].get(searchIndex);
                        cnvList.add(cnv);
                    }
                }
            }

        } else {
            RNABoundaryIndex rbi = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[chromID].get(tailIndex);
            if (rbi.mRNAIndexList.isEmpty()) {
                return cnvList;
            }
            for (int i = 0; i < rbi.mRNAIndexList.size(); i++) {
                int searchIndex = rbi.mRNAIndexList.getQuick(i);
                RefDup cnv = (RefDup) chromosomes[chromID].get(searchIndex);
                cnvList.add(cnv);
            }
        }
        return cnvList;
    }

    private int binarySearchStartPos(int pos, int left, int right, int chromID) {
        if (left > right) {
            return -left - 1;
        }
        int middle = (left + right) / 2;
// the search is based on start posistion

        if (((RNASegment) chromosomemRNABoundaryIndexList[chromID].get(middle)).start == pos) {
            return middle;
        } else if (((RNASegment) chromosomemRNABoundaryIndexList[chromID].get(middle)).start > pos) {
            return binarySearchStartPos(pos, left, middle - 1, chromID);
        } else {
            return binarySearchStartPos(pos, middle + 1, right, chromID);
        }
    }

    private int binarySearchStartPos(int pos, int left, int right, List mRNABounds) {
        if (left > right) {
            return -left - 1;
        }
        int middle = (left + right) / 2;
// the search is based on start posistion

        if (((RNASegment) mRNABounds.get(middle)).start == pos) {
            return middle;
        } else if (((RNASegment) mRNABounds.get(middle)).start > pos) {
            return binarySearchStartPos(pos, left, middle - 1, mRNABounds);
        } else {
            return binarySearchStartPos(pos, middle + 1, right, mRNABounds);
        }
    }

    class RNASegment {

        protected int start;
        protected int end;
        protected int mRNAIndex;
        protected char strand;

        public RNASegment(int start, int end, int mRNAIndex, char strand) {
            this.start = start;
            this.end = end;

            this.mRNAIndex = mRNAIndex;
            this.strand = strand;
        }

        public int getEnd() {
            return end;
        }

        public int getmRNAIndex() {
            return mRNAIndex;
        }

        public int getStart() {
            return start;
        }
    }

    class GeneRNASegmentComparator implements Comparator<RNASegment> {

        @Override
        public int compare(RNASegment arg0, RNASegment arg1) {
            int result = -1;
            if (arg0.start == arg1.start) {
                result = arg0.end - arg1.end;
            } else {
                result = arg0.start - arg1.start;
            }
            return result;
        }
    }

    public void sortmRNAMakeIndexonChromosomes() {
        geneChromPosMap.clear();

        int size = 0;
        List<RNASegment> rNAPositionIndexList = new ArrayList<RNASegment>();

        for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
            if (chromosomes[i] != null) {
                Collections.sort(chromosomes[i], new SeqSegmentComparator());
                size = chromosomes[i].size();
                rNAPositionIndexList.clear();
                if (chromosomemRNABoundaryIndexList[i] == null) {
                    chromosomemRNABoundaryIndexList[i] = new ArrayList<RNABoundaryIndex>();
                } else {
                    chromosomemRNABoundaryIndexList[i].clear();
                }
                for (int j = 0; j < size; j++) {
                    RefmRNA mrna = (RefmRNA) chromosomes[i].get(j);
                    geneChromPosMap.put(mrna.getRefID() + ":" + STAND_CHROM_NAMES[i] + ":" + mrna.codingStart + ":" + mrna.codingEnd, new int[]{i, j});
                    if (mrna.getStrand() == '-') {
                        RNASegment gns = new RNASegment(mrna.start - donwstreamDis, mrna.end + upstreamDis, j, mrna.getStrand());
                        rNAPositionIndexList.add(gns);

                        chromosomemRNABoundaryIndexList[i].add(new RNABoundaryIndex(mrna.start - donwstreamDis));
                        chromosomemRNABoundaryIndexList[i].add(new RNABoundaryIndex(mrna.end + upstreamDis));
                    } else {
                        //default is +
                        RNASegment gns = new RNASegment(mrna.start - donwstreamDis, mrna.end + donwstreamDis, j, mrna.getStrand());
                        rNAPositionIndexList.add(gns);
                        chromosomemRNABoundaryIndexList[i].add(new RNABoundaryIndex(mrna.start - upstreamDis));
                        chromosomemRNABoundaryIndexList[i].add(new RNABoundaryIndex(mrna.end + donwstreamDis));
                    }
                }

                // the chromosomemRNABoundaryIndexList will be used to map a variant
                Collections.sort(rNAPositionIndexList, new GeneRNASegmentComparator());
                Collections.sort(chromosomemRNABoundaryIndexList[i], new RNABoundaryIndexComparator());
                int mrnaSize = rNAPositionIndexList.size();

                int boundSize = chromosomemRNABoundaryIndexList[i].size();
                for (int j = 0; j < boundSize; j++) {
                    RNABoundaryIndex rbi = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[i].get(j);

                    int pos = rbi.position;
                    int genStartIndex = binarySearchStartPos(pos, 0, mrnaSize - 1, rNAPositionIndexList);
                    if (genStartIndex < 0) {
                        genStartIndex = -genStartIndex - 1;
                    }
                    if (genStartIndex >= mrnaSize) {
                        genStartIndex = mrnaSize - 1;
                    }
                    while (genStartIndex < mrnaSize && rNAPositionIndexList.get(genStartIndex).start <= pos) {
                        genStartIndex++;
                    }
                    if (genStartIndex >= mrnaSize) {
                        genStartIndex = mrnaSize - 1;
                    }
                    //consider the max distance. and it will further categrize this according to 
                    //strand information of transcripts
                    //there are overlapped genes 
                    //can you provide some tail information to avoid exploring always from the 0 index
                    for (int searchIndex = 0; searchIndex <= genStartIndex; searchIndex++) {
                        RNASegment region = rNAPositionIndexList.get(searchIndex);
                        if (region.strand == '-') {
                            if (region.start - donwstreamDis <= pos && region.end + upstreamDis >= pos) {
                                rbi.addIndexes(region.mRNAIndex);
                            }
                        } else {
                            //default is +
                            if (region.start - upstreamDis <= pos && region.end + donwstreamDis >= pos) {
                                rbi.addIndexes(region.mRNAIndex);
                            }
                        }
                    }
                }
            }
        }
        hasNotSorted = false;
    }

    public void sortCNVMakeIndexonChromosomes() {
        geneChromPosMap.clear();

        int size = 0;
        List<RNASegment> rNAPositionIndexList = new ArrayList<RNASegment>();

        for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
            if (chromosomes[i] != null) {
                Collections.sort(chromosomes[i], new SeqSegmentComparator());
                size = chromosomes[i].size();
                rNAPositionIndexList.clear();
                if (chromosomemRNABoundaryIndexList[i] == null) {
                    chromosomemRNABoundaryIndexList[i] = new ArrayList<RNABoundaryIndex>();
                } else {
                    chromosomemRNABoundaryIndexList[i].clear();
                }
                for (int j = 0; j < size; j++) {
                    RefCNV cnv = (RefCNV) chromosomes[i].get(j);
                    geneChromPosMap.put(cnv.getDescription() + ":" + STAND_CHROM_NAMES[i] + ":" + cnv.getStart() + ":" + cnv.getEnd(), new int[]{i, j});
                    //default is +
                    RNASegment gns = new RNASegment(cnv.start - donwstreamDis, cnv.end + donwstreamDis, j, '+');
                    rNAPositionIndexList.add(gns);
                    chromosomemRNABoundaryIndexList[i].add(new RNABoundaryIndex(cnv.start - upstreamDis));
                    chromosomemRNABoundaryIndexList[i].add(new RNABoundaryIndex(cnv.end + donwstreamDis));
                }

                // the chromosomemRNABoundaryIndexList will be used to map a variant
                Collections.sort(rNAPositionIndexList, new GeneRNASegmentComparator());
                Collections.sort(chromosomemRNABoundaryIndexList[i], new RNABoundaryIndexComparator());
                int mrnaSize = rNAPositionIndexList.size();

                int boundSize = chromosomemRNABoundaryIndexList[i].size();
                for (int j = 0; j < boundSize; j++) {
                    RNABoundaryIndex rbi = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[i].get(j);

                    int pos = rbi.position;
                    int genStartIndex = binarySearchStartPos(pos, 0, mrnaSize - 1, rNAPositionIndexList);
                    if (genStartIndex < 0) {
                        genStartIndex = -genStartIndex - 1;
                    }
                    if (genStartIndex >= mrnaSize) {
                        genStartIndex = mrnaSize - 1;
                    }
                    while (genStartIndex < mrnaSize && rNAPositionIndexList.get(genStartIndex).start <= pos) {
                        genStartIndex++;
                    }
                    if (genStartIndex >= mrnaSize) {
                        genStartIndex = mrnaSize - 1;
                    }
                    //consider the max distance. and it will further categrize this according to 
                    //strand information of transcripts
                    //there are overlapped genes 
                    //can you provide some tail information to avoid exploring always from the 0 index
                    for (int searchIndex = 0; searchIndex <= genStartIndex; searchIndex++) {
                        RNASegment region = rNAPositionIndexList.get(searchIndex);
                        if (region.strand == '-') {
                            if (region.start - donwstreamDis <= pos && region.end + upstreamDis >= pos) {
                                rbi.addIndexes(region.mRNAIndex);
                            }
                        } else {
                            //default is +
                            if (region.start - upstreamDis <= pos && region.end + donwstreamDis >= pos) {
                                rbi.addIndexes(region.mRNAIndex);
                            }
                        }
                    }
                }
            }
        }
        hasNotSorted = false;
    }

    public void sortDupMakeIndexonChromosomes() {
        geneChromPosMap.clear();

        int size = 0;
        List<RNASegment> rNAPositionIndexList = new ArrayList<RNASegment>();

        for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
            if (chromosomes[i] != null) {
                Collections.sort(chromosomes[i], new SeqSegmentComparator());
                size = chromosomes[i].size();
                rNAPositionIndexList.clear();
                if (chromosomemRNABoundaryIndexList[i] == null) {
                    chromosomemRNABoundaryIndexList[i] = new ArrayList<RNABoundaryIndex>();
                } else {
                    chromosomemRNABoundaryIndexList[i].clear();
                }
                for (int j = 0; j < size; j++) {
                    RefDup cnv = (RefDup) chromosomes[i].get(j);
                    geneChromPosMap.put(cnv.getDescription() + ":" + STAND_CHROM_NAMES[i] + ":" + cnv.getStart() + ":" + cnv.getEnd(), new int[]{i, j});
                    //default is +
                    RNASegment gns = new RNASegment(cnv.start - donwstreamDis, cnv.end + donwstreamDis, j, '+');
                    rNAPositionIndexList.add(gns);
                    chromosomemRNABoundaryIndexList[i].add(new RNABoundaryIndex(cnv.start - upstreamDis));
                    chromosomemRNABoundaryIndexList[i].add(new RNABoundaryIndex(cnv.end + donwstreamDis));
                }

                // the chromosomemRNABoundaryIndexList will be used to map a variant
                Collections.sort(rNAPositionIndexList, new GeneRNASegmentComparator());
                Collections.sort(chromosomemRNABoundaryIndexList[i], new RNABoundaryIndexComparator());
                int mrnaSize = rNAPositionIndexList.size();

                int boundSize = chromosomemRNABoundaryIndexList[i].size();
                for (int j = 0; j < boundSize; j++) {
                    RNABoundaryIndex rbi = (RNABoundaryIndex) chromosomemRNABoundaryIndexList[i].get(j);

                    int pos = rbi.position;
                    int genStartIndex = binarySearchStartPos(pos, 0, mrnaSize - 1, rNAPositionIndexList);
                    if (genStartIndex < 0) {
                        genStartIndex = -genStartIndex - 1;
                    }
                    if (genStartIndex >= mrnaSize) {
                        genStartIndex = mrnaSize - 1;
                    }
                    while (genStartIndex < mrnaSize && rNAPositionIndexList.get(genStartIndex).start <= pos) {
                        genStartIndex++;
                    }
                    if (genStartIndex >= mrnaSize) {
                        genStartIndex = mrnaSize - 1;
                    }
                    //consider the max distance. and it will further categrize this according to 
                    //strand information of transcripts
                    //there are overlapped genes 
                    //can you provide some tail information to avoid exploring always from the 0 index
                    for (int searchIndex = 0; searchIndex <= genStartIndex; searchIndex++) {
                        RNASegment region = rNAPositionIndexList.get(searchIndex);
                        if (region.strand == '-') {
                            if (region.start - donwstreamDis <= pos && region.end + upstreamDis >= pos) {
                                rbi.addIndexes(region.mRNAIndex);
                            }
                        } else {
                            //default is +
                            if (region.start - upstreamDis <= pos && region.end + donwstreamDis >= pos) {
                                rbi.addIndexes(region.mRNAIndex);
                            }
                        }
                    }
                }
            }
        }
        hasNotSorted = false;
    }
}
