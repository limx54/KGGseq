/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import cern.colt.list.BooleanArrayList;
import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenLongObjectHashMap;
import cern.jet.stat.Probability;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.GZIPInputStream;
import javax.swing.tree.DefaultMutableTreeNode;
import org.apache.log4j.Logger;
import org.cobi.bayes.Bayes;
import org.cobi.kggseq.Constants;
import static org.cobi.kggseq.Constants.STAND_CHROM_NAMES;
import org.cobi.kggseq.GlobalManager;
import org.cobi.kggseq.entity.AnnotationSummarySet;
import org.cobi.kggseq.entity.Chromosome;
import org.cobi.kggseq.entity.CombOrders;
import org.cobi.kggseq.entity.FiltrationSummarySet;
import org.cobi.kggseq.entity.Gene;
import org.cobi.kggseq.entity.GeneFeature;
import org.cobi.kggseq.entity.Genome;
import org.cobi.kggseq.entity.Individual;
import org.cobi.kggseq.entity.IntSet;
import org.cobi.kggseq.entity.IntSetComparator1;
import org.cobi.kggseq.entity.PPIGraph;
import org.cobi.kggseq.entity.GeneSet;
import org.cobi.kggseq.entity.RNABoundaryIndex;
import org.cobi.kggseq.entity.RefCNV;
import org.cobi.kggseq.entity.RefDup;
import org.cobi.kggseq.entity.ReferenceGenome;
import org.cobi.kggseq.entity.RegressionParams;
import org.cobi.kggseq.entity.Variant;
import org.cobi.kggseq.entity.mRNA;
import org.cobi.randomforests.MyRandomForest;
import org.cobi.util.file.LocalFileFunc;
import org.cobi.util.net.NCBIRetriever;
import org.cobi.util.stat.ContingencyTable;
import org.cobi.util.stat.HardyWeinbergCalculation;
import org.cobi.util.stat.MultipleTestingMethod;
import org.cobi.util.text.BGZFInputStream;
import org.cobi.util.text.LocalFile;
import org.cobi.util.text.StringArrayStringComparator;
import org.cobi.util.text.Util;
import org.cobi.util.thread.Task;

/**
 *
 * @author mxli
 */
public class VariantAnnotator implements Constants {

    private static final Logger LOG = Logger.getLogger(VariantAnnotator.class);
    private int indexChrom = 0;
    private int indexPosition = 1;
    private int indexREF = 2;
    private int indexALT = 3;
    Map<String, Integer> chromNameIndexMap = new HashMap<String, Integer>();
    List<String> pubMedFilter;
    Map<String, Byte> geneFeatureMap = new HashMap<String, Byte>();

    public VariantAnnotator() {
        for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
            chromNameIndexMap.put(STAND_CHROM_NAMES[i], i);
            //variantPositionIndexMap[t] = new OpenIntIntHashMap();
        }
        for (int i = 0; i < VAR_FEATURE_NAMES.length; i++) {
            geneFeatureMap.put(VAR_FEATURE_NAMES[i], (byte) i);
        }
        pubMedFilter = new ArrayList<String>();
        pubMedFilter.add("gene");
        pubMedFilter.add("genes");
        pubMedFilter.add("mRNA");
        pubMedFilter.add("protein");
        pubMedFilter.add("proteins");
        pubMedFilter.add("transcription");
        pubMedFilter.add("transcript");
        pubMedFilter.add("transcripts");
        pubMedFilter.add("expressed");
        pubMedFilter.add("expression");
        pubMedFilter.add("expressions");
        pubMedFilter.add("locus");
        pubMedFilter.add("loci");
        pubMedFilter.add("SNP");
    }

    public Map<String, String[]> readGeneNames(String dbPath) throws Exception {
        List<String[]> geneItems = new ArrayList<String[]>();
        int[] indices = new int[3];
        /*
         * indices[0] = 1;//Approved\tSymbol indices[1] = 2;//Approved Name
         * indices[2] = 4;//Previous Symbols // indices[3] = 5; //Synonyms
         * indices[3] = 7; //Accession Numbers indices[4] = 8; //RefSeq IDs
         */

        indices[0] = 1;// Approved\tSymbol
        indices[1] = 7;// Approved Name
        indices[2] = 8;// Previous Symbols
        // indices[3] = 5; //Synonyms
        // indices[3] = 7; //Accession Numbers
        // indices[4] = 8; //RefSeq IDs

        LocalFile.retrieveData(dbPath, geneItems, indices, "\t");
        Map<String, String[]> geneSymMap = new HashMap<String, String[]>();
        List<String> tempList = new ArrayList<String>();

        for (String[] item : geneItems) {
            // ignore the full name of genes
            // tempList.add(item[1]);
            for (int i = 2; i < item.length; i++) {
                if (item[i].trim().isEmpty()) {
                    continue;
                }
                String[] cells = item[i].split(",");
                for (int t = 0; t < cells.length; t++) {
                    tempList.add(cells[t].trim());
                }
            }
            geneSymMap.put(item[0], tempList.toArray(new String[tempList.size()]));
            tempList.clear();
        }
        return geneSymMap;
    }

    //put this function in stack to speed up the parsing
    private void tokenize(String string, char delimiter, int maxIndex, String[] temp) {
        int wordCount = 0;
        int i = 0;
        int j = string.indexOf(delimiter);

        while (j >= 0) {
            temp[wordCount] = string.substring(i, j);
            if (wordCount >= maxIndex) {
                wordCount++;
                break;
            }
            wordCount++;
            i = j + 1;
            j = string.indexOf(delimiter, i);

        }
        if (wordCount <= maxIndex) {
            if (i < string.length()) {
                temp[wordCount++] = string.substring(i);
            }
        }
    }

    private int tokenize(String string, char delimiter, String[] temp) {
        int wordCount = 0;
        int i = 0;
        int j = string.indexOf(delimiter);

        while (j >= 0) {
            temp[wordCount++] = string.substring(i, j);
            i = j + 1;
            j = string.indexOf(delimiter, i);
        }

        if (i < string.length()) {
            temp[wordCount++] = string.substring(i);
        }
        return wordCount;
    }

    public void readExonicScoreNSFPNucleotideMerge(Chromosome chromosome, String resourcePath, String refGenomeVer,
            int[] scoreIndexes, int[] predicIndex, FiltrationSummarySet ass, boolean needProgressionIndicator) {
        int indexPos = 0;
        int indexref = 1;
        int indexalt = 2;
        int indexaaref = 6;
        int indexaaalt = 7;

        int maxColNum = indexPos;
        maxColNum = Math.max(maxColNum, indexPos);
        maxColNum = Math.max(maxColNum, indexref);
        maxColNum = Math.max(maxColNum, indexalt);
        int varFeatureNum = ass.getAvailableFeatureIndex();

        String currentLine = null;
        // String currChr = null;
        int currFilePostion = 0;

        char ref, alt;

        float[] scores = new float[scoreIndexes.length];
        String[] predicResults = new String[predicIndex.length];
        String missingLabel = ".";
        StringBuilder tmpBuffer = new StringBuilder();
        Variant[] vars = null;
        long lineCounter = 0;
        int varAssignedScore = 0;

        // for development of new method
        DoubleArrayList[] scoreLists = new DoubleArrayList[5];
        for (int i = 0; i < scoreLists.length; i++) {
            scoreLists[i] = new DoubleArrayList();
        }
        boolean fullMatch;

        // int varFeatureNum = genome.getVariantFeatureLabels().size();
        File rsFile = null;
        int unmatchedNum = 0;
        try {
            String chrName = chromosome.getName();
            if (chromosome == null || chromosome.variantList == null || chromosome.variantList.isEmpty()) {
                return;
            }

            rsFile = new File(resourcePath + chrName + ".gz");
            if (!rsFile.exists()) {
                LOG.warn(rsFile.getCanonicalPath() + " does not exist! Scores on this chromosome are ignored!");
                return;
            }
            // System.out.print(" Chromosome " +
            // Options.REF_CHROM_NAMES[t]);

            BufferedReader br = LocalFileFunc.getBufferedReader(rsFile.getCanonicalPath());

            lineCounter++;
            // skip to the head line
            currentLine = br.readLine();
            String[] cells = Util.tokenize(currentLine, '\t', maxColNum);
            String[] fullCells = Util.tokenize(currentLine, '\t');

            String[] tmpStrs;

            boolean hasNoScore;

            int varIndex = 0;
            // ensure all variants are sorted according to the
            // refStartPosition; otherwise itwil get stucked.
            Variant var = null;
            int varListSize = chromosome.variantList.size();
            int varPos = -1;
            boolean needNewRow = true;

            int varIndex1 = 0;
            int varPos1 = 0;
            var = chromosome.variantList.get(varIndex);
            varPos = var.refStartPosition;

            while (varIndex < varListSize) {
                // System.out.println(currentLine);
                if (needNewRow) {
                    // StringTokenizer st = new
                    // StringTokenizer(currentLine.trim(), "\t");
                    currentLine = br.readLine();
                    if (currentLine == null) {
                        break;
                    }
                    tokenize(currentLine, '\t', maxColNum, cells);
                    if (cells[indexPos].equals(missingLabel)) {
                        currFilePostion = -1;
                        continue;
                    } else {
                        currFilePostion = Util.parseInt(cells[indexPos]);
                    }
                    lineCounter++;
                    if (needProgressionIndicator && lineCounter % 50000 == 0) {
                        String prog = String.valueOf(lineCounter);
                        System.out.print(prog);
                        char[] backSpaces = new char[prog.length()];
                        Arrays.fill(backSpaces, '\b');
                        System.out.print(backSpaces);
                    }
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
                    var = chromosome.variantList.get(varIndex);
                    varPos = var.refStartPosition;
                    continue;
                }

                //ingore indel
                if (var.isIndel) {
                    varIndex++;
                    if (varIndex >= varListSize) {
                        break;
                    }
                    var = chromosome.variantList.get(varIndex);
                    varPos = var.refStartPosition;
                    continue;
                }

                // in the variant genome, there is only one unique variant
                varIndex1 = varIndex;

                //this alogrithm does not work for indels
                //it must be equal and there may be multiple positions equal
                ref = cells[indexref].charAt(0);
                alt = cells[indexalt].charAt(0);
                unmatchedNum = 0;
                //The list may have variants with the same coordinates 
                do {
                    var = chromosome.variantList.get(varIndex1);
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
                        tokenize(currentLine, '\t', fullCells);
                        Arrays.fill(scores, Float.NaN);
                        Arrays.fill(predicResults, null);

                        for (int iCol = 0; iCol < scoreIndexes.length; iCol++) {
                            if (fullCells[scoreIndexes[iCol]].equals(missingLabel)) {

                                scores[iCol] = Float.NaN;
                            } else if (!fullCells[scoreIndexes[iCol]].contains(";")) {
                                scores[iCol] = Util.parseFloat(fullCells[scoreIndexes[iCol]]);
                            } else {
                                // just use the first one
                                // System.out.println(tmpBuffer.toString());
                                tmpStrs = fullCells[scoreIndexes[iCol]].split(";");
                                hasNoScore = true;
                                for (String tmpStr : tmpStrs) {
                                    if (!tmpStr.equals(missingLabel)) {
                                        scores[iCol] = Util.parseFloat(tmpStr);
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
                            if (fullCells[predicIndex[iCol]].equals(missingLabel)) {

                                predicResults[iCol] = ".";
                            } else {

                                predicResults[iCol] = fullCells[predicIndex[iCol]];
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
                    var = chromosome.variantList.get(varIndex);
                    varPos = var.refStartPosition;
                }
                //otherwise only move file's index
                needNewRow = true;
            }
            br.close();
            ass.increaseCount(0, varAssignedScore);
            for (Variant var1 : chromosome.variantList) {
                if (var1.scores1 == null) {

                    var1.scores1 = new float[scores.length];
                    Arrays.fill(var1.scores1, Float.NaN);
                    for (int k = 0; k < predicResults.length; k++) {
                        var1.setFeatureValue(varFeatureNum + k, null);
                    }
                }
            }

        } catch (Exception e) {
            LOG.error("Exception at line \"" + currentLine + "\" of file " + rsFile.getName());
            e.printStackTrace();
        }
    }

    public void readExonicScoreNSFPNucleotideMergeMultiThread(Chromosome chromosome, String resourcePath, String refGenomeVer,
            int[] scoreIndexes, int[] predicIndex, FiltrationSummarySet ass, boolean needProgressionIndicator, int threadNum) {
        int indexPos = 0;
        int indexref = 1;
        int indexalt = 2;
        int indexaaref = 6;
        int indexaaalt = 7;

        int maxColNum = indexPos;
        maxColNum = Math.max(maxColNum, indexPos);
        maxColNum = Math.max(maxColNum, indexref);
        maxColNum = Math.max(maxColNum, indexalt);
        int varFeatureNum = ass.getAvailableFeatureIndex();

        String currentLine = null;
        // String currChr = null;
        int currFilePostion = 0;

        char ref, alt;

        float[] scores = new float[scoreIndexes.length];
        String[] predicResults = new String[predicIndex.length];
        String missingLabel = ".";
        StringBuilder tmpBuffer = new StringBuilder();
        Variant[] vars = null;
        long lineCounter = 0;
        int varAssignedScore = 0;

        // for development of new method
        DoubleArrayList[] scoreLists = new DoubleArrayList[5];
        for (int i = 0; i < scoreLists.length; i++) {
            scoreLists[i] = new DoubleArrayList();
        }
        boolean fullMatch;

        // int varFeatureNum = genome.getVariantFeatureLabels().size();
        File rsFile = null;
        int unmatchedNum = 0;
        try {
            String chrName = chromosome.getName();
            if (chromosome == null || chromosome.variantList == null || chromosome.variantList.isEmpty()) {
                return;
            }

            rsFile = new File(resourcePath + chrName + ".gz");
            if (!rsFile.exists()) {
                LOG.warn(rsFile.getCanonicalPath() + " does not exist! Scores on this chromosome are ignored!");
                return;
            }
            // System.out.print(" Chromosome " +
            // Options.REF_CHROM_NAMES[t]);
            int varSize = chromosome.variantList.size();
            int[] bounderies = partitionEvenBlock(varSize, threadNum);
            int actualThreadNum = bounderies.length - 1;
            String chromName = chromosome.getName();
            BGZFInputStream bfIndexes = new BGZFInputStream(rsFile.getCanonicalPath(), actualThreadNum, true);
            if (!bfIndexes.checkIndex()) {
                bfIndexes.adjustPos();
                bfIndexes.buildIndex(rsFile.getCanonicalPath(), -1, 0, true);
            }

            bfIndexes.readIndex(false, null);
            int indexCHROM = -1, indexPOS = 0;
            ExecutorService exec = Executors.newFixedThreadPool(actualThreadNum);
            final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
            int runningThread = 0;
            boolean hasHead = true;

            for (int f = 1; f < bounderies.length; f++) {
                long[] pos = bfIndexes.findIndex(chromosome.variantList.get(bounderies[f - 1]).refStartPosition, chromosome.variantList.get(bounderies[f] - 1).refStartPosition);
                if (pos[0] == pos[1]) {
                    for (int t = bounderies[f - 1]; t <= bounderies[f] - 1; t++) {
                        for (int s = 0; s < predicIndex.length; s++) {
                            chromosome.variantList.get(t).setFeatureValue(varFeatureNum + s, ".");
                        }
                        chromosome.variantList.get(t).scores1 = new float[scores.length];
                        Arrays.fill(chromosome.variantList.get(t).scores1, Float.NaN);
                    }
                    continue;
                }
                final BGZFInputStream bf = new BGZFInputStream(rsFile.getCanonicalPath(), 1, true);
                bf.adjustPos(pos[0], pos[1]);
                bf.creatSpider(pos[0] != 0);
                if (hasHead && pos[0] == 0) {
                    //skip the head line
                    bf.spider[0].readLine();
                }
                VarAnnotTask varTaks = new VarAnnotTask(bf.spider[0], chromosome.variantList.subList(bounderies[f - 1], bounderies[f]), indexCHROM, indexPOS, predicIndex, scoreIndexes, ass, 0) {
                    @Override
                    protected void fireTaskComplete() {
                        try {
                            bf.spider[0].closeInputStream();
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                };
                serv.submit(varTaks);
                runningThread++;
            }
            for (int s = 0; s < runningThread; s++) {
                Future<String> task = serv.take();
                String infor = task.get();
                // System.out.println(infor);
            }
            exec.shutdown();

            ass.increaseCount(0, varAssignedScore);
            for (Variant var1 : chromosome.variantList) {
                if (var1.scores1 == null) {
                    var1.scores1 = new float[scores.length];
                    Arrays.fill(var1.scores1, Float.NaN);
                    for (int k = 0; k < predicResults.length; k++) {
                        var1.setFeatureValue(varFeatureNum + k, null);
                    }
                }
            }

        } catch (Exception e) {
            LOG.error("Exception at line \"" + currentLine + "\" of file " + rsFile.getName());
            e.printStackTrace();
        }
    }

    //Important comments: Acually, the postitions of variants cannot be retreived by a positions because 1) a variant may be splitted into multiple ones according
    //to their alternative alleles and 2) the flexible annotation of Indels. This the the reason why we cannot use merged searching alogrithm to speeed up (Sadly!).
    //The good news is that there is a very fast hashmap OpenIntIntHashMap from Colt and I used it
    public int hardFilterByANNOVARefFormat(Chromosome chromosome, int chrID, AnnotationSummarySet ass, boolean needProgressionIndicator) {
        indexChrom = 0;
        indexPosition = 1;
        indexREF = 2;
        indexALT = 3;
        List<Variant> varList = chromosome.variantList;
        BufferedReader br = ass.getBr();
        StringBuilder newLine = ass.getLastLine();

        //No matched SNPs means null 
        //This is much faster than Set<Integer> 
        OpenIntIntHashMap excludeIndexes = new OpenIntIntHashMap();
        String currentLine = null;
        try {
            int maxColNum = indexChrom;
            maxColNum = Math.max(maxColNum, indexPosition);
            maxColNum = Math.max(maxColNum, indexREF);
            maxColNum = Math.max(maxColNum, indexALT);

            int lineCounter = 0;
            boolean incomplete = true;

            int filePosition = -1;
            StringBuilder tmpBuffer = new StringBuilder();
            String ref;
            String alt;
            String mafStr = null;

            String[] alts;
            //  String[] mafStrs;
            int[] varIndexes = null;
            boolean isDeleion = false;
            boolean isInsertion = false;

            char[] backSpaces = null;
            int delNum = 0;
            StringBuilder sb = new StringBuilder();

            String[] cells = null;
            if (newLine.length() == 0) {
                currentLine = br.readLine();
                if (currentLine == null) {
                    return 0;
                }
            } else {
                currentLine = newLine.toString();
                newLine.delete(0, newLine.length());
            }
            int fileChrID = 0;
            do {
                /*
                 if (currentLine.indexOf("BP") >= 0 || currentLine.indexOf("bp") >= 0) {
                 continue;
                 }
                 */
                lineCounter++;
                if (needProgressionIndicator && lineCounter % 50000 == 0) {
                    String prog = String.valueOf(lineCounter);
                    System.out.print(prog);
                    backSpaces = new char[prog.length()];
                    Arrays.fill(backSpaces, '\b');
                    System.out.print(backSpaces);
                }

                //StringTokenizer st = new StringTokenizer(currentLine.trim());
                cells = Util.tokenize(currentLine, '\t');
                if (cells.length < 2) {
                    cells = Util.tokenizeIngoreConsec(currentLine, ' ');
                }

                fileChrID = chromNameIndexMap.get(cells[indexChrom]);
                if (chrID < fileChrID) {
                    newLine.append(currentLine);
                    break;
                } else if (chrID > fileChrID) {
                    continue;
                }

                filePosition = Util.parseInt(cells[indexPosition]);
                //initialize varaibles
                incomplete = true;
                ref = cells[indexREF];
                alt = cells[indexALT];

                alts = alt.split(",");
                int tmpPos = 0;

                //once the variant is in db, it at least has a zero freq
                for (int s = 0; s < alts.length; s++) {
                    if (alts[s] == null || alts[s].isEmpty()) {
                        continue;
                    }

                    isDeleion = false;
                    isInsertion = false;
                    tmpPos = filePosition;
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

                    varIndexes = chromosome.lookupVariantIndexes(tmpPos);
                    if (varIndexes == null) {
                        continue;
                    }

                    // System.out.println(fileChr);

                    /*
                     if (chrID == null) {
                     //System.out.println("Unrecognized chromosome name: " + fileChr + " at line " + lineCounter + ": " + currentLine);
                     } else {
                     }
                     * 
                     */
                    for (int index : varIndexes) {
                        Variant var = varList.get(index);
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
                                            excludeIndexes.put(index, 0);
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
                                            excludeIndexes.put(index, 0);
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
                                    excludeIndexes.put(index, 0);
                                    break;
                                }
                            }
                        }
                    }
                }
            } while ((currentLine = br.readLine()) != null);

            int leftVarNum = varList.size();
            List<Variant> tmpVariantList = new ArrayList<Variant>();
            for (int j = 0; j < leftVarNum; j++) {
                if (!excludeIndexes.containsKey(j)) {
                    tmpVariantList.add(varList.get(j));
                }
            }
            ass.setAnnotNum(excludeIndexes.size() + ass.getAnnotNum());

            varList.clear();
            varList.addAll(tmpVariantList);
            tmpVariantList.clear();
            leftVarNum = varList.size();
            ass.setLeftNum(leftVarNum + ass.getLeftNum());

            ass.setTotalNum(lineCounter + ass.getTotalNum());
            // String info = leftVarNum + " variant(s) are left after hard filtering in database " + dbName + "!";
            //   LOG.info(info);
            if (needProgressionIndicator) {
                backSpaces = new char[7];
                Arrays.fill(backSpaces, '\b');
                System.out.print(backSpaces);
            }

            return leftVarNum;
        } catch (Exception ex) {
            if (currentLine != null) {
                System.err.println("Errors in a row: " + currentLine);
            }
            ex.printStackTrace();
        }
        return 0;
    }

    private void tokenizeIngoreConsec(String string, char delimiter, String[] temp) {
        int wordCount = 0;
        int i = 0;
        int j = string.indexOf(delimiter);

        while (j >= 0) {
            if (i < j) {
                temp[wordCount++] = string.substring(i, j);
            }
            i = j + 1;
            j = string.indexOf(delimiter, i);
        }

        if (i < string.length()) {
            temp[wordCount++] = string.substring(i);
        }
    }

    public void markByANNOVARefFormat(Chromosome chromosome, int chrID, AnnotationSummarySet ass, boolean needProgressionIndicator) {
        indexChrom = 0;
        indexPosition = 1;
        indexREF = 2;
        indexALT = 3;
        int indexMAF = 4;
        List<Variant> varList = chromosome.variantList;
        BufferedReader br = ass.getBr();
        StringBuilder newLine = ass.getLastLine();

        if (varList.isEmpty()) {
            return;
        }
        int feautreNum = ass.getAvailableFeatureIndex();

        //No matched SNPs means null 
        String missingVal = "N";
        for (Variant var : varList) {
            var.setFeatureValue(feautreNum, missingVal);
        }
        String currentLine = null;
        try {
            String varChrom = chromosome.getName();
            int maxColNum = indexChrom;
            maxColNum = Math.max(maxColNum, indexPosition);
            maxColNum = Math.max(maxColNum, indexREF);
            maxColNum = Math.max(maxColNum, indexALT);
            maxColNum = Math.max(maxColNum, indexMAF);

            int lineCounter = 0;

            int filePosition = -1;
            StringBuilder tmpBuffer = new StringBuilder();
            String ref;
            String alt;
            String mafStr = null;
            float maf;
            String[] alts;
            String[] mafStrs;
            int[] varIndex = null;
            boolean isDeleion = false;
            boolean isInsertion = false;

            /*
             File CHAIN_FILE = new File("./resources/hg19ToHg18.over.chain");
             LiftOver liftOver = new LiftOver(CHAIN_FILE);
             int failtedMapVarNum = 0;
             * 
             */
            if (needProgressionIndicator) {

            }
            char[] backSpaces = null;
            int delNum = 0;
            int existVarNum = 0;
            StringBuilder sb = new StringBuilder();
            boolean hitOnce = false;

            String[] cells = null;
            if (newLine.length() == 0) {
                currentLine = br.readLine();
                if (currentLine == null) {
                    return;
                }
            } else {
                currentLine = newLine.toString();
                newLine.delete(0, newLine.length());
            }
            int fileChrID;
            //get the column of a row
            cells = Util.tokenize(currentLine, '\t');
            boolean ingoreConsec = false;
            if (cells.length < 2) {
                cells = Util.tokenizeIngoreConsec(currentLine, ' ');
                ingoreConsec = true;
            }
            do {
                /*
                 if (currentLine.indexOf("BP") >= 0 || currentLine.indexOf("bp") >= 0) {
                 continue;
                 }*/
                lineCounter++;
                if (needProgressionIndicator && lineCounter % 50000 == 0) {
                    String prog = String.valueOf(lineCounter);
                    System.out.print(prog);
                    backSpaces = new char[prog.length()];
                    Arrays.fill(backSpaces, '\b');
                    System.out.print(backSpaces);
                }

                //StringTokenizer st = new StringTokenizer(currentLine.trim());
                if (ingoreConsec) {
                    tokenizeIngoreConsec(currentLine, ' ', cells);
                } else {
                    tokenize(currentLine, '\t', cells);
                }

                //initialize varaibles
                mafStrs = null;

                fileChrID = chromNameIndexMap.get(cells[indexChrom]);
                if (chrID < fileChrID) {
                    newLine.append(currentLine);
                    break;
                } else if (chrID > fileChrID) {
                    continue;
                }
                filePosition = Util.parseInt(cells[indexPosition]);
                ref = cells[indexREF];
                alt = cells[indexALT];
                if (cells.length > indexMAF) {
                    mafStr = cells[indexMAF];
                }

                alts = alt.split(",");
                if (mafStr != null) {
                    mafStrs = mafStr.split(",", -1);
                }

                //  System.err.println(currentLine);
                hitOnce = false;
                int tmpPos = 0;
                //once the variant is in db, it at least has a zero freq
                for (int s = 0; s < alts.length; s++) {
                    if (alts[s] == null || alts[s].isEmpty()) {
                        continue;
                    }
                    maf = Float.NaN;
                    if (mafStrs != null && s < mafStrs.length) {
                        //this missing score is denoted by .
                        if (mafStrs[s] != null && !mafStrs[s].isEmpty() && !mafStrs[s].equals(".")) {
                            maf = Util.parseFloat(mafStrs[s]);
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

                    varIndex = chromosome.lookupVariantIndexes(tmpPos);
                    if (varIndex == null) {
                        continue;
                    }

                    // System.out.println(fileChr);
                    for (int index : varIndex) {
                        Variant var = varList.get(index);
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
                }
                if (hitOnce) {
                    existVarNum++;
                }
            } while ((currentLine = br.readLine()) != null);

            ass.setLeftNum(existVarNum + ass.getLeftNum());
            ass.setAnnotNum(existVarNum + ass.getAnnotNum());
            ass.setTotalNum(lineCounter + ass.getTotalNum());
            if (needProgressionIndicator) {
                backSpaces = new char[7];
                Arrays.fill(backSpaces, '\b');
                System.out.print(backSpaces);
            }
        } catch (Exception ex) {
            if (currentLine != null) {
                System.err.println("Errors in a row: " + currentLine);
            }
            ex.printStackTrace();
        }
    }

    public void markByANNOVARefFormatThread(String filePath, Chromosome chromosome, AnnotationSummarySet ass, int threadNum) {
        int indexCHROM = 0, indexPOS = 1;
        List<Variant> varList = chromosome.variantList;

        if (varList.isEmpty()) {
            return;
        }
        int feautreNum = ass.getAvailableFeatureIndex();

        String currentLine = null;
        try {
            int lineCounter = 0;
            char[] backSpaces = null;
            int delNum = 0;
            int existVarNum = 0;
            StringBuilder sb = new StringBuilder();

            String chrName = chromosome.getName();
            if (chromosome == null || chromosome.variantList == null || chromosome.variantList.isEmpty()) {
                return;
            }

            File rsFile = new File(filePath);
            if (!rsFile.exists()) {
                LOG.warn(rsFile.getCanonicalPath() + " does not exist! Scores on this chromosome are ignored!");
                return;
            }
            // System.out.println(" Chromosome " + chromosome.getName() + " " + rsFile.getName());
            // Options.REF_CHROM_NAMES[t]);
            int varSize = chromosome.variantList.size();
            int[] bounderies = partitionEvenBlock(varSize, threadNum);
            int actualThreadNum = bounderies.length - 1;
            String chromName = chromosome.getName();
            BGZFInputStream bfIndexes = new BGZFInputStream(rsFile.getCanonicalPath(), actualThreadNum, true);
            if (!bfIndexes.checkIndex()) {
                bfIndexes.adjustPos();
                bfIndexes.buildIndex(rsFile.getCanonicalPath(), 0, 1, false);
            }
            bfIndexes.readIndex(false, chromName);

            //No matched SNPs means null 
            String missingVal = "N";
            for (Variant var : chromosome.variantList) {
                var.setFeatureValue(feautreNum, missingVal);
                var.hasBeenAcced = false;
            }

            ExecutorService exec = Executors.newFixedThreadPool(actualThreadNum);
            final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
            int runningThread = 0;
            boolean hasHead = false;
            for (int f = 1; f < bounderies.length; f++) {
                long[] pos = bfIndexes.findIndex(chromosome.variantList.get(bounderies[f - 1]).refStartPosition, chromosome.variantList.get(bounderies[f] - 1).refStartPosition);
                if (pos[0] == pos[1]) {
                    for (int t = bounderies[f - 1]; t <= bounderies[f] - 1; t++) {
                        chromosome.variantList.get(t).setFeatureValue(feautreNum, ".");
                    }
                    continue;
                }
                final BGZFInputStream bf = new BGZFInputStream(rsFile.getCanonicalPath(), 1, true);
                bf.adjustPos(pos[0], pos[1]);

                bf.creatSpider(pos[0] != 0);
                if (hasHead && pos[0] == 0) {
                    //skip the head line
                    bf.spider[0].readLine();
                }
                VarAnnotTask varTaks = new VarAnnotTask(bf.spider[0], indexCHROM, indexPOS, chromosome, ass, 1) {
                    @Override
                    protected void fireTaskComplete() {
                        try {
                            bf.spider[0].closeInputStream();
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                };
                serv.submit(varTaks);
                runningThread++;
            }

            for (int s = 0; s < runningThread; s++) {
                Future<String> task = serv.take();
                String infor = task.get();
                // System.out.println(infor);
            }
            exec.shutdown();

            ass.setLeftNum(existVarNum + ass.getLeftNum());
            ass.setAnnotNum(existVarNum + ass.getAnnotNum());
            ass.setTotalNum(lineCounter + ass.getTotalNum());
        } catch (Exception ex) {
            if (currentLine != null) {
                System.err.println("Errors in a row: " + currentLine);
            }
            ex.printStackTrace();
        }
    }

    public void markByVCFAlleleFormatThread(String filePath, int[] altFreqColIndexes, Chromosome chromosome, AnnotationSummarySet ass, int threadNum) {
        int indexCHROM = 0, indexPOS = 1;
        List<Variant> varList = chromosome.variantList;

        if (varList.isEmpty()) {
            return;
        }
        int feautreNum = ass.getAvailableFeatureIndex();

        String currentLine = null;
        try {
            int lineCounter = 0;

            int existVarNum = 0;

            if (chromosome == null || chromosome.variantList == null || chromosome.variantList.isEmpty()) {
                return;
            }

            File rsFile = new File(filePath);
            if (!rsFile.exists()) {
                LOG.warn(rsFile.getCanonicalPath() + " does not exist! Scores on this chromosome are ignored!");
                return;
            }
            // System.out.println(" Chromosome " + chromosome.getName() + " " + rsFile.getName());
            // Options.REF_CHROM_NAMES[t]);
            int varSize = chromosome.variantList.size();
            int[] bounderies = partitionEvenBlock(varSize, threadNum);
            int actualThreadNum = bounderies.length - 1;
            String chromName = chromosome.getName();
            BGZFInputStream bfIndexes = new BGZFInputStream(rsFile.getCanonicalPath(), actualThreadNum, true);
            if (!bfIndexes.checkIndex()) {
                bfIndexes.adjustPos();
                bfIndexes.buildIndex(rsFile.getCanonicalPath(), 0, 1, false);
            }
            bfIndexes.readIndex(false, chromName);

            //No matched SNPs means null 
            String missingVal = "N";
            for (Variant var : chromosome.variantList) {
                for (int i = 0; i < altFreqColIndexes.length; i++) {
                    var.setFeatureValue(feautreNum + i, missingVal);
                }
                var.hasBeenAcced = false;
            }

            ExecutorService exec = Executors.newFixedThreadPool(actualThreadNum);
            final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
            int runningThread = 0;
            boolean hasHead = false;
            for (int f = 1; f < bounderies.length; f++) {
                long[] pos = bfIndexes.findIndex(chromosome.variantList.get(bounderies[f - 1]).refStartPosition, chromosome.variantList.get(bounderies[f] - 1).refStartPosition);
                if (pos[0] == pos[1]) {
                    for (int t = bounderies[f - 1]; t <= bounderies[f] - 1; t++) {
                        chromosome.variantList.get(t).setFeatureValue(feautreNum, ".");
                    }
                    continue;
                }
                final BGZFInputStream bf = new BGZFInputStream(rsFile.getCanonicalPath(), 1, true);
                bf.adjustPos(pos[0], pos[1]);

                bf.creatSpider(pos[0] != 0);
                if (hasHead && pos[0] == 0) {
                    //skip the head line
                    bf.spider[0].readLine();
                }
                VarAnnotTask varTaks = new VarAnnotTask(bf.spider[0], indexCHROM, indexPOS, chromosome, ass, 3) {
                    @Override
                    protected void fireTaskComplete() {
                        try {
                            bf.spider[0].closeInputStream();
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                };
                varTaks.setAltFreqColIndexes(altFreqColIndexes);
                serv.submit(varTaks);
                runningThread++;
            }

            for (int s = 0; s < runningThread; s++) {
                Future<String> task = serv.take();
                String infor = task.get();
                // System.out.println(infor);
            }
            exec.shutdown();

            ass.setLeftNum(existVarNum + ass.getLeftNum());
            ass.setAnnotNum(existVarNum + ass.getAnnotNum());
            ass.setTotalNum(lineCounter + ass.getTotalNum());
        } catch (Exception ex) {
            if (currentLine != null) {
                System.err.println("Errors in a row: " + currentLine);
            }
            ex.printStackTrace();
        }
    }

    public void readNocodeScoreNucleotideMultiThread(String filePath, Chromosome chromosome, AnnotationSummarySet ass, int threadNum) {
        int indexCHROM = 0, indexPOS = 2;
        List<Variant> varList = chromosome.variantList;

        if (varList.isEmpty()) {
            return;
        }
        int feautreNum = ass.getAvailableFeatureIndex();

        String currentLine = null;
        try {
            int lineCounter = 0;
            int delNum = 0;
            int existVarNum = 0;
            StringBuilder sb = new StringBuilder();

            String chrName = chromosome.getName();
            if (chromosome == null || chromosome.variantList == null || chromosome.variantList.isEmpty()) {
                return;
            }

            File rsFile = new File(filePath);
            if (!rsFile.exists()) {
                LOG.warn(rsFile.getCanonicalPath() + " does not exist! Scores on this chromosome are ignored!");
                return;
            }
            // System.out.println(" Chromosome " + chromosome.getName() + " " + rsFile.getName());
            // Options.REF_CHROM_NAMES[t]);
            int varSize = chromosome.variantList.size();
            int[] bounderies = partitionEvenBlock(varSize, threadNum);
            int actualThreadNum = bounderies.length - 1;
            String chromName = chromosome.getName();
            BGZFInputStream bfIndexes = new BGZFInputStream(rsFile.getCanonicalPath(), actualThreadNum, true);
            if (!bfIndexes.checkIndex()) {
                bfIndexes.adjustPos();
                bfIndexes.buildIndex(rsFile.getCanonicalPath(), 0, 1, false);
            }
            bfIndexes.readIndex(false, chromName);

            //No matched SNPs means null 
            String missingVal = "N";
            for (Variant var : chromosome.variantList) {
                var.setFeatureValue(feautreNum, missingVal);
                var.hasBeenAcced = false;
            }
            int scoreNum = 10;

            ExecutorService exec = Executors.newFixedThreadPool(actualThreadNum);
            final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
            int runningThread = 0;
            boolean hasHead = false;
            for (int f = 1; f < bounderies.length; f++) {
                long[] pos = bfIndexes.findIndex(chromosome.variantList.get(bounderies[f - 1]).refStartPosition, chromosome.variantList.get(bounderies[f] - 1).refStartPosition);
                if (pos[0] == pos[1]) {
                    for (int t = bounderies[f - 1]; t <= bounderies[f] - 1; t++) {
                        chromosome.variantList.get(t).scores2 = new float[scoreNum];
                        Arrays.fill(chromosome.variantList.get(t).scores1, Float.NaN);
                    }
                    continue;
                }
                final BGZFInputStream bf = new BGZFInputStream(rsFile.getCanonicalPath(), 1, true);
                // System.out.println(chromosome.variantList.get(bounderies[f - 1]).refStartPosition+" : "+chromosome.variantList.get(bounderies[f] - 1).refStartPosition+" : "+pos[0] + " : " + pos[1]);
                bf.adjustPos(pos[0], pos[1]);

                bf.creatSpider(pos[0] != 0);
                if (hasHead && pos[0] == 0) {
                    //skip the head line
                    bf.spider[0].readLine();
                }
                VarAnnotTask varTaks = new VarAnnotTask(bf.spider[0], chromosome.variantList.subList(bounderies[f - 1], bounderies[f]), indexCHROM, indexPOS, ass, 2) {
                    @Override
                    protected void fireTaskComplete() {
                        try {
                            bf.spider[0].closeInputStream();
                        } catch (Exception ex) {
                            ex.printStackTrace();
                        }
                    }
                };
                serv.submit(varTaks);
                runningThread++;
            }

            for (int s = 0; s < runningThread; s++) {
                Future<String> task = serv.take();
                String infor = task.get();
                // System.out.println(infor);
            }
            exec.shutdown();

            ass.setLeftNum(existVarNum + ass.getLeftNum());
            ass.setAnnotNum(existVarNum + ass.getAnnotNum());
            ass.setTotalNum(lineCounter + ass.getTotalNum());
        } catch (Exception ex) {
            if (currentLine != null) {
                System.err.println("Errors in a row: " + currentLine);
            }
            ex.printStackTrace();
        }
    }

    public void riskPredictionRareDiseaseBest(Chromosome chromosome, List<CombOrders> combOrderList, boolean filterNonDisMut, List<String> names,
            FiltrationSummarySet dbNSFPMendelPred) throws Exception {

        int triedTime = 0;
        int maxTriedTime = 1;
        StringBuilder tmpStrB = new StringBuilder();

        // humvar predict
        double[] priors = new double[]{0.05, 0.01, 0};
        double[] mafBins = new double[]{0.01, 0.02, 0.03};
        double prior;

        List<Integer> paramIndexes = new ArrayList<Integer>();
        Set<Integer> dataIndexes = new HashSet<Integer>();
        int combListSize = combOrderList.size();
        List<Variant> tempVarList = new ArrayList<Variant>();
        int totalVarNum = 0;
        Set<String> geneSymbSet = new HashSet<String>();

        int counterDis = 0;
        int counterCon = 0;

        double sum = 0;
        RegressionParams tmpRP = null;
        RegressionParams bestRP = null;

        double tmpP, bestP;
        List<Integer> bestPParamIndexes = new ArrayList<Integer>();

        boolean isDeleteriousness = false;
        int featureNum = dbNSFPMendelPred.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            if (var.scores1 != null) {
                dataIndexes.clear();
                for (int j = 0; j < var.scores1.length; j++) {
                    if (!Float.isNaN(var.scores1[j])) {
                        dataIndexes.add(j);
                    }
                }

                bestRP = null;
                tmpRP = null;
                bestP = 0;
                isDeleteriousness = false;
                tmpStrB.delete(0, tmpStrB.length());
                triedTime = 0;

                for (int t = combListSize - 1; t >= 0; t--) {
                    CombOrders cmbOrder = combOrderList.get(t);
                    paramIndexes.clear();
                    if (dataIndexes.containsAll(cmbOrder.indexes)) {
                        paramIndexes.addAll(cmbOrder.indexes);

                        tmpRP = cmbOrder.rp;
                        sum = tmpRP.coef[0];
                        for (int j = 1; j < tmpRP.coef.length; j++) {
                            sum += (tmpRP.coef[j] * var.scores1[paramIndexes.get(j - 1)]);
                        }
                        // calculate the conditional probablity with MAF
                        // if (Float.isNaN(var.altAF) || var.altAF <= 0.01)
                        {
                            prior = ((1 - priors[0]) / priors[0]) * tmpRP.sampleCase2CtrRatio;
                            tmpP = 1 + prior * Math.exp(-sum);
                            tmpP = 1 / tmpP;
                        }
                        /*
                         * else if (var.altAF <= 0.03) { prior = ((1 -
                         * priors[1]) / priors[1]) *
                         * tmpRP.sampleCase2CtrRatio; tmpP = 1 + prior *
                         * Math.exp(-sum); tmpP = 1 / tmpP; } else { tmpP =
                         * 0; }
                         */

                        if (bestP <= tmpP) {
                            bestP = tmpP;
                            bestRP = tmpRP;
                            bestPParamIndexes.clear();
                            bestPParamIndexes.addAll(combOrderList.get(t).indexes);
                        }

                        if (tmpP >= tmpRP.optimalCutoff) {
                            var.setFeatureValue(featureNum, String.valueOf(tmpP));
                            var.setFeatureValue(featureNum + 1, "Y");
                            for (Integer ind : paramIndexes) {
                                tmpStrB.append(names.get(ind));
                                tmpStrB.append(';');
                            }
                            tmpStrB.deleteCharAt(tmpStrB.length() - 1);
                            tmpStrB.append(':');
                            tmpStrB.append(Util.doubleToString(tmpRP.optimalCutoff, 4));
                            tmpStrB.append(':');
                            tmpStrB.append(Util.doubleToString(tmpRP.truePositiveRate, 3));
                            tmpStrB.append(':');
                            tmpStrB.append(Util.doubleToString(tmpRP.trueNegativeRate, 3));
                            var.setFeatureValue(featureNum + 2, tmpStrB.toString());
                            counterDis += 1;
                            if (filterNonDisMut) {
                                tempVarList.add(var);
                            }
                            if (var.geneSymb != null) {
                                geneSymbSet.add(var.geneSymb);
                            }

                            isDeleteriousness = true;
                            // to reduce false positive only test once
                            break;
                        } else {
                            triedTime++;
                        }

                        if (triedTime >= maxTriedTime) {
                            // to reduce false positive only test once
                            break;
                        }
                    }
                }

                if (!isDeleteriousness) {
                    if (tmpRP == null) {
                        // keep varaint have no risk scores which may be
                        // safer
                        tempVarList.add(var);
                        var.setFeatureValue(featureNum, null);
                        var.setFeatureValue(featureNum + 1, null);
                        var.setFeatureValue(featureNum + 2, null);

                    } else {
                        // noly filter it for missense variants
                        if (var.smallestFeatureID != 6) {
                            tempVarList.add(var);
                        }
                        var.setFeatureValue(featureNum, String.valueOf(bestP));
                        var.setFeatureValue(featureNum + 1, "N");
                        if (!bestPParamIndexes.isEmpty()) {
                            for (Integer ind : bestPParamIndexes) {
                                tmpStrB.append(names.get(ind));
                                tmpStrB.append(';');
                            }
                            tmpStrB.deleteCharAt(tmpStrB.length() - 1);
                        }
                        tmpStrB.append(':');
                        tmpStrB.append(Util.doubleToString(bestRP.optimalCutoff, 4));
                        tmpStrB.append(':');
                        tmpStrB.append(Util.doubleToString(bestRP.truePositiveRate, 3));
                        tmpStrB.append(':');
                        tmpStrB.append(Util.doubleToString(bestRP.trueNegativeRate, 3));
                        var.setFeatureValue(featureNum + 2, tmpStrB.toString());
                        counterCon += 1;
                    }
                }
            } else {
                // keep varaint have no risk scores which may be safer
                tempVarList.add(var);
                var.setFeatureValue(featureNum, null);
                var.setFeatureValue(featureNum + 1, null);
                var.setFeatureValue(featureNum + 2, null);
            }
        }

        if (filterNonDisMut) {
            chromosome.variantList.clear();
            chromosome.variantList.addAll(tempVarList);
            totalVarNum += tempVarList.size();
            tempVarList.clear();

        }

        chromosome.buildVariantIndexMap();
        dbNSFPMendelPred.increaseCount(0, counterDis);
        dbNSFPMendelPred.increaseCount(1, geneSymbSet.size());
        dbNSFPMendelPred.increaseCount(2, counterCon);
        dbNSFPMendelPred.increaseCount(3, totalVarNum);
    }

    public void riskPredictionRareDiseaseFixParam(Chromosome chromosome, CombOrders cmbOrder, boolean filterNonDisMut, List<String> names, FiltrationSummarySet dbNSFPMendelPred) throws Exception {
        // humvar predict
        double[] priors = new double[]{0.05, 0.01, 0.0001};
        double[] mafBins = new double[]{0.01, 0.02, 0.03};
        double prior;

        List<Integer> paramIndexes = new ArrayList<Integer>();
        Set<Integer> dataIndexes = new HashSet<Integer>();

        List<Variant> tempVarList = new ArrayList<Variant>();
        int totalVarNum = 0;
        Set<String> geneSymbSet = new HashSet<String>();

        int counterDis = 0;
        int counterCon = 0;

        double sum = 0;
        RegressionParams tmpRP = null;
        RegressionParams bestRP = null;

        double tmpP, bestP;
        List<Integer> bestPParamIndexes = new ArrayList<Integer>();
        int featureNum = dbNSFPMendelPred.getAvailableFeatureIndex();
        boolean isDeleteriousness = false;
        StringBuilder tmpStrB = new StringBuilder();
        if (chromosome == null) {
            return;
        }
        for (Variant var : chromosome.variantList) {
            if (var.scores1 != null) {
                dataIndexes.clear();
                for (int j = 0; j < var.scores1.length; j++) {
                    if (!Float.isNaN(var.scores1[j])) {
                        dataIndexes.add(j);
                    }
                }

                tmpRP = null;
                bestRP = null;
                bestP = 0;
                isDeleteriousness = false;
                tmpStrB.delete(0, tmpStrB.length());

                paramIndexes.clear();
                if (dataIndexes.containsAll(cmbOrder.indexes)) {
                    paramIndexes.addAll(cmbOrder.indexes);
                    tmpRP = cmbOrder.rp;
                    sum = tmpRP.coef[0];
                    for (int j = 1; j < tmpRP.coef.length; j++) {
                        sum += (tmpRP.coef[j] * var.scores1[paramIndexes.get(j - 1)]);
                    }
                    // calculate the conditional probablity with MAF
                    // if (Float.isNaN(var.altAF) || var.altAF <= 0.01)
                    {
                        prior = ((1 - priors[0]) / priors[0]) * tmpRP.sampleCase2CtrRatio;
                        tmpP = 1 + prior * Math.exp(-sum);
                        tmpP = 1 / tmpP;
                    }
                    /*
                     * else if (var.altAF <= 0.03) { prior = ((1 -
                     * priors[1]) / priors[1]) *
                     * tmpRP.sampleCase2CtrRatio; tmpP = 1 + prior *
                     * Math.exp(-sum); tmpP = 1 / tmpP; } else { tmpP =
                     * 0; }
                     */

                    if (bestP <= tmpP) {
                        bestP = tmpP;
                        bestRP = tmpRP;
                        bestPParamIndexes.clear();
                        bestPParamIndexes.addAll(cmbOrder.indexes);
                    }

                    if (tmpP >= tmpRP.optimalCutoff) {
                        var.setFeatureValue(featureNum, String.valueOf(tmpP));
                        var.setFeatureValue(featureNum + 1, "Y");
                        for (Integer ind : paramIndexes) {
                            tmpStrB.append(names.get(ind));
                            tmpStrB.append(';');
                        }
                        tmpStrB.deleteCharAt(tmpStrB.length() - 1);
                        tmpStrB.append(':');
                        tmpStrB.append(Util.doubleToString(tmpRP.optimalCutoff, 4));
                        tmpStrB.append(':');
                        tmpStrB.append(Util.doubleToString(tmpRP.truePositiveRate, 3));
                        tmpStrB.append(':');
                        tmpStrB.append(Util.doubleToString(tmpRP.trueNegativeRate, 3));
                        var.setFeatureValue(featureNum + 2, tmpStrB.toString());
                        counterDis += 1;
                        if (filterNonDisMut) {
                            tempVarList.add(var);
                        }
                        if (var.geneSymb != null) {
                            geneSymbSet.add(var.geneSymb);
                        }
                        isDeleteriousness = true;
                    }
                }

                if (!isDeleteriousness) {
                    if (tmpRP == null) {
                        // keep varaint have no risk scores which may be
                        // safer
                        tempVarList.add(var);
                        var.setFeatureValue(featureNum, null);
                        var.setFeatureValue(featureNum + 1, null);
                        var.setFeatureValue(featureNum + 2, null);
                    } else {
                        // noly filter it for missense variants
                        if (var.smallestFeatureID != 6) {
                            tempVarList.add(var);
                        }
                        var.setFeatureValue(featureNum, String.valueOf(bestP));
                        var.setFeatureValue(featureNum + 1, "N");
                        if (!bestPParamIndexes.isEmpty()) {
                            for (Integer ind : bestPParamIndexes) {
                                tmpStrB.append(names.get(ind));
                                tmpStrB.append(';');
                            }
                            tmpStrB.deleteCharAt(tmpStrB.length() - 1);
                        }
                        tmpStrB.append(':');
                        tmpStrB.append(Util.doubleToString(bestRP.optimalCutoff, 4));
                        tmpStrB.append(':');
                        tmpStrB.append(Util.doubleToString(bestRP.truePositiveRate, 3));
                        tmpStrB.append(':');
                        tmpStrB.append(Util.doubleToString(bestRP.trueNegativeRate, 3));
                        var.setFeatureValue(featureNum + 2, tmpStrB.toString());
                        counterCon += 1;
                    }
                }
            } else {
                // keep varaint have no risk scores which may be safer
                tempVarList.add(var);
                var.setFeatureValue(featureNum, null);
                var.setFeatureValue(featureNum + 1, null);
                var.setFeatureValue(featureNum + 2, null);
            }
        }

        if (filterNonDisMut) {
            chromosome.variantList.clear();
            chromosome.variantList.addAll(tempVarList);
            totalVarNum += tempVarList.size();
            tempVarList.clear();
        }

        chromosome.buildVariantIndexMap();
        dbNSFPMendelPred.increaseCount(0, counterDis);
        dbNSFPMendelPred.increaseCount(1, geneSymbSet.size());
        dbNSFPMendelPred.increaseCount(2, counterCon);
        dbNSFPMendelPred.increaseCount(3, totalVarNum);

    }

    public void pseudogeneAnnotationGene(Chromosome chromosome, AnnotationSummarySet ass, String dbPath) throws Exception {
//        genome.addmRNAFeatureLabel("GeneDescription");
//        genome.addmRNAFeatureLabel("Pseudogenes");

        List<String[]> geneItems = new ArrayList<String[]>();
        int[] indices = new int[3];
        indices[0] = 1;
        indices[1] = 2;
        indices[2] = 3;
        LocalFile.retrieveData(dbPath, geneItems, indices, "\t");
        StringArrayStringComparator sacmp = new StringArrayStringComparator(0);
        Collections.sort(geneItems, sacmp);
        String lastGeneSymb = null;
        String[] lastGeneFeature = null;
        String[] keys = new String[3];
        int intCount = 0;

//        Chromosome[] chroms = genome.getChromosomes();
//        for (int t = 0; t < chroms.length; t++) {
        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        lastGeneFeature = null;
        for (mRNA mrna : chromosome.mRNAList) {
            if (mrna.geneSymb != null) {
                if (lastGeneSymb != null && mrna.geneSymb.equals(lastGeneSymb)) {
                    mrna.addFeatureValue(lastGeneFeature[0]);
                    mrna.addFeatureValue(lastGeneFeature[1]);
                } else {
                    keys[0] = mrna.geneSymb;
                    lastGeneFeature = searchRelevantGeneSymbols(geneItems, sacmp, keys, 0);
                    lastGeneSymb = mrna.geneSymb;
                    mrna.addFeatureValue(lastGeneFeature[0]);
                    mrna.addFeatureValue(lastGeneFeature[1]);
                }
                intCount++;
            } else {
                mrna.addFeatureValue(".");
                mrna.addFeatureValue(".");
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + intCount);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - intCount);
        //       }
    }

    public String[] searchRelevantGeneSymbols(List<String[]> geneItems, StringArrayStringComparator sacmp, String[] key, int keyIndex) throws Exception {
        String[] geneDescrip = new String[2];
        StringBuilder sb = new StringBuilder();
        int index = Collections.binarySearch(geneItems, key, sacmp);
        if (index < 0) {
            return geneDescrip;
        } else {
            int startIndex = index - 1;
            while (startIndex >= 0 && geneItems.get(startIndex)[keyIndex].startsWith(key[keyIndex])) {
                startIndex--;
            }
            startIndex++;
            int endIndex = index + 1;
            while (endIndex < geneItems.size() && geneItems.get(endIndex)[keyIndex].startsWith(key[keyIndex])) {
                endIndex++;
            }
            for (int i = startIndex; i < endIndex; i++) {
                if (geneItems.get(i)[1].contains("pseudogene")) {
                    sb.append(geneItems.get(i)[0]);
                    sb.append(", ");
                }
            }
            geneDescrip[0] = geneItems.get(index)[1] + " (" + geneItems.get(index)[2] + ")";
            if (sb.length() > 0) {
                geneDescrip[1] = sb.substring(0, sb.length() - 2);
            } else {
                geneDescrip[1] = ".";
            }
            return geneDescrip;
        }

    }

    public void pseudogeneAnnotationVar(Chromosome chromosome, AnnotationSummarySet ass, String dbPath) throws Exception {
//        genome.addVariantFeatureLabel("GeneDescription");
//        genome.addVariantFeatureLabel("Pseudogenes");

        List<String[]> geneItems = new ArrayList<String[]>();
        int[] indices = new int[3];
        indices[0] = 1;
        indices[1] = 2;
        indices[2] = 3;
        LocalFile.retrieveData(dbPath, geneItems, indices, "\t");//time-consuming
        StringArrayStringComparator sacmp = new StringArrayStringComparator(0);
        Collections.sort(geneItems, sacmp);
        String lastGeneSymb = null;
        String[] lastGeneFeature = null;
        String[] keys = new String[3];
        int intCount = 0;

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        lastGeneFeature = null;
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                if (lastGeneSymb != null && var.geneSymb.equals(lastGeneSymb)) {
                    var.setFeatureValue(ass.getAvailableFeatureIndex(), lastGeneFeature[0]);
                    var.setFeatureValue(ass.getAvailableFeatureIndex() + 1, lastGeneFeature[1]);
                } else {
                    keys[0] = var.geneSymb;
                    lastGeneFeature = searchRelevantGeneSymbols(geneItems, sacmp, keys, 0);
                    lastGeneSymb = var.geneSymb;
                    var.setFeatureValue(ass.getAvailableFeatureIndex(), lastGeneFeature[0]);
                    var.setFeatureValue(ass.getAvailableFeatureIndex() + 1, lastGeneFeature[1]);
                }
                intCount++;
            } else {
                var.setFeatureValue(ass.getAvailableFeatureIndex(), ".");
                var.setFeatureValue(ass.getAvailableFeatureIndex() + 1, ".");
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + intCount);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - intCount);
    }

    public void omimGeneAnnotationGene(Chromosome chromosome, AnnotationSummarySet ass, String dbPath) throws Exception {

        List<String[]> geneItems = new ArrayList<String[]>();
        int[] indices = new int[3];
        indices[0] = 0;
        indices[1] = 1;
        indices[2] = 2;
        LocalFile.retrieveData(dbPath, geneItems, indices, "[|]");
        Map<String, IntArrayList> geneRowMap = new HashMap<String, IntArrayList>();
        int size = geneItems.size();
        for (int i = 0; i < size; i++) {
            String[] item = geneItems.get(i);
            String[] genes = item[1].split(",");
            for (String gs : genes) {
                IntArrayList indexes = geneRowMap.get(gs);
                if (indexes == null) {
                    indexes = new IntArrayList();
                    indexes.add(i);
                    geneRowMap.put(gs, indexes);
                } else {
                    indexes.add(i);
                }
            }
        }

        String lastGeneSymb = null;
        String[] lastGeneFeature = null;

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        lastGeneFeature = null;
        int intCount = 0;
        for (mRNA mrna : chromosome.mRNAList) {
            if (mrna.geneSymb != null) {
                if (mrna.geneSymb != null) {
                    if (lastGeneSymb != null && mrna.geneSymb.equals(lastGeneSymb)) {
                        mrna.addFeatureValue(lastGeneFeature[0]);
                        mrna.addFeatureValue(lastGeneFeature[1]);
                    } else {
                        Arrays.fill(lastGeneFeature, null);
                        IntArrayList indexes = geneRowMap.get(mrna.geneSymb);
                        if (indexes != null) {
                            lastGeneFeature[0] = geneItems.get(indexes.getQuick(0))[0];
                            lastGeneFeature[1] = geneItems.get(indexes.getQuick(0))[2];
                            for (int t = 1; t < indexes.size(); t++) {
                                lastGeneFeature[0] += ("|" + geneItems.get(indexes.getQuick(t))[0]);
                                // lastGeneFeature[1] += ("|" +
                                // geneItems.get(indexes.getQuick(t))[2]);
                            }
                        }

                        lastGeneSymb = mrna.geneSymb;
                        mrna.addFeatureValue(lastGeneFeature[0]);
                        mrna.addFeatureValue(lastGeneFeature[1]);
                    }
                    intCount++;
                } else {
                    mrna.addFeatureValue(".");
                    mrna.addFeatureValue(".");
                }
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + intCount);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - intCount);
    }

    public void omimGeneAnnotationVar(Chromosome chromosome, AnnotationSummarySet ass, String dbPath) throws Exception {
        List<String[]> geneItems = new ArrayList<String[]>();
        int[] indices = new int[3];
        indices[0] = 0;
        indices[1] = 1;
        indices[2] = 2;
        LocalFile.retrieveData(dbPath, geneItems, indices, "[|]");
        Map<String, IntArrayList> geneRowMap = new HashMap<String, IntArrayList>();
        int size = geneItems.size();
        for (int i = 0; i < size; i++) {
            String[] item = geneItems.get(i);
            String[] genes = item[1].split(",");
            for (String gs : genes) {
                IntArrayList indexes = geneRowMap.get(gs);
                if (indexes == null) {
                    indexes = new IntArrayList();
                    indexes.add(i);
                    geneRowMap.put(gs, indexes);
                } else {
                    indexes.add(i);
                }
            }
        }

        String lastGeneSymb = null;
        String[] lastGeneFeature = null;
        lastGeneFeature = new String[2];
        int hitNum = 0;

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        Arrays.fill(lastGeneFeature, null);
        int avialbleIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                if (lastGeneSymb != null && var.geneSymb.equals(lastGeneSymb)) {
                    var.setFeatureValue(avialbleIndex, lastGeneFeature[0]);
                    var.setFeatureValue(avialbleIndex + 1, lastGeneFeature[1]);
                } else {
                    Arrays.fill(lastGeneFeature, null);
                    IntArrayList indexes = geneRowMap.get(var.geneSymb);
                    if (indexes != null) {
                        lastGeneFeature[0] = geneItems.get(indexes.getQuick(0))[0];
                        lastGeneFeature[1] = geneItems.get(indexes.getQuick(0))[2];
                        for (int t = 1; t < indexes.size(); t++) {
                            lastGeneFeature[0] += ("|" + geneItems.get(indexes.getQuick(t))[0]);
                            // lastGeneFeature[1] += ("|" +
                            // geneItems.get(indexes.getQuick(t))[2]);
                        }
                    }
                    lastGeneSymb = var.geneSymb;
                    var.setFeatureValue(avialbleIndex, lastGeneFeature[0]);
                    var.setFeatureValue(avialbleIndex + 1, lastGeneFeature[1]);
                    hitNum++;
                }
            } else {
                var.setFeatureValue(avialbleIndex, ".");
                var.setFeatureValue(avialbleIndex + 1, ".");
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - hitNum);
        // LOG.info(hitNum + " sequence variant(s) are highlighted by OMIM information!");
    }

    public void cosmicGeneAnnotationVar(Chromosome chromosome, AnnotationSummarySet ass, Map<String, Map<String, Integer>> cosmicGeneMut) throws Exception {
        String lastGeneSymb = null;
        String lastGeneFeature = ".";

        int hitNum = 0;

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;

        int avialbleIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                if (lastGeneSymb != null && var.geneSymb.equals(lastGeneSymb)) {
                    var.setFeatureValue(avialbleIndex, lastGeneFeature);
                } else {
                    lastGeneSymb = var.geneSymb;
                    Map<String, Integer> mutInfor = cosmicGeneMut.get(lastGeneSymb);
                    if (mutInfor != null) {
                        lastGeneFeature = mutInfor.toString();
                    } else {
                        lastGeneFeature = ".";
                    }
                    var.setFeatureValue(avialbleIndex, lastGeneFeature);
                    hitNum++;
                }
            } else {
                var.setFeatureValue(avialbleIndex, ".");
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - hitNum);
        // LOG.info(hitNum + " sequence variant(s) are highlighted by OMIM information!");
    }

    public void superDupAnnotation(Chromosome chromosome, AnnotationSummarySet ass, ReferenceGenome refGenome) throws Exception {
        int leftVariantNum = 0;
        int gainNum = 0;
        int lossNum = 0;
        int sampleSize = 0;
        StringBuilder sb = new StringBuilder();

        if (chromosome == null) {
            return;
        }
        int avialbleIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            List<RefDup> cnvList = refGenome.getDupFeature(chromosome.getName(), var, new RNABoundaryIndex(0));
            if (cnvList == null || cnvList.isEmpty()) {
                var.setFeatureValue(avialbleIndex, ".");
                continue;
            }

            for (RefDup refCNC : cnvList) {
                sb.append(refCNC.getJcKScore());
                sb.append(", ");
            }
            leftVariantNum++;
            var.setFeatureValue(avialbleIndex, sb.substring(0, sb.length() - 2));

            sb.delete(0, sb.length());
        }

        ass.setAnnotNum(ass.getAnnotNum() + leftVariantNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - leftVariantNum);

    }

    public void cnvAnnotation(Chromosome chromosome, AnnotationSummarySet ass, ReferenceGenome refGenome) throws Exception {
        int leftVariantNum = 0;
        int gainNum = 0;
        int lossNum = 0;
        int sampleSize = 0;

        StringBuilder sb = new StringBuilder();

        if (chromosome == null) {
            return;
        }
        int avialbleIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            List<RefCNV> cnvList = refGenome.getCNVFeature(chromosome.getName(), var, new RNABoundaryIndex(0));
            if (cnvList == null || cnvList.isEmpty()) {
                var.setFeatureValue(avialbleIndex, ".");
                var.setFeatureValue(avialbleIndex + 1, ".");
                var.setFeatureValue(avialbleIndex + 2, ".");
                var.setFeatureValue(avialbleIndex + 3, ".");
                continue;
            }
            gainNum = 0;
            lossNum = 0;
            sampleSize = 0;
            for (RefCNV refCNC : cnvList) {
                sampleSize += refCNC.getSampleSize();
                gainNum += refCNC.getObservedGains();
                lossNum += refCNC.getObservedLosses();
                sb.append(refCNC.getDescription());
                sb.append(", ");
            }
            leftVariantNum++;
            var.setFeatureValue(avialbleIndex, sb.substring(0, sb.length() - 2));
            var.setFeatureValue(avialbleIndex + 1, String.valueOf(sampleSize));
            var.setFeatureValue(avialbleIndex + 2, String.valueOf(lossNum));
            var.setFeatureValue(avialbleIndex + 3, String.valueOf(gainNum));

            sb.delete(0, sb.length());
        }

        ass.setAnnotNum(ass.getAnnotNum() + leftVariantNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - leftVariantNum);
    }

    public void canidateGeneExploreGene(Chromosome chromosome, AnnotationSummarySet ass, Set<String> candiGeneSet) {
        int hitNum = 0;
        if (chromosome == null) {
            return;
        }
        for (mRNA mrna : chromosome.mRNAList) {
            if (mrna.geneSymb != null) {
                if (candiGeneSet.contains(mrna.geneSymb)) {
                    mrna.addFeatureValue("Y");
                    hitNum++;
                } else {
                    mrna.addFeatureValue("N");
                }
            } else {
                mrna.addFeatureValue(".");
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - hitNum);
    }

    public void canidateGeneExploreVar(Chromosome chromosome, AnnotationSummarySet ass, Set<String> candiGeneSet) {
        int hitNum = 0;
        if (chromosome == null) {
            return;
        }
        int afI = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                if (candiGeneSet.contains(var.geneSymb)) {
                    var.setFeatureValue(afI, "Y");
                    hitNum++;
                } else {
                    var.setFeatureValue(afI, "N");
                }
            } else {
                var.setFeatureValue(afI, ".");
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - hitNum);
    }

    public void canidateGenePPIExploreGene(Chromosome chromosome, AnnotationSummarySet ass, Set<String> candiGeneSet, PPIGraph ppiTree, int ppiDepth) throws Exception {
        String lastGeneSymb = null;
        String lastGenePPIFeature = null;

        StringBuilder sb = new StringBuilder();
        DefaultMutableTreeNode proteinTree = null;

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        lastGenePPIFeature = "";
        int hitNum = 0;
        for (mRNA mrna : chromosome.mRNAList) {
            if (mrna.getSymbol() != null) {
                if (lastGeneSymb != null && mrna.getSymbol().equals(lastGeneSymb)) {//Maybe isoforms?
                    mrna.addFeatureValue(lastGenePPIFeature);
                    if (lastGenePPIFeature.length() > 0) {
                        hitNum++;
                    }
                } else {
                    sb.delete(0, sb.length());

                    proteinTree = new DefaultMutableTreeNode(mrna.getSymbol());
                    ppiTree.constructPPITree(proteinTree, 0, ppiDepth);
                    ppiTree.trimPPITree(proteinTree, candiGeneSet);
                    if (proteinTree.getChildCount() > 0) {
                        ppiTree.convertTree2Text(proteinTree, sb);
                        // System.out.println(sb);
                    }
                    if (sb.length() > 0) {
                        lastGenePPIFeature = sb.substring(0, sb.length() - 1);
                        hitNum++;
                    } else {
                        lastGenePPIFeature = "";
                    }
                    lastGeneSymb = mrna.getSymbol();
                    mrna.addFeatureValue(lastGenePPIFeature);
                }
            } else {
                mrna.addFeatureValue(null);
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - hitNum);
    }

    public void canidateGenePPIExploreVar(Chromosome chromosome, AnnotationSummarySet ass, Set<String> candiGeneSet, PPIGraph ppiTree, int ppiDepth) throws Exception {
        String lastGeneSymb = null;
        String lastGenePPIFeature = null;

        StringBuilder sb = new StringBuilder();
        DefaultMutableTreeNode proteinTree = null;
        int hitNum = 0;

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        lastGenePPIFeature = "";

        int avialbleIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                if (lastGeneSymb != null && var.geneSymb.equals(lastGeneSymb)) {
                    var.setFeatureValue(avialbleIndex, lastGenePPIFeature);
                    if (lastGenePPIFeature.length() > 0) {
                        hitNum++;
                    }
                } else {
                    sb.delete(0, sb.length());

                    proteinTree = new DefaultMutableTreeNode(var.geneSymb);
                    ppiTree.constructPPITree(proteinTree, 0, ppiDepth);
                    ppiTree.trimPPITree(proteinTree, candiGeneSet);
                    if (proteinTree.getChildCount() > 0) {
                        ppiTree.convertTree2Text(proteinTree, sb);
                        // System.out.println(sb);
                    }
                    if (sb.length() > 0) {
                        lastGenePPIFeature = sb.substring(0, sb.length() - 1);
                        hitNum++;
                    } else {
                        lastGenePPIFeature = "";
                    }
                    lastGeneSymb = var.geneSymb;
                    var.setFeatureValue(avialbleIndex, lastGenePPIFeature);
                }
            } else {
                var.setFeatureValue(avialbleIndex, null);
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - hitNum);
        //LOG.info(hitNum + " sequence variant(s) are highlighted by PPI information!");
    }

    public void canidateGenePathwayExploreGene(Chromosome chromosome, AnnotationSummarySet ass, Set<String> candiGeneSet, Map<String, GeneSet> mappedPathes) throws Exception {

        String lastGeneSymb = null;
        // protein-coding gene \t19061
        int popuSize = 19061;
        int subPopuSize = 0;
        int sampleSize = 0;
        int subSampleSize = 0;

        String lastGenePathwayFeature = null;
        StringBuilder sb = new StringBuilder();
        List<PathwayP> pathwayPList = new ArrayList<PathwayP>();

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        lastGenePathwayFeature = "";
        int hitNum = 0;
        for (mRNA mrna : chromosome.mRNAList) {
            if (mrna.geneSymb != null) {
                if (lastGeneSymb != null && mrna.geneSymb.equals(lastGeneSymb)) {
                    mrna.addFeatureValue(lastGenePathwayFeature);
                    if (lastGenePathwayFeature.length() > 0) {
                        hitNum++;
                    }
                } else {
                    // find pahtway information
                    pathwayPList.clear();
                    for (Map.Entry<String, GeneSet> mPathway : mappedPathes.entrySet()) {
                        HashSet<String> pathwayGenes = mPathway.getValue().getGeneSymbols();
                        if (pathwayGenes.contains(mrna.geneSymb)) {
                            sb.delete(0, sb.length());
                            sb.append(mPathway.getKey());
                            sb.append('#');
                            sb.append(pathwayGenes.size());
                            sb.append(":(");
                            sb.append(mrna.geneSymb);
                            if (!candiGeneSet.contains(mrna.geneSymb)) {
                                subSampleSize = 0;
                            } else {
                                subSampleSize = 1;
                            }

                            Iterator<String> itSeed = candiGeneSet.iterator();
                            while (itSeed.hasNext()) {
                                String gene1 = itSeed.next();
                                if (pathwayGenes.contains(gene1)) {
                                    sb.append(',');
                                    sb.append(gene1);
                                }
                            }
                            sb.append(") ");
                            subPopuSize = pathwayGenes.size();
                            sampleSize = candiGeneSet.size();
                            if (!candiGeneSet.contains(mrna.geneSymb)) {
                                sampleSize++;
                            }

                            double p = MultipleTestingMethod.hypergeometricEnrichmentTest(popuSize, subPopuSize, sampleSize, subSampleSize);
                            pathwayPList.add(new PathwayP(p, sb.toString()));
                        }
                    }

                    sb.delete(0, sb.length());
                    if (!pathwayPList.isEmpty()) {
                        Collections.sort(pathwayPList, new PathwayPComparator());
                        for (PathwayP pp : pathwayPList) {
                            sb.append(pp.pathway);
                            sb.append(" p=").append(Util.formatPValue(pp.p));
                            sb.append("; ");
                        }
                    }

                    if (sb.length() > 0) {
                        lastGenePathwayFeature = sb.substring(0, sb.length() - 2);
                        hitNum++;
                    } else {
                        lastGenePathwayFeature = "";
                    }
                    lastGeneSymb = mrna.geneSymb;
                    mrna.addFeatureValue(lastGenePathwayFeature);
                }
            } else {
                mrna.addFeatureValue(null);
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - hitNum);
    }

    public void canidateGenePathwayExploreVar(Chromosome chromosome, AnnotationSummarySet ass, Set<String> candiGeneSet, Map<String, GeneSet> mappedPathes) throws Exception {
        // protein-coding gene \t19061
        int popuSize = 19061;

        int subPopuSize = 0;
        int sampleSize = 0;
        int subSampleSize = 0;
        String lastGeneSymb = null;
        String lastGenePathwayFeature = null;
        StringBuilder sb = new StringBuilder();
        List<PathwayP> pathwayPList = new ArrayList<PathwayP>();
        int hitNum = 0;

        if (chromosome == null) {
            return;
        }

        lastGeneSymb = null;
        lastGenePathwayFeature = "";

        int avialbleIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                if (lastGeneSymb != null && var.geneSymb.equals(lastGeneSymb)) {
                    var.setFeatureValue(avialbleIndex, lastGenePathwayFeature);
                    if (lastGenePathwayFeature.length() > 0) {
                        hitNum++;
                    }
                } else {
                    // find pahtway information
                    pathwayPList.clear();

                    for (Map.Entry<String, GeneSet> mPathway : mappedPathes.entrySet()) {
                        HashSet<String> pathwayGenes = mPathway.getValue().getGeneSymbols();
                        if (pathwayGenes.contains(var.geneSymb)) {
                            sb.delete(0, sb.length());
                            sb.append(mPathway.getKey());
                            sb.append('#');
                            sb.append(pathwayGenes.size());
                            sb.append(":(");
                            sb.append(var.geneSymb);

                            if (!candiGeneSet.contains(var.geneSymb)) {
                                subSampleSize = 0;
                            } else {
                                subSampleSize = 1;
                            }
                            Iterator<String> itSeed = candiGeneSet.iterator();
                            while (itSeed.hasNext()) {
                                String gene1 = itSeed.next();
                                if (pathwayGenes.contains(gene1)) {
                                    sb.append(',');
                                    sb.append(gene1);
                                    subSampleSize++;
                                }
                            }
                            sb.append(") ");
                            subPopuSize = pathwayGenes.size();
                            sampleSize = candiGeneSet.size();
                            if (!candiGeneSet.contains(var.geneSymb)) {
                                sampleSize++;
                            }
                            double p = MultipleTestingMethod.hypergeometricEnrichmentTest(popuSize, subPopuSize, sampleSize, subSampleSize);
                            pathwayPList.add(new PathwayP(p, sb.toString()));
                        }
                    }
                    sb.delete(0, sb.length());
                    if (!pathwayPList.isEmpty()) {
                        Collections.sort(pathwayPList, new PathwayPComparator());
                        for (PathwayP pp : pathwayPList) {
                            sb.append(pp.pathway);
                            sb.append(" p=").append(Util.formatPValue(pp.p));
                            sb.append("; ");
                        }
                    }

                    if (sb.length() > 0) {
                        lastGenePathwayFeature = sb.substring(0, sb.length() - 2);
                        hitNum++;
                    } else {
                        lastGenePathwayFeature = "";
                    }
                    lastGeneSymb = var.geneSymb;
                    var.setFeatureValue(avialbleIndex, lastGenePathwayFeature);
                }
            } else {
                var.setFeatureValue(avialbleIndex, null);
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - hitNum);
        //LOG.info(hitNum + " sequence variant(s) are highlighted by GeneSet/Pathway information!");
    }

    public void ibdRegionExplore(Chromosome chromosome, AnnotationSummarySet ass, List<String[]> regionItems) {

        List<Variant> tmpList = new ArrayList<Variant>();
//        int hitNum = 0;
        if (chromosome == null) {
            return;
        }

        int avialbleIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            String[] ibdGrams = searchIbdRegions(regionItems, chromosome.getName(), var.refStartPosition);
            if (ibdGrams != null) {
                var.setFeatureValue(avialbleIndex, ibdGrams[0]);
                var.setFeatureValue(avialbleIndex + 1, ibdGrams[1]);
                tmpList.add(var);
            } else {
                var.setFeatureValue(avialbleIndex, null);
                var.setFeatureValue(avialbleIndex + 1, null);
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpList.size());
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpList);
        tmpList.clear();
//        hitNum += chromosome.variantList.size();
        chromosome.buildVariantIndexMap();
//        genome.buildVariantIndexMapOnChromosomes();
//        genome.setVarNum(hitNum);
//        LOG.info(hitNum + " variant(s) are left after filtered by IBD Region filtering!");
    }

    public void exploreLongIBSRegion(Chromosome chromosome, AnnotationSummarySet ass, List<Variant> allVariants, int minIBS, boolean isPhased, List<Individual> subjectList, int[] pedEncodeGytIDMap, int[] caeSetID) {
        int caseNum = caeSetID.length;
        if (caseNum <= 1) {
            String info = "The number of patients is only " + caseNum + ".  Identical by state (IBS) checking (--ibs-check) function doese not work!";
            LOG.info(info);
            return;
        }
        IntSet tmpIntSet = new IntSet(0, 0, (byte) 0);

        String endMak = "";
        String startMak = "";
        StringBuilder info = new StringBuilder();
        int[] alleles1 = null;
        int[] alleles2 = null;
        IntSetComparator1 intsetCmp = new IntSetComparator1();
        OpenIntIntHashMap sharedAlleles = new OpenIntIntHashMap();
        boolean hasShared = true;
        int index1 = 0, index2 = 0;
        int begIndex = 0;
        int subID = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();
        int hitNum = 0;

        if (chromosome == null) {
            return;
        }
        endMak = "";
        startMak = "";
        OpenIntIntHashMap posIDMap = new OpenIntIntHashMap();

        int base = 0;
        int alleleNum = 0;
        int varNum = allVariants.size();
        BooleanArrayList isValid = new BooleanArrayList(varNum);
        boolean[] bits = new boolean[32];
        int startIndex;
        int subNum = pedEncodeGytIDMap.length;
        for (Variant var : allVariants) {
            sharedAlleles.clear();
            begIndex = 0;
            alleleNum = var.getAltAlleles().length + 1;
            if (isPhased) {
                base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
            } else {
                base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
            }

            hasShared = true;
            do {
                subID = caeSetID[begIndex];
                subID = pedEncodeGytIDMap[subID];
                if (subID < 0) {
                    begIndex++;
                    continue;
                }
                if (var.compressedGtyLabel >= 0) {
                    if (var.compressedGtyLabel == 0) {
                        Arrays.fill(bits, 0, base, false);
                    } else {
                        startIndex = subID;
                        for (int i = 0; i < base; i++) {
                            if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = false;
                            } else if (startIndex < var.compressedGty[0]) {
                                bits[i] = false;
                            } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = true;
                            } else if (startIndex == var.compressedGty[0]) {
                                bits[i] = true;
                            } else {
                                bits[i] = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                            }
                            startIndex += subNum;
                        }
                    }

                    if (isPhased) {
                        alleles1 = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                    } else {
                        alleles1 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                    }
                } else if (isPhased) {
                    alleles1 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                } else {
                    alleles1 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                }

                if (alleles1 == null) {
                    begIndex++;
                } else {
                    break;
                }
            } while (begIndex < caseNum);

            if (alleles1 == null) {
                continue;
            }
            sharedAlleles.put(alleles1[0], 0);
            sharedAlleles.put(alleles1[1], 0);

            for (int j = begIndex + 1; j < caseNum; j++) {
                subID = caeSetID[j];
                subID = pedEncodeGytIDMap[subID];
                if (subID < 0) {
                    continue;
                }
                if (var.compressedGtyLabel >= 0) {
                    if (var.compressedGtyLabel == 0) {
                        Arrays.fill(bits, 0, base, false);
                    } else {
                        startIndex = subID;
                        for (int i = 0; i < base; i++) {
                            if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = false;
                            } else if (startIndex < var.compressedGty[0]) {
                                bits[i] = false;
                            } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = true;
                            } else if (startIndex == var.compressedGty[0]) {
                                bits[i] = true;
                            } else {
                                bits[i] = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                            }
                            startIndex += subNum;
                        }
                    }

                    if (isPhased) {
                        alleles2 = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                    } else {
                        alleles2 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                    }
                } else if (isPhased) {
                    alleles2 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                } else {
                    alleles2 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                }

                if (alleles2 == null) {
                    continue;
                }
                if (!sharedAlleles.containsKey(alleles2[0]) && !sharedAlleles.containsKey(alleles2[1])) {
                    hasShared = false;
                    break;
                }
            }
            isValid.add(hasShared);
            posIDMap.put(var.refStartPosition, isValid.size() - 1);
            //System.out.println(var.refStartPosition + "\t" + hasShared);
        }

        index1 = 0;
        index2 = 0;
        tmpVarList.clear();
        int regionLen = 0;
        int availbleIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            tmpIntSet.index1 = var.refStartPosition;
            index1 = posIDMap.get(tmpIntSet.index1);

            if (index1 <= 0) {
                var.setFeatureValue(availbleIndex, null);
                var.setFeatureValue(availbleIndex + 1, null);
                continue;
            }
            index1--;
            index2 = index1 + 1;
            while (index2 < varNum && isValid.getQuick(index2)) {
                index2++;
            }
            while (index1 >= 0 && isValid.getQuick(index1)) {
                index1--;
            }
            if (index2 >= varNum) {
                index2 = varNum - 1;
                endMak = "EndAtTail";
            }
            if (index1 < 0) {
                index1 = 0;
                startMak = "StartAtTail";
            }
            info.delete(0, info.length());
            info.append("chr").append(STAND_CHROM_NAMES[chromosome.getId()]).append(":").append(allVariants.get(index1).refStartPosition).append("-").append(allVariants.get(index2).refStartPosition)
                    .append("#Var:").append(index2 - index1 + 1);

            if (startMak.length() > 0) {
                info.append(startMak);
            }
            if (endMak.length() > 0) {
                info.append(endMak);
            }
            var.setFeatureValue(availbleIndex, info.toString());
            regionLen = allVariants.get(index2).refStartPosition - allVariants.get(index1).refStartPosition + 1;
            var.setFeatureValue(availbleIndex + 1, String.valueOf(regionLen));
            if (regionLen >= minIBS) {
                tmpVarList.add(var);
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVarList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVarList.size());

        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarList);
        hitNum += tmpVarList.size();
        tmpVarList.clear();

        chromosome.buildVariantIndexMap();
        //testingGenome.buildVariantIndexMapOnChromosomes();
        //testingGenome.setVarNum(hitNum);
        //LOG.info(hitNum + " variant(s) are left after filtered by IBS filtering!");
    }

    public void exploreLongHomozygosityRegion(Chromosome chromosome, AnnotationSummarySet ass, List<Variant> allVariants, int minIBS, boolean isPhased, List<Individual> subjectList, int[] pedEncodeGytIDMap, int[] caeSetID, int[] controlSetID) {
        int caseNum = caeSetID.length;
        if (caseNum <= 0) {
            String info = "There is no patients. Identical by state (IBS) checking (--ibs-check) function doese not work!";
            LOG.info(info);
            return;
        }
        int controlNum = controlSetID.length;
        IntSet tmpIntSet = new IntSet(0, 0, (byte) 0);

        String endMak = "";
        String startMak = "";
        StringBuilder info = new StringBuilder();

        int[] alleles1 = null;
        int[] alleles2 = null;
        IntSetComparator1 intsetCmp = new IntSetComparator1();

        boolean hasShared = true;
        int index1 = 0, index2 = 0;
        int begIndex = 0;
        int hitNum = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();
        int subID = 0;

        if (chromosome == null) {
            return;
        }
        endMak = "";
        startMak = "";
        OpenIntIntHashMap posIDMap = new OpenIntIntHashMap();

        int base = 0;
        int alleleNum = 0;
        int varNum = allVariants.size();
        BooleanArrayList isValid = new BooleanArrayList(varNum);
        int startIndex;
        boolean[] bits = new boolean[32];
        int subNum = pedEncodeGytIDMap.length;
        for (Variant var : allVariants) {
            begIndex = 0;
            if (begIndex >= caseNum) {
                continue;
            }
            alleleNum = var.getAltAlleles().length + 1;
            if (isPhased) {
                base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
            } else {
                base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
            }
            hasShared = true;
            do {
                subID = caeSetID[begIndex];
                subID = pedEncodeGytIDMap[subID];
                if (subID < 0) {
                    begIndex++;
                    continue;
                }
                if (var.compressedGtyLabel >= 0) {
                    if (var.compressedGtyLabel == 0) {
                        Arrays.fill(bits, 0, base, false);
                    } else {
                        startIndex = subID;
                        for (int i = 0; i < base; i++) {
                            if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = false;
                            } else if (startIndex < var.compressedGty[0]) {
                                bits[i] = false;
                            } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = true;
                            } else if (startIndex == var.compressedGty[0]) {
                                bits[i] = true;
                            } else {
                                bits[i] = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                            }
                            startIndex += subNum;
                        }
                    }
                    if (isPhased) {
                        alleles1 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                    } else {
                        alleles1 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                    }
                } else if (isPhased) {
                    alleles1 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                } else {
                    alleles1 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                }

                if (alleles1 == null) {
                    begIndex++;
                } else {
                    break;
                }
            } while (begIndex < caseNum);

            if (alleles1 == null) {
                isValid.add(true);
                posIDMap.put(var.refStartPosition, isValid.size() - 1);
                continue;
            }
            if (alleles1[0] != alleles1[1]) {
                isValid.add(false);
                posIDMap.put(var.refStartPosition, isValid.size() - 1);
                continue;
            }

            for (int j = begIndex + 1; j < caseNum; j++) {
                subID = caeSetID[j];
                subID = pedEncodeGytIDMap[subID];
                if (subID < 0) {
                    continue;
                }
                if (var.compressedGtyLabel >= 0) {
                    if (var.compressedGtyLabel == 0) {
                        Arrays.fill(bits, 0, base, false);
                    } else {
                        startIndex = subID;
                        for (int i = 0; i < base; i++) {
                            if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = false;
                            } else if (startIndex < var.compressedGty[0]) {
                                bits[i] = false;
                            } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = true;
                            } else if (startIndex == var.compressedGty[0]) {
                                bits[i] = true;
                            } else {
                                bits[i] = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                            }
                            startIndex += subNum;
                        }
                    }
                    if (isPhased) {
                        alleles2 = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                    } else {
                        alleles2 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                    }
                } else if (isPhased) {
                    alleles2 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                } else {
                    alleles2 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                }

                if (alleles2 == null) {
                    continue;
                }
                if (alleles2[0] != alleles2[1] || alleles2[0] != alleles1[0]) {
                    hasShared = false;
                    break;
                }
            }
            for (int j = 0; j < controlNum; j++) {
                subID = controlSetID[j];
                subID = pedEncodeGytIDMap[subID];
                if (subID < 0) {
                    continue;
                }
                if (var.compressedGtyLabel >= 0) {
                    if (var.compressedGtyLabel == 0) {
                        Arrays.fill(bits, 0, base, false);
                    } else {
                        startIndex = subID;
                        for (int i = 0; i < base; i++) {
                            if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = false;
                            } else if (startIndex < var.compressedGty[0]) {
                                bits[i] = false;
                            } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                                bits[i] = true;
                            } else if (startIndex == var.compressedGty[0]) {
                                bits[i] = true;
                            } else {
                                bits[i] = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                            }
                            startIndex += subNum;
                        }
                    }
                    if (isPhased) {
                        alleles2 = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                    } else {
                        // alleles2 = subjectList.get(subID).markerGtySetArray.getUnphasedGtyBool(pair.index2, pair.len);
                        alleles2 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                    }
                } else if (isPhased) {
                    alleles2 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                } else {
                    // alleles2 = subjectList.get(subID).markerGtySetArray.getUnphasedGtyBool(pair.index2, pair.len);
                    alleles2 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, subNum);
                }

                if (alleles2 == null) {
                    continue;
                }
                if (alleles2[0] == alleles2[1] && alleles2[0] == alleles1[0]) {
                    hasShared = false;
                    break;
                }
            }
            isValid.add(hasShared);
            posIDMap.put(var.refStartPosition, isValid.size() - 1);
        }
        index1 = 0;
        index2 = 0;
        tmpVarList.clear();
        int regionLen = 0;
        int availableIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            tmpIntSet.index1 = var.refStartPosition;
            index1 = posIDMap.get(tmpIntSet.index1);
            if (index1 <= 0) {
                var.setFeatureValue(availableIndex, null);
                var.setFeatureValue(availableIndex + 1, null);
                continue;
            }
            index1--;
            index2 = index1 + 1;
            while (index2 < varNum && isValid.getQuick(index2)) {
                index2++;
            }
            while (index1 >= 0 && isValid.getQuick(index1)) {
                index1--;
            }
            if (index2 >= varNum) {
                index2 = varNum - 1;
                endMak = "EndAtTail";
            }
            if (index1 < 0) {
                index1 = 0;
                startMak = "StartAtTail";
            }
            info.delete(0, info.length());
            info.append("chr").append(STAND_CHROM_NAMES[chromosome.getId()]).append(":").append(allVariants.get(index1).refStartPosition).append("-").append(allVariants.get(index2).refStartPosition)
                    .append("#Var:").append(index2 - index1 + 1);

            if (startMak.length() > 0) {
                info.append(startMak);
            }
            if (endMak.length() > 0) {
                info.append(endMak);
            }
            var.setFeatureValue(availableIndex, info.toString());
            regionLen = allVariants.get(index2).refStartPosition - allVariants.get(index1).refStartPosition + 1;
            var.setFeatureValue(availableIndex + 1, String.valueOf(regionLen));
            if (regionLen >= minIBS) {
                tmpVarList.add(var);
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVarList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVarList.size());

        hitNum += tmpVarList.size();
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarList);
        tmpVarList.clear();

        chromosome.buildVariantIndexMap();

    }

    // model 0: allelic; 1: dominant; 2: recessive; 3: genotypic
    public void assocTestVar(Chromosome chromosome, AnnotationSummarySet ass, DoubleArrayList[] pValueList, Genome genome) throws Exception {

        // AffectedRefHomGtyNum\tAffectedHetGtyNum\tAffectedAltHomGtyNum\tUnaffectedRefHomGtyNum\tUnaffectedHetGtyNum\tUnaffectedAltHomGtyNum
//        Chromosome[] chroms = genome.getChromosomes();
        int[] caseGtys = new int[3];
        int[] controlGtys = new int[3];

        boolean hasLess5 = false;

        long[][] countsInt = new long[2][3];
        double[] oddsValues = new double[3];
        double[] pValues = new double[4];

        if (chromosome == null) {
            return;
        }
        int avialbleIndex = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            Arrays.fill(caseGtys, 0);
            caseGtys[0] = var.getAffectedRefHomGtyNum();
            caseGtys[1] = var.getAffectedHetGtyNum();
            caseGtys[2] = var.getAffectedAltHomGtyNum();
            Arrays.fill(controlGtys, 0);
            controlGtys[0] = var.getUnaffectedRefHomGtyNum();
            controlGtys[1] = var.getUnaffectedHetGtyNum();
            controlGtys[2] = var.getUnaffectedAltHomGtyNum();
            /*
             * % wild risk % ___________ % case | A | B |Rs1 % |_____|_____|
             * % control | C | D |Rs2 % |_____|_____|____ % Cs1 Cs2 N
             */
            hasLess5 = false;
            // model 0: allelic;
            for (int k = 0; k < 2; k++) {
                Arrays.fill(countsInt[k], 0);
            }
            countsInt[0][0] = caseGtys[0] * 2 + caseGtys[1];
            countsInt[0][1] = caseGtys[2] * 2 + caseGtys[1];
            countsInt[1][0] = controlGtys[0] * 2 + controlGtys[1];
            countsInt[1][1] = controlGtys[2] * 2 + controlGtys[1];

            for (int k = 0; k < 2; k++) {
                for (int j = 0; j < 2; j++) {
                    if (countsInt[k][j] < 0) {
                        countsInt[k][j] = 0;
                    }
                    if (countsInt[k][j] < 5) {
                        hasLess5 = true;
                    }
                }
            }

            oddsValues[0] = ((double) countsInt[0][1] / (countsInt[1][1] == 0 ? 0.5f : countsInt[1][1]))
                    * ((double) countsInt[1][0] / (countsInt[0][0] == 0 ? 0.5f : countsInt[0][0]));
            if (hasLess5) {
                pValues[0] = ContingencyTable.fisherExact22(countsInt, 2, 2, 2);
            } else {
                pValues[0] = ContingencyTable.pearsonChiSquared22(countsInt);
                pValues[0] = Probability.chiSquareComplemented(1, pValues[0]);
            }
            if (!Double.isNaN(pValues[0])) {
                pValueList[0].add(pValues[0]);
            }

            hasLess5 = false;
            // 1: dominant;
            countsInt[0][0] = caseGtys[0];
            countsInt[0][1] = caseGtys[2] + caseGtys[1];
            countsInt[1][0] = controlGtys[0];
            countsInt[1][1] = controlGtys[2] + controlGtys[1];

            for (int k = 0; k < 2; k++) {
                for (int j = 0; j < 2; j++) {
                    if (countsInt[k][j] < 0) {
                        countsInt[k][j] = 0;
                    }
                    if (countsInt[k][j] < 5) {
                        hasLess5 = true;
                    }
                }
            }

            oddsValues[1] = ((double) countsInt[0][1] / (countsInt[1][1] == 0 ? 0.5f : countsInt[1][1]))
                    * ((double) countsInt[1][0] / (countsInt[0][0] == 0 ? 0.5f : countsInt[0][0]));
            if (hasLess5) {
                pValues[1] = ContingencyTable.fisherExact22(countsInt, 2, 2, 2);
            } else {
                pValues[1] = ContingencyTable.pearsonChiSquared22(countsInt);
                pValues[1] = Probability.chiSquareComplemented(1, pValues[1]);
            }

            if (!Double.isNaN(pValues[1])) {
                pValueList[1].add(pValues[1]);
            }
            hasLess5 = false;
            // 2: recessive; 3: genotypic
            countsInt[0][0] = caseGtys[0] + caseGtys[1];
            countsInt[0][1] = caseGtys[2];
            countsInt[1][0] = controlGtys[0] + controlGtys[1];
            countsInt[1][1] = controlGtys[2];

            for (int k = 0; k < 2; k++) {
                for (int j = 0; j < 2; j++) {
                    if (countsInt[k][j] < 0) {
                        countsInt[k][j] = 0;
                    }
                    if (countsInt[k][j] < 5) {
                        hasLess5 = true;
                    }
                }
            }

            oddsValues[2] = ((double) countsInt[0][1] / (countsInt[1][1] == 0 ? 0.5f : countsInt[1][1]))
                    * ((double) countsInt[1][0] / (countsInt[0][0] == 0 ? 0.5f : countsInt[0][0]));
            if (hasLess5) {
                pValues[2] = ContingencyTable.fisherExact22(countsInt, 2, 2, 2);
            } else {
                pValues[2] = ContingencyTable.pearsonChiSquared22(countsInt);
                pValues[2] = Probability.chiSquareComplemented(1, pValues[2]);
            }
            if (!Double.isNaN(pValues[2])) {
                pValueList[2].add(pValues[2]);
            }

            countsInt[0][0] = caseGtys[0];
            countsInt[0][1] = caseGtys[1];
            countsInt[0][2] = caseGtys[2];
            countsInt[1][0] = controlGtys[0];
            countsInt[1][1] = controlGtys[1];
            countsInt[1][2] = controlGtys[2];

            pValues[3] = ContingencyTable.chiSquareTest(countsInt);
            pValues[3] = Probability.chiSquareComplemented(2, pValues[3]);
            /*
             for (int k = 0; k < 2; k++) {
             for (int j = 0; j < 3; j++) {
             if (countsInt[k][j] < 0) {
             countsInt[k][j] = 0;
             }
             if (countsInt[k][j] < 5) {
             hasLess5 = true;
             }
             }
             }
             if (hasLess5) {
             pValues[3] = ContingencyTable.fisherExact22(countsInt, 2, 3, 2);
             } else {
             pValues[3] = ContingencyTable.chiSquareTest(countsInt);
             pValues[3] = Probability.chiSquareComplemented(2, pValues[3]);
             }*/

            for (int t = 0; t < 3; t++) {
                if (Double.isNaN(pValues[t])) {
                    var.setFeatureValue(avialbleIndex + t * 2, ".");
                } else {
                    var.setFeatureValue(avialbleIndex + t * 2, String.valueOf(pValues[t]));
                }

                if (Double.isNaN(oddsValues[t])) {
                    var.setFeatureValue(avialbleIndex + t * 2 + 1, ".");
                } else {
                    var.setFeatureValue(avialbleIndex + t * 2 + 1, String.valueOf(oddsValues[t]));
                }
            }
            if (!Double.isNaN(pValues[3])) {
                pValueList[3].add(pValues[3]);
                var.setFeatureValue(avialbleIndex + 6, String.valueOf(pValues[3]));
            } else {
                var.setFeatureValue(avialbleIndex + 6, ".");
            }
        }
    }

    public void pubMedIDGeneExploreGene(Chromosome chromosome, AnnotationSummarySet ass, List<String> pubmedMeshList, Map<String, String[]> geneNamesMap, Map<String, String> genePubMedID) throws Exception {
        NCBIRetriever ncbiRetriever = new NCBIRetriever();
        String lastGeneSymb = null;
        String lastGeneFeature = null;
        int account = 1;
        List<String> geneNames = new ArrayList<String>();
        String[] names;

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        lastGeneFeature = null;
        for (mRNA mrna : chromosome.mRNAList) {
            if (mrna.geneSymb != null) {
                if (lastGeneSymb != null && mrna.geneSymb.equals(lastGeneSymb)) {
                    mrna.addFeatureValue(lastGeneFeature);
                } else {
                    lastGeneFeature = genePubMedID.get(mrna.geneSymb);
                    if (lastGeneFeature == null) {
                        geneNames.clear();

                        names = geneNamesMap.get(mrna.geneSymb);
                        if (names != null && names.length > 0) {
                            geneNames.addAll(Arrays.asList(names));
                        }

                        geneNames.add(mrna.geneSymb);
                        LOG.info(account + ": Searching NCBI PubMed for " + pubmedMeshList.toString() + " and " + geneNames.toString());
                        while ((lastGeneFeature = ncbiRetriever.pubMedIDESearch(pubmedMeshList, geneNames, pubMedFilter)) == null) {
                            // System.out.print("reconnecting...");
                        }//Will the pubMedFilter be dropped?
                        // System.out.println();
                        account++;
                    }

                    lastGeneSymb = mrna.geneSymb;
                    mrna.addFeatureValue(lastGeneFeature);
                }
            } else {
                mrna.addFeatureValue(null);
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + account);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - account);
    }

    public void pubMedIDIdeogramExploreVar(Chromosome chromosome, AnnotationSummarySet ass, List<String> pubmedMeshList, List<String[]> ideogramItems, String refGenomeVersion, boolean ignoreNonPathogenic, int driverPredicIndex) throws Exception {
        if (!pubmedMeshList.isEmpty() && pubmedMeshList.get(0).equals("ANY")) {
            return;
        }

        NCBIRetriever ncbiRetriever = new NCBIRetriever();
        int account = 1;

        List<String> ideoGrams = null;
        Map<String, String> historyResultMap = new HashMap<String, String>();
        StringBuilder result = new StringBuilder();
        String predicType;

        if (chromosome == null) {
            return;
        }
        int featureNum = ass.getAvailableFeatureIndex();
        int hitNum = 0;
        for (Variant var : chromosome.variantList) {

            //only did it for missens variant
            if (var.smallestFeatureID == 6) {
                if (ignoreNonPathogenic && driverPredicIndex >= 0) {
                    predicType = var.getFeatureValues()[driverPredicIndex];
                    if (predicType != null && predicType.equals("N")) {
                        var.setFeatureValue(featureNum, null);
                        continue;
                    }
                }
            }

            result.delete(0, result.length());
            ideoGrams = searchIdeogramRegions(ideogramItems, chromosome.getName(), var.refStartPosition);
            for (String reg : ideoGrams) {
                // ignore regions which are two broad
                /*
                 * if (reg.indexOf(".") < 0) { continue; }
                 */
                String pubIDs = historyResultMap.get(reg);
                if (pubIDs == null) {
                    if (ignoreNonPathogenic && driverPredicIndex >= 0) {
                        predicType = var.getFeatureValues()[driverPredicIndex];
                        if (predicType != null && predicType.equals("N")) {
                            var.setFeatureValue(featureNum, null);
                            continue;
                        }
                    }

                    LOG.info(account + ": Searching NCBI PubMed for " + pubmedMeshList.toString() + " and " + reg);
                    while ((pubIDs = ncbiRetriever.pubMedIDESearch(pubmedMeshList, reg)) == null) {
                        // System.out.print("reconnecting...");
                    }

                    // System.out.println(pubIDs);
                    if (pubIDs == null) {
                        pubIDs = "";
                    }
                    // System.out.println();
                    historyResultMap.put(reg, pubIDs);
                    account++;
                }
                if (pubIDs.trim().length() > 0) {
                    result.append(reg).append(":[").append(pubIDs).append("] ");
                }
            }
            if (result != null && result.length() > 0) {
                hitNum++;
            }
            var.setFeatureValue(featureNum, result.toString());
        }

        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - hitNum);
    }

    public void pubMedIDGeneExploreVar(Chromosome chromosome, AnnotationSummarySet ass, List<String> pubmedMeshList, boolean ignoreNonPathogenic, Map<String, String[]> geneNamesMap, Map<String, String> genePubMedID, int driverPredicIndex) throws Exception {
        NCBIRetriever ncbiRetriever = new NCBIRetriever();

        String lastGeneSymb = null;
        String lastGeneFeature = null;
        int account = 1;
        String predicType = null;
        List<String> geneNames = new ArrayList<String>();

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        lastGeneFeature = null;
        int featureNum = ass.getAvailableFeatureIndex();
        int hitNum = 0;
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                //only did it for missens variant
                if (var.smallestFeatureID == 6) {
                    if (ignoreNonPathogenic && driverPredicIndex >= 0) {
                        predicType = var.getFeatureValues()[driverPredicIndex];
                        if (predicType != null && predicType.equals("N")) {
                            var.setFeatureValue(featureNum, null);
                            continue;
                        }
                    }
                }
                if (lastGeneSymb != null && var.geneSymb.equals(lastGeneSymb)) {
                    var.setFeatureValue(featureNum, lastGeneFeature);
                } else {
                    lastGeneFeature = genePubMedID.get(var.geneSymb);
                    if (lastGeneFeature == null) {
                        geneNames.clear();
                        // unforturnately after I add the alais, there are
                        // too many false postive hits
                        // so I finally withraw this function
                        String[] names = geneNamesMap.get(var.geneSymb);
                        if (names != null && names.length > 0) {
                            geneNames.addAll(Arrays.asList(names));
                        }

                        geneNames.add(var.geneSymb);

                        LOG.info(account + ": Searching NCBI PubMed for " + pubmedMeshList.toString() + " and " + geneNames.toString());
                        while ((lastGeneFeature = ncbiRetriever.pubMedIDESearch(pubmedMeshList, geneNames, pubMedFilter)) == null) {
                            // System.out.print("reconnecting...");
                        }//Will pubMiedFilter be dropped?
                        account++;
                        // System.out.println();
                    }

                    lastGeneSymb = var.geneSymb;
                    var.setFeatureValue(featureNum, lastGeneFeature);
                }
                if (lastGeneFeature != null && lastGeneFeature.length() > 0) {
                    hitNum++;
                }
            } else {
                var.setFeatureValue(featureNum, null);
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - hitNum);
    }

    // Read ChrPosRs File
    public void readChrPosRs(OpenIntIntHashMap[] altMapID, File resourceFileRSID)
            throws FileNotFoundException, IOException {

        BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(resourceFileRSID))));
        String line = null;
        while ((line = reader.readLine()) != null) {
            String[] array = line.split("\t");
            int id = Integer.parseInt(array[0]);
            int chr = 0;
            if (array[1].equals("AltOnly") || array[1].equals("NotOn") || array[1].equals("PAR") || array[1].equals("Multi") || array[1].equals("") || array[2].equals("")) {
                continue;
            }
            if (array[1].matches("\\d+")) {
                chr = Integer.parseInt(array[1]);
            } else if (array[1].equals("X")) {
                chr = 23;
            } else if (array[1].equals("Y")) {
                chr = 24;
            } else if (array[1].equals("XY")) {
                chr = 25;
            } else if (array[1].equals("MT")) {
                chr = 26;
            } else if (array[1].equals("Un")) {
                chr = 27;
            }
            int index = chr - 1;
            int pos = Integer.parseInt(array[2]);
            // pos in dbsnp SNPChrPosOnRef is 0-based, So we changed it into 1-based and save in hashMap
            pos++;
            altMapID[index].put(pos, id);
        }
        reader.close();
    }

    public void addRSID(Chromosome chr, AnnotationSummarySet ass, OpenIntIntHashMap altMapID) {
        int total = 0;
        int annotNum = 0;
        String chrNam = chr.getName();
        for (Variant var : chr.variantList) {
            int pos = var.refStartPosition;
            total++;
            if (altMapID.containsKey(pos)) {
                annotNum++;
                var.label = "rs" + altMapID.get(pos);
            } else {
                var.label = chrNam + ":" + pos;
            }
        }
        ass.setAnnotNum(ass.getAnnotNum() + annotNum);
        ass.setTotalNum(ass.getTotalNum() + total);
        ass.setLeftNum(ass.getLeftNum() + total - annotNum);
    }

    public void assignSelectedVar2Genes(List<Variant> variantList, String chrName, Map<String, List<Variant>> geneVars) {
        String geneFeatureAnnot;
        StringBuilder sb = new StringBuilder();
        int index;
        Set<String> geneSet = new HashSet<String>();
        Map<String, Set<String>> unifyingSet = new HashMap<String, Set<String>>();

        StringBuilder id = new StringBuilder();
        for (Variant var : variantList) {
            String strGene = var.geneSymb;
            if (strGene == null) {
                continue;
            }
            geneFeatureAnnot = var.getRefGeneAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            geneFeatureAnnot = var.getgEncodeAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            geneFeatureAnnot = var.getKnownGeneAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            geneFeatureAnnot = var.getEnsemblGeneAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }

            if (sb.length() > 0) {
                String[] items = Util.tokenize(sb.substring(0, sb.length() - 1), ';');
                for (String item : items) {
                    index = item.indexOf(':');
                    if (index > 0) {
                        geneSet.add(item.substring(0, index));
                    }
                }
                sb.delete(0, sb.length());
            }
            if (var.label == null || var.label.equals(".")) {
                id.append(chrName).append(':').append(var.refStartPosition).append(":").append(var.getRefAllele());
                for (String altA : var.getAltAlleles()) {
                    id.append(":");
                    id.append(altA);
                }
                var.label = id.toString();
                id.delete(0, id.length());
            }

            for (String strGene1 : geneSet) {
                List<Variant> vars = geneVars.get(strGene1);
                Set<String> idSet = unifyingSet.get(strGene1);

                if (vars == null) {
                    vars = new ArrayList<Variant>();
                    geneVars.put(strGene1, vars);

                    idSet = new HashSet<String>();
                    unifyingSet.put(strGene1, idSet);
                }
                if (!idSet.contains(var.label)) {
                    vars.add(var);
                    idSet.add(var.label);
                }
            }

            geneSet.clear();
        }
        unifyingSet.clear();

    }

    public void summarizeVarCountsBySubject(Map<String, List<Variant>> geneVars, OpenLongObjectHashMap wahBit, List<Individual> subjectList, int[] pedEncodeGytIDMap, boolean isPhasedGty, Map<String, Integer> geneAlleCount) throws Exception {
        if (geneVars == null) {
            return;
        }
        int[] gtys = null;
        int alleleNum, base = 2;
        int subID = -1;
        int gtyID = 0;
        int count;

        boolean[] bits = new boolean[32];
        int startIndex;
        int subNum = subjectList.size();
        for (Individual indiv : subjectList) {
            if (indiv == null) {
                continue;
            }

            subID++;
            gtyID = pedEncodeGytIDMap[subID];
            if (gtyID < 0) {
                continue;
            }

            for (Map.Entry<String, List<Variant>> items : geneVars.entrySet()) {
                List<Variant> vars = items.getValue();
                count = 0;
                for (Variant var : vars) {
                    alleleNum = var.getAltAlleles().length + 1;
                    if (isPhasedGty) {
                        base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
                        if (var.compressedGtyLabel >= 0) {
                            if (var.compressedGtyLabel == 0) {
                                Arrays.fill(bits, 0, base, false);
                            } else {
                                startIndex = gtyID;
                                for (int i = 0; i < base; i++) {
                                    if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                                        bits[i] = false;
                                    } else if (startIndex < var.compressedGty[0]) {
                                        bits[i] = false;
                                    } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                                        bits[i] = true;
                                    } else if (startIndex == var.compressedGty[0]) {
                                        bits[i] = true;
                                    } else {
                                        bits[i] = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                                    }
                                    startIndex += subNum;
                                }
                            }
                            gtys = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, gtyID);
                        } else {
                            gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, gtyID, subNum);
                        }

                    } else {
                        base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
                        if (var.compressedGtyLabel >= 0) {
                            if (var.compressedGtyLabel == 0) {
                                Arrays.fill(bits, 0, base, false);
                            } else {
                                startIndex = gtyID;
                                for (int i = 0; i < base; i++) {
                                    if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                                        bits[i] = false;
                                    } else if (startIndex < var.compressedGty[0]) {
                                        bits[i] = false;
                                    } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                                        bits[i] = true;
                                    } else if (startIndex == var.compressedGty[0]) {
                                        bits[i] = true;
                                    } else {
                                        bits[i] = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                                    }
                                    startIndex += subNum;
                                }
                            }
                            gtys = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, gtyID);
                        } else {
                            gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, gtyID, subNum);
                        }

                    }

                    if (gtys != null) {
                        if (gtys[0] != 0) {
                            count++;
                        }
                        if (gtys[1] != 0) {
                            count++;
                        }
                    }
                }
                geneAlleCount.put(items.getKey(), count);
            }
        }

    }

    public void assignSelectedVar2Genesets(Map<String, GeneSet> dbPathwaySet, Map<String, List<Variant>> geneVars,
            Map<String, List<Variant>> genesetVars, byte chrID, List<Variant> invlovledVariantListst) {
        Map<String, Set<String>> unifyingSet = new HashMap<String, Set<String>>();
        Set<String> idSetAll = new HashSet<String>();
        for (Map.Entry<String, GeneSet> item : dbPathwaySet.entrySet()) {
            GeneSet refGeneSet = item.getValue();
            Set<String> geneSet = refGeneSet.getGeneSymbols();

            List<Variant> vars1 = genesetVars.get(item.getKey());
            Set<String> idSet = unifyingSet.get(item.getKey());
            if (vars1 == null) {
                vars1 = new ArrayList<Variant>();
                genesetVars.put(item.getKey(), vars1);
            }
            if (idSet == null) {
                idSet = new HashSet<String>();
                unifyingSet.put(item.getKey(), idSet);
            }
            for (String gen : geneSet) {
                List<Variant> vars = geneVars.get(gen);
                if (vars == null) {
                    continue;
                }

                for (Variant var : vars) {
                    if (!idSet.contains(var.label)) {
                        var.chrID = chrID;
                        vars1.add(var);
                        idSet.add(var.label);
                    }
                    if (!idSetAll.contains(var.label)) {
                        idSetAll.add(var.label);
                        invlovledVariantListst.add(var);
                    }
                }
            }
        }
        idSetAll.clear();
        unifyingSet.clear();

    }

    public void casecontrolUniqueModelFilterVar(Chromosome chromosome, AnnotationSummarySet ass, boolean[] uniqueFilters) {
        // AffectedRefHomGtyNum\tAffectedHetGtyNum\tAffectedAltHomGtyNum\tUnaffectedRefHomGtyNum\tUnaffectedHetGtyNum\tUnaffectedAltHomGtyNum
        int hardFilteringNum = 0;
//        Set<Byte> caseSharedHomoAllele = new HashSet<Byte>();
//        Set<Byte> caseSharedHeteAllele = new HashSet<Byte>();
//        Set<Byte> controlSharedHomoAllele = new HashSet<Byte>();
//        Set<Byte> caseSharedAllele = new HashSet<Byte>();
//        Set<Byte> controlSharedAllele = new HashSet<Byte>();

        byte[] gtys = null;
        int subID = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();
        int hitNum = 0;

        if (chromosome == null) {
            return;
        }
        tmpVarList.clear();
        for (Variant var : chromosome.variantList) {
//            caseSharedHomoAllele.clear();
//            caseSharedHeteAllele.clear();
//            controlSharedHomoAllele.clear();
//            caseSharedAllele.clear();
//            controlSharedAllele.clear();

            // "case-unique"
            if (uniqueFilters[0]) {
                if ((var.getUnaffectedHetGtyNum() <= 0 && var.getUnaffectedAltHomGtyNum() <= 0 && (var.getAffectedHetGtyNum() > 0 || var.getAffectedAltHomGtyNum() > 0))
                        || (var.getUnaffectedHetGtyNum() <= 0 && var.getUnaffectedRefHomGtyNum() <= 0 && (var.getAffectedHetGtyNum() > 0 || var.getAffectedRefHomGtyNum() > 0))) // if
                // (hetA
                // +
                // homA
                // == 0
                // ||
                // hetU
                // +
                // homU
                // > 0)
                {
                    tmpVarList.add(var);
                } else {
                    hardFilteringNum++;
                }
                // //"control-unique"
            } else if (uniqueFilters[1]) {
                if ((var.getAffectedHetGtyNum() <= 0 && var.getAffectedAltHomGtyNum() <= 0 && (var.getUnaffectedHetGtyNum() > 0 || var.getUnaffectedAltHomGtyNum() > 0))
                        || (var.getAffectedHetGtyNum() <= 0 && var.getAffectedRefHomGtyNum() <= 0 && (var.getUnaffectedHetGtyNum() > 0 || var.getUnaffectedRefHomGtyNum() > 0))) // if
                // (hetA
                // +
                // homA
                // == 0
                // ||
                // hetU
                // +
                // homU
                // > 0)
                {
                    tmpVarList.add(var);
                } else {
                    hardFilteringNum++;
                }
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVarList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVarList.size());
//        hitNum += tmpVarList.size();
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarList);
        tmpVarList.clear();

//        genome.setVarNum(hitNum);
//        genome.buildVariantIndexMapOnChromosomes();
        chromosome.buildVariantIndexMap();
//        StringBuilder message = new StringBuilder();
//        if (hardFilteringNum > 0) {
//            // message.append(hardFilteringNum).append(" variants are ignored by genotype-based hard-filtering;\n");
//        }
//        message.append(hitNum).append(" variant(s) are left after filtration according to ").append(hardFilterModel).append("!");
//        LOG.info(message.toString());
    }

    public void overlappedGeneExploreVar(Chromosome chromosome, AnnotationSummarySet ass, List<Individual> subjectList, int[] pedEncodeGytIDMap, boolean ignoreNonPathogenic, IntArrayList caseSubIDs, IntArrayList controlSubIDs, int pathogenicPredicIndex, Genome uniqueGenome) {
        String lastGeneSymb = null;
        String predicType = null;
        int caseNum = caseSubIDs.size();
        int controlNum = controlSubIDs.size();
        boolean isPhased = uniqueGenome.isIsPhasedGty();

        List<Variant> geneVars = new ArrayList<Variant>();
        List<Variant> tmpVarListIndiv = new ArrayList<Variant>();
        List<Variant> tmpVarListGene = new ArrayList<Variant>();
        List<Variant> tmpVarListChrom = new ArrayList<Variant>();
        OpenIntIntHashMap haveMutSubIDSet = new OpenIntIntHashMap();
        OpenIntIntHashMap caseUniqeAlleles = new OpenIntIntHashMap();
//        int leftVarNum = 0;

        if (chromosome == null) {
            return;
        }
        if (chromosome.variantList == null || chromosome.variantList.isEmpty()) {
            return;
        }
        lastGeneSymb = null;
        geneVars.clear();
        int subID = 0;
        int[] gty;
        boolean[] bits = new boolean[32];
        Map<String, List<Variant>> disGeneVarMap = new HashMap<String, List<Variant>>();
        List<Variant> geneDisVars = new ArrayList<Variant>();
        int startIndex;
        int subNum = subjectList.size();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                //only allow ignoreNonPathogenic to workfor misssense variants
                if (var.smallestFeatureID == 6 && ignoreNonPathogenic && pathogenicPredicIndex >= 0) {
//                    predicType = var.getFeatureValues().get(pathogenicPredicIndex);
                    predicType = var.getFeatureValues()[pathogenicPredicIndex];
                    if (predicType != null && predicType.equals("N")) {
                        continue;
                    }
                }

                geneDisVars = disGeneVarMap.get(var.geneSymb);
                if (geneDisVars == null) {
                    geneDisVars = new ArrayList<Variant>();
                    disGeneVarMap.put(var.geneSymb, geneDisVars);
                }
                geneDisVars.add(var);
            }
        }

        int effectiveCase = 0;
        for (int j = 0; j < caseNum; j++) {
            subID = caseSubIDs.getQuick(j);
            subID = pedEncodeGytIDMap[subID];
            if (subID < 0) {
                continue;
            }
            effectiveCase++;
        }

        for (Map.Entry<String, List<Variant>> item : disGeneVarMap.entrySet()) {
            lastGeneSymb = item.getKey();
            geneVars = item.getValue();
            tmpVarListGene.clear();
            haveMutSubIDSet.clear();
            caseUniqeAlleles.clear();

            for (Variant tVar : geneVars) {
                tmpVarListIndiv.clear();
                int alleleNum = tVar.getAltAlleles().length + 1;
                int base = 0;
                if (isPhased) {
                    base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
                } else {
                    base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
                }
                for (int j = 0; j < caseNum; j++) {
//                            Individual mIndiv = subjectList.get(caseSubIDs.getQuick(j));
                    subID = caseSubIDs.getQuick(j);
                    subID = pedEncodeGytIDMap[subID];
                    if (subID < 0) {
                        continue;
                    }
                    if (tVar.compressedGtyLabel >= 0) {
                        if (tVar.compressedGtyLabel == 0) {
                            Arrays.fill(bits, 0, base, false);
                        } else {
                            startIndex = subID;
                            for (int i = 0; i < base; i++) {
                                if (startIndex > tVar.compressedGty[tVar.compressedGty.length - 1]) {
                                    bits[i] = false;
                                } else if (startIndex < tVar.compressedGty[0]) {
                                    bits[i] = false;
                                } else if (startIndex == tVar.compressedGty[tVar.compressedGty.length - 1]) {
                                    bits[i] = true;
                                } else if (startIndex == tVar.compressedGty[0]) {
                                    bits[i] = true;
                                } else {
                                    bits[i] = (Arrays.binarySearch(tVar.compressedGty, startIndex) >= 0);
                                }
                                startIndex += subNum;
                            }
                        }
                        if (isPhased) {
//                                gty = mIndiv.markerGtySetArray.getPhasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                            gty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                        } else {
//                                gty = mIndiv.markerGtySetArray.getUnphasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                            gty = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                        }
                    } else if (isPhased) {
//                                gty = mIndiv.markerGtySetArray.getPhasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                        gty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, subID, subNum);
                    } else {
//                                gty = mIndiv.markerGtySetArray.getUnphasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                        gty = BinaryGtyProcessor.getUnphasedGtyAt(tVar.encodedGty, alleleNum, base, subID, subNum);
                    }

                    if (gty == null) {
                        continue;
                    }
                    if (gty[0] != gty[1]) {
                        caseUniqeAlleles.put(gty[0], 0);
                        caseUniqeAlleles.put(gty[1], 0);
                    } else {
                        caseUniqeAlleles.put(gty[0], 0);
                    }
                }

                if (caseUniqeAlleles.isEmpty()) {
                    continue;
                }
                // deliberately remove reference allele
                caseUniqeAlleles.removeKey(0);
                for (int j = 0; j < controlNum; j++) {
                    subID = controlSubIDs.getQuick(j);
                    subID = pedEncodeGytIDMap[subID];
                    if (subID < 0) {
                        continue;
                    }
                    if (tVar.compressedGtyLabel >= 0) {
                        if (tVar.compressedGtyLabel == 0) {
                            Arrays.fill(bits, 0, base, false);
                        } else {
                            startIndex = subID;
                            for (int i = 0; i < base; i++) {
                                if (startIndex > tVar.compressedGty[tVar.compressedGty.length - 1]) {
                                    bits[i] = false;
                                } else if (startIndex < tVar.compressedGty[0]) {
                                    bits[i] = false;
                                } else if (startIndex == tVar.compressedGty[tVar.compressedGty.length - 1]) {
                                    bits[i] = true;
                                } else if (startIndex == tVar.compressedGty[0]) {
                                    bits[i] = true;
                                } else {
                                    bits[i] = (Arrays.binarySearch(tVar.compressedGty, startIndex) >= 0);
                                }
                                startIndex += subNum;
                            }
                        }
//                            Individual mIndiv = subjectList.get(controlSubIDs.getQuick(j));
                        if (isPhased) {
//                                gty = mIndiv.markerGtySetArray.getPhasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                            gty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                        } else {
//                                gty = mIndiv.markerGtySetArray.getUnphasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                            gty = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                        }
                    } else if (isPhased) {
//                                gty = mIndiv.markerGtySetArray.getPhasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                        gty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, subID, subNum);
                    } else {
//                                gty = mIndiv.markerGtySetArray.getUnphasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                        gty = BinaryGtyProcessor.getUnphasedGtyAt(tVar.encodedGty, alleleNum, base, subID, subNum);
                    }

                    if (gty == null) {
                        continue;
                    }
                    if (gty[0] != gty[1]) {
                        caseUniqeAlleles.removeKey(gty[0]);
                        caseUniqeAlleles.removeKey(gty[1]);
                    } else {
                        caseUniqeAlleles.removeKey(gty[0]);
                    }
                    if (caseUniqeAlleles.isEmpty()) {
                        break;
                    }
                }
                // If this variants have no case-unqie variants
                if (caseUniqeAlleles.isEmpty()) {
                    continue;
                }
                tmpVarListGene.add(tVar);
                for (int j = 0; j < caseNum; j++) {
                    subID = caseSubIDs.getQuick(j);
                    subID = pedEncodeGytIDMap[subID];
                    if (subID < 0) {
                        continue;
                    }
                    if (tVar.compressedGtyLabel >= 0) {
                        if (tVar.compressedGtyLabel == 0) {
                            Arrays.fill(bits, 0, base, false);
                        } else {
                            startIndex = subID;
                            for (int i = 0; i < base; i++) {
                                if (startIndex > tVar.compressedGty[tVar.compressedGty.length - 1]) {
                                    bits[i] = false;
                                } else if (startIndex < tVar.compressedGty[0]) {
                                    bits[i] = false;
                                } else if (startIndex == tVar.compressedGty[tVar.compressedGty.length - 1]) {
                                    bits[i] = true;
                                } else if (startIndex == tVar.compressedGty[0]) {
                                    bits[i] = true;
                                } else {
                                    bits[i] = (Arrays.binarySearch(tVar.compressedGty, startIndex) >= 0);
                                }
                                startIndex += subNum;
                            }
                        }
                        if (isPhased) {
//                                gty = mIndiv.markerGtySetArray.getPhasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                            gty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                        } else {
//                                gty = mIndiv.markerGtySetArray.getUnphasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                            gty = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                        }
                    } else if (isPhased) {
//                                gty = mIndiv.markerGtySetArray.getPhasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                        gty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, subID, subNum);
                    } else {
//                                gty = mIndiv.markerGtySetArray.getUnphasedGtyBool(tVar.genotypeIndex, tVar.getAltAlleles().length + 1);
                        gty = BinaryGtyProcessor.getUnphasedGtyAt(tVar.encodedGty, alleleNum, base, subID, subNum);
                    }

                    if (gty == null) {
                        continue;
                    }
                    if (caseUniqeAlleles.containsKey(gty[0])) {
                        haveMutSubIDSet.put(caseSubIDs.getQuick(j), 0);
                    } else if (caseUniqeAlleles.containsKey(gty[1])) {
                        haveMutSubIDSet.put(caseSubIDs.getQuick(j), 0);
                    }
                }
            }
            // require every one has a least one case-unique
            // variants at a gene
            if (effectiveCase == haveMutSubIDSet.size()) {
                tmpVarListChrom.addAll(tmpVarListGene);
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVarListChrom.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + tmpVarListChrom.size());

        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarListChrom);
        chromosome.setHasNotOrderVariantList(true);
        tmpVarListChrom.clear();
//        leftVarNum += (chromosome.variantList.size());

        chromosome.buildVariantIndexMap();
//        genome.buildVariantIndexMapOnChromosomes();

//        StringBuilder info = new StringBuilder();
//        info.append(leftVarNum).append(" variant(s) are left after filtered by the unique variants on gene level.\n");
//
//        LOG.info(info);
//        genome.setVarNum(leftVarNum);
    }

    public int[] partitionEvenBlock(int totalSnpSize, int blockNum) throws Exception {
        int[] bigBlockIndexes;
        int intervalLen = totalSnpSize / blockNum;
        if (intervalLen == 0) {
            bigBlockIndexes = new int[totalSnpSize + 1];
            for (int i = 0; i < bigBlockIndexes.length; i++) {
                bigBlockIndexes[i] = i;
            }
        } else {
            bigBlockIndexes = new int[blockNum + 1];
            Arrays.fill(bigBlockIndexes, 0);
            for (int s = 1; s < blockNum; s++) {
                bigBlockIndexes[s] = s * intervalLen;
            }
            bigBlockIndexes[blockNum] = totalSnpSize;
        }

        return bigBlockIndexes;
    }

    private void formatBialleleicUnphasedBits(List<Variant> variantList, int totalPedSubjectNum, int[][] bits1,
            int[][] bits2, int[][] bits3, int[][] bits4, double[] mean1, int[] sum12, boolean[] hasMissingGty) {
        int varSize = variantList.size();
        int firstBitByteEnd;

        final int gtyLen = 32;
        int unitNum;
        int firstTailLen = totalPedSubjectNum % 8;
        int offsetTmp = 0;
        if (totalPedSubjectNum % gtyLen != 0) {
            unitNum = totalPedSubjectNum / gtyLen + 1;
        } else {
            unitNum = totalPedSubjectNum / gtyLen;
        }

        if (totalPedSubjectNum % 8 != 0) {
            firstBitByteEnd = totalPedSubjectNum / 8 + 1;
        } else {
            firstBitByteEnd = totalPedSubjectNum / 8;
        }

        int base = 2;
        int secondBitByteStart = (totalPedSubjectNum) / 8;
        int secondBitBytEnd = totalPedSubjectNum * 2 / 8;
        if ((totalPedSubjectNum * 2) % 8 != 0) {
            secondBitBytEnd += 1;
        }

        int thirdBitByteStart = (totalPedSubjectNum * 2) / 8;
        int thirdBitBytEnd = (totalPedSubjectNum * 3) / 8;
        if ((totalPedSubjectNum * 3) % 8 != 0) {
            thirdBitBytEnd += 1;
        }

        int wordNum = firstBitByteEnd / 4;
        int maskInt = totalPedSubjectNum % gtyLen;
        int c = 0;
        if (maskInt != 0) {
            c = (0x80000000 >> (maskInt - 1));
            c = ~c;
            //System.out.println(Integer.toBinaryString(c));
        }
        int nonNumIndiv = 0;
        int startIndex;
        boolean isOne;
        int subID;
        int index1, index2;
        for (int i = 0; i < varSize; i++) {
            Variant var = variantList.get(i);
            if (var.getAltAlleles().length > 1) {
                continue;
            }

            bits1[i] = new int[unitNum];
            bits2[i] = new int[unitNum];
            bits3[i] = new int[unitNum];
            bits4[i] = new int[unitNum];

            Arrays.fill(bits1[i], 0);
            Arrays.fill(bits2[i], 0);
            Arrays.fill(bits3[i], 0);
            Arrays.fill(bits4[i], 0);
            if (var.compressedGtyLabel >= 0) {
                for (subID = 0; subID < totalPedSubjectNum; subID++) {
                    startIndex = subID;
                    if (var.compressedGtyLabel == 0) {
                        isOne = false;
                    } else if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = false;
                    } else if (startIndex < var.compressedGty[0]) {
                        isOne = false;
                    } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = true;
                    } else if (startIndex == var.compressedGty[0]) {
                        isOne = true;
                    } else {
                        isOne = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                    }
                    if (isOne) {
                        index1 = subID / gtyLen;
                        index2 = subID % gtyLen;
                        bits1[i][index1] |= GlobalManager.intOpers[index2];
                    }
                }

                for (subID = 0; subID < totalPedSubjectNum; subID++) {
                    startIndex = subID + totalPedSubjectNum;
                    if (var.compressedGtyLabel == 0) {
                        isOne = false;
                    } else if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = false;
                    } else if (startIndex < var.compressedGty[0]) {
                        isOne = false;
                    } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = true;
                    } else if (startIndex == var.compressedGty[0]) {
                        isOne = true;
                    } else {
                        isOne = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                    }
                    if (isOne) {
                        index1 = subID / gtyLen;
                        index2 = subID % gtyLen;
                        bits3[i][index1] |= GlobalManager.intOpers[index2];
                    }
                }

                if (maskInt != 0) {
                    //System.out.println(Integer.toBinaryString(bits1[t][wordNum]));
                    bits1[i][wordNum] |= c;
                    //System.out.println(Integer.toBinaryString(bits1[t][wordNum]));
                }
                if (maskInt != 0) {
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                    bits3[i][wordNum] |= c;
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                }
            } else {
                byte[] bytes = var.encodedGty;
                for (int j = 0; j < firstBitByteEnd; j++) {
                    wordNum = j >>> 2;
                    offsetTmp = j - (wordNum << 2);
                    bits1[i][wordNum] |= ((bytes[j] & 0xFF) << (24 - (offsetTmp << 3)));
                    //  System.out.println(String.format("%8s", Integer.toBinaryString(bytes[j] & 0xFF)).replace(' ', '0'));
                }
                if (maskInt != 0) {
                    //System.out.println(Integer.toBinaryString(bits1[t][wordNum]));
                    bits1[i][wordNum] |= c;
                    //System.out.println(Integer.toBinaryString(bits1[t][wordNum]));
                }

                /*
            for (int j = 0; j < unitNum; j++) {
                System.out.println(String.format("%32s", Integer.toBinaryString(bits1[t][j])).replace(' ', '0'));
            }
                 */
                for (int j = secondBitByteStart; j < secondBitBytEnd; j++) {
                    wordNum = (j - secondBitByteStart) >>> 2;
                    if (wordNum == unitNum) {
                        wordNum--;
                        // System.out.println(String.format("%8s", Integer.toBinaryString(bytes[j] & 0xFF)).replace(' ', '0'));
                        break;
                    }
                    offsetTmp = (j - secondBitByteStart) - (wordNum << 2);
                    bits3[i][wordNum] |= (((bytes[j] & 0xFF) << (24 - (offsetTmp << 3) + firstTailLen)));
                    // System.out.println(String.format("%32s", Integer.toBinaryString(bits3[t][wordNum])).replace(' ', '0'));
                    if (offsetTmp == 3 && firstTailLen != 0) {
                        // c = (0x80000000 >> (32 - 8 + firstTailLen - 1));
                        //System.out.println(String.format("%8s", Integer.toBinaryString(c)).replace(' ', '0'));
                        bits3[i][wordNum] |= (((bytes[j + 1] & 0xFF)) >>> (8 - firstTailLen));
                        //byte kk = (byte) ((bytes[j + 1] & 0xFF) & c);
                        //System.out.println(String.format("%8s", Integer.toBinaryString(kk).replace(' ', '0')));
                    }
                    // System.out.println(String.format("%8s", Integer.toBinaryString(bytes[j] & 0xFF)).replace(' ', '0'));
                }

                //it is necessary when they are more than two bits
                if (maskInt != 0) {
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                    bits3[i][wordNum] |= c;
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                }
                /*
            for (int j = 0; j < unitNum; j++) {
                System.out.println(String.format("%32s", Integer.toBinaryString(bits3[t][j])).replace(' ', '0'));
            }
            System.out.println();
                 */
            }
            nonNumIndiv = 0;
            for (int j = 0; j < unitNum; j++) {
                bits4[i][j] = ~(bits1[i][j] & bits3[i][j]);
                bits1[i][j] = bits1[i][j] & bits4[i][j];
                bits3[i][j] = bits3[i][j] & bits4[i][j];
                bits2[i][j] = (bits1[i][j] | bits3[i][j]);
                mean1[i] += ((Integer.bitCount(bits1[i][j]) << 1) + Integer.bitCount(bits3[i][j]));
                sum12[i] += ((Integer.bitCount(bits1[i][j]) << 2) + Integer.bitCount(bits3[i][j]));
                nonNumIndiv += Integer.bitCount(bits4[i][j]);
                //System.out.println(String.format("%32s", Integer.toBinaryString(bits4[i][j])).replace(' ', '0'));
            }
            hasMissingGty[i] = nonNumIndiv != totalPedSubjectNum;
        }

    }

    private void formatBialleleicPhasedBits(List<Variant> variantList, int totalPedSubjectNum, int[][] bits1,
            int[][] bits2, int[][] bits3, double[] mean1, boolean[] hasMissingGty) {
        int varSize = variantList.size();
        int firstBitByteEnd;

        final int gtyLen = 32;
        int unitNum;
        int firstTailLen = totalPedSubjectNum % 8;
        int offsetTmp = 0;
        if (totalPedSubjectNum % gtyLen != 0) {
            unitNum = totalPedSubjectNum / gtyLen + 1;
        } else {
            unitNum = totalPedSubjectNum / gtyLen;
        }

        if (totalPedSubjectNum % 8 != 0) {
            firstBitByteEnd = totalPedSubjectNum / 8 + 1;
        } else {
            firstBitByteEnd = totalPedSubjectNum / 8;
        }

        int secondTailLen = totalPedSubjectNum * 2 % 8;
        int base = 2;
        int secondBitByteStart = (totalPedSubjectNum) / 8;
        int secondBitBytEnd = totalPedSubjectNum * 2 / 8;
        if ((totalPedSubjectNum * 2) % 8 != 0) {
            secondBitBytEnd += 1;
        }

        int thirdBitByteStart = (totalPedSubjectNum * 2) / 8;
        int thirdBitBytEnd = totalPedSubjectNum * 3 / 8;
        if ((totalPedSubjectNum * 3) % 8 != 0) {
            thirdBitBytEnd += 1;
        }

        int wordNum = firstBitByteEnd / 4;
        int maskInt = totalPedSubjectNum % gtyLen;
        int c = 0;
        if (maskInt != 0) {
            c = (0x80000000 >> (maskInt - 1));
            c = ~c;
            //System.out.println(Integer.toBinaryString(c));
        }
        int nonNumIndiv = 0;
        int startIndex;
        boolean isOne;
        int subID;
        int index1, index2;
        for (int i = 0; i < varSize; i++) {
            Variant var = variantList.get(i);
            if (var.getAltAlleles().length > 1) {
                continue;
            }

            bits1[i] = new int[unitNum];
            bits2[i] = new int[unitNum];
            bits3[i] = new int[unitNum];
            Arrays.fill(bits1[i], 0);
            Arrays.fill(bits2[i], 0);
            Arrays.fill(bits3[i], 0);

            if (var.compressedGtyLabel >= 0) {
                for (subID = 0; subID < totalPedSubjectNum; subID++) {
                    startIndex = subID;
                    if (var.compressedGtyLabel == 0) {
                        isOne = false;
                    } else if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = false;
                    } else if (startIndex < var.compressedGty[0]) {
                        isOne = false;
                    } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = true;
                    } else if (startIndex == var.compressedGty[0]) {
                        isOne = true;
                    } else {
                        isOne = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                    }
                    if (isOne) {
                        index1 = subID / gtyLen;
                        index2 = subID % gtyLen;
                        bits3[i][index1] |= GlobalManager.intOpers[index2];
                    }
                }
                for (subID = 0; subID < totalPedSubjectNum; subID++) {
                    startIndex = subID + totalPedSubjectNum;
                    if (var.compressedGtyLabel == 0) {
                        isOne = false;
                    } else if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = false;
                    } else if (startIndex < var.compressedGty[0]) {
                        isOne = false;
                    } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = true;
                    } else if (startIndex == var.compressedGty[0]) {
                        isOne = true;
                    } else {
                        isOne = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                    }
                    if (isOne) {
                        index1 = subID / gtyLen;
                        index2 = subID % gtyLen;
                        bits1[i][index1] |= GlobalManager.intOpers[index2];
                    }
                }

                for (subID = 0; subID < totalPedSubjectNum; subID++) {
                    startIndex = subID + totalPedSubjectNum + totalPedSubjectNum;
                    if (var.compressedGtyLabel == 0) {
                        isOne = false;
                    } else if (startIndex > var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = false;
                    } else if (startIndex < var.compressedGty[0]) {
                        isOne = false;
                    } else if (startIndex == var.compressedGty[var.compressedGty.length - 1]) {
                        isOne = true;
                    } else if (startIndex == var.compressedGty[0]) {
                        isOne = true;
                    } else {
                        isOne = (Arrays.binarySearch(var.compressedGty, startIndex) >= 0);
                    }
                    if (isOne) {
                        index1 = subID / gtyLen;
                        index2 = subID % gtyLen;
                        bits2[i][index1] |= GlobalManager.intOpers[index2];
                    }
                }
                if (maskInt != 0) {
                    //System.out.println(Integer.toBinaryString(c));
                    bits1[i][wordNum] |= c;
                    //System.out.println(Integer.toBinaryString(bits1[t][wordNum]));
                }
                if (maskInt != 0) {
                    //System.out.println(Integer.toBinaryString(bits1[t][wordNum]));
                    bits2[i][wordNum] |= c;
                    //System.out.println(Integer.toBinaryString(bits1[t][wordNum]));
                }
                if (maskInt != 0) {
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                    bits3[i][wordNum] |= c;
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                }
            } else {
                byte[] bytes = var.encodedGty;
                for (int j = 0; j < firstBitByteEnd; j++) {
                    wordNum = j >>> 2;
                    offsetTmp = j - (wordNum << 2);
                    bits3[i][wordNum] |= ((bytes[j] & 0xFF) << (24 - (offsetTmp << 3)));
                    //  System.out.println(String.format("%8s", Integer.toBinaryString(bytes[j] & 0xFF)).replace(' ', '0'));
                }
                if (maskInt != 0) {
                    //System.out.println(Integer.toBinaryString(bits1[t][wordNum]));
                    bits3[i][wordNum] |= c;
                    //System.out.println(Integer.toBinaryString(bits1[t][wordNum]));
                }

                /*
            for (int j = 0; j < unitNum; j++) {
                System.out.println(String.format("%32s", Integer.toBinaryString(bits1[t][j])).replace(' ', '0'));
            }
                 */
                for (int j = secondBitByteStart; j < secondBitBytEnd; j++) {
                    wordNum = (j - secondBitByteStart) >>> 2;
                    if (wordNum == unitNum) {
                        wordNum--;
                        // System.out.println(String.format("%8s", Integer.toBinaryString(bytes[j] & 0xFF)).replace(' ', '0'));
                        break;
                    }
                    offsetTmp = (j - secondBitByteStart) - (wordNum << 2);
                    bits1[i][wordNum] |= (((bytes[j] & 0xFF) << (24 - (offsetTmp << 3) + firstTailLen)));
                    // System.out.println(String.format("%32s", Integer.toBinaryString(bits3[t][wordNum])).replace(' ', '0'));
                    if (offsetTmp == 3 && firstTailLen != 0) {
                        // c = (0x80000000 >> (32 - 8 + firstTailLen - 1));
                        //System.out.println(String.format("%8s", Integer.toBinaryString(c)).replace(' ', '0'));
                        bits1[i][wordNum] |= (((bytes[j + 1] & 0xFF)) >>> (8 - firstTailLen));
                        //byte kk = (byte) ((bytes[j + 1] & 0xFF) & c);
                        //System.out.println(String.format("%8s", Integer.toBinaryString(kk).replace(' ', '0')));
                    }
                    // System.out.println(String.format("%8s", Integer.toBinaryString(bytes[j] & 0xFF)).replace(' ', '0'));
                }

                //it is necessary when they are more than two bits
                if (maskInt != 0) {
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                    bits1[i][wordNum] |= c;
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                }

                for (int j = thirdBitByteStart; j < thirdBitBytEnd; j++) {
                    wordNum = (j - thirdBitByteStart) >>> 2;
                    if (wordNum == unitNum) {
                        wordNum--;
                        // System.out.println(String.format("%8s", Integer.toBinaryString(bytes[j] & 0xFF)).replace(' ', '0'));
                        break;
                    }
                    offsetTmp = (j - thirdBitByteStart) - (wordNum << 2);
                    bits2[i][wordNum] |= (((bytes[j] & 0xFF) << (24 - (offsetTmp << 3) + secondTailLen)));
                    // System.out.println(String.format("%32s", Integer.toBinaryString(bits3[t][wordNum])).replace(' ', '0'));
                    if (offsetTmp == 3 && firstTailLen != 0) {
                        // c = (0x80000000 >> (32 - 8 + firstTailLen - 1));
                        //System.out.println(String.format("%8s", Integer.toBinaryString(c)).replace(' ', '0'));
                        bits2[i][wordNum] |= (((bytes[j + 1] & 0xFF)) >>> (8 - secondTailLen));
                        //byte kk = (byte) ((bytes[j + 1] & 0xFF) & c);
                        //System.out.println(String.format("%8s", Integer.toBinaryString(kk).replace(' ', '0')));
                    }
                    // System.out.println(String.format("%8s", Integer.toBinaryString(bytes[j] & 0xFF)).replace(' ', '0'));
                }

                //it is necessary when they are more than two bits
                if (maskInt != 0) {
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                    bits2[i][wordNum] |= c;
                    // System.out.println(Integer.toBinaryString(bits3[t][wordNum]));
                }

                /*
            for (int j = 0; j < unitNum; j++) {
                System.out.println(String.format("%32s", Integer.toBinaryString(bits3[t][j])).replace(' ', '0'));
            }
            System.out.println();
                 */
            }
            nonNumIndiv = 0;
            for (int j = 0; j < unitNum; j++) {
                bits3[i][j] = ~bits3[i][j];
                bits1[i][j] = (bits1[i][j] & bits3[i][j]);
                bits2[i][j] = (bits2[i][j] & bits3[i][j]);
                mean1[i] += (Integer.bitCount(bits1[i][j]) + Integer.bitCount(bits2[i][j]));
                nonNumIndiv += Integer.bitCount(bits3[i][j]);
                //System.out.println(String.format("%32s", Integer.toBinaryString(bits1[i][j])).replace(' ', '0'));
                //System.out.println(String.format("%32s", Integer.toBinaryString(bits2[i][j])).replace(' ', '0'));
                //System.out.println(String.format("%32s", Integer.toBinaryString(bits3[i][j])).replace(' ', '0'));
            }
            hasMissingGty[i] = nonNumIndiv != totalPedSubjectNum;
            mean1[i] /= (nonNumIndiv << 1);
        }

    }

    private double getCorrelation(DoubleArrayList xVect, DoubleArrayList yVect) {
        double meanX = 0.0, meanY = 0.0;
        for (int i = 0; i < xVect.size(); i++) {
            meanX += xVect.getQuick(i);
            meanY += yVect.getQuick(i);
        }

        meanX /= xVect.size();
        meanY /= yVect.size();

        double sumXY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
        for (int i = 0; i < xVect.size(); i++) {
            sumXY += ((xVect.getQuick(i) - meanX) * (yVect.getQuick(i) - meanY));
            sumX2 += Math.pow(xVect.getQuick(i) - meanX, 2.0);
            sumY2 += Math.pow(yVect.getQuick(i) - meanY, 2.0);
        }

        return (sumXY / (Math.sqrt(sumX2) * Math.sqrt(sumY2)));
    }

    private void formatUnphasedGenotypesCalCorr(List<Variant> variantList, int[] pedEncodeGytIDMap) {
        int varSize = variantList.size();
        int base = 2;
        int alleleNum = 2;
        int gtyID;
        int[] gtys;
        int[][] genotypes = new int[varSize][];
        boolean[][] genotypeMissing = new boolean[varSize][];
        int totalPedSubjectNum = pedEncodeGytIDMap.length;
        int nonMssingNum = 0;
        for (int t = 0; t < varSize; t++) {
            Variant var = variantList.get(t);
            if (var.getAltAlleles().length > 1) {
                continue;
            }
            genotypeMissing[t] = new boolean[totalPedSubjectNum];
            genotypes[t] = new int[totalPedSubjectNum];
            Arrays.fill(genotypes[t], 0);
            Arrays.fill(genotypeMissing[t], false);
            nonMssingNum = 0;
            for (int i = 0; i < totalPedSubjectNum; i++) {
                gtyID = pedEncodeGytIDMap[i];
                gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, gtyID, totalPedSubjectNum);
                if (gtys == null) {
                    genotypeMissing[t][i] = true;
                    continue;
                } else if (gtys[0] == 1 && gtys[1] == 1) {
                    genotypes[t][i] = 2;
                    // mean1[t] += 2;
                    //  sum12[t] += 4;
                } else if (gtys[0] == 1 || gtys[1] == 1) {
                    genotypes[t][i] = 1;
                    //  mean1[t]++;
                    //  sum12[t]++;
                }
                nonMssingNum++;
            }
            // mean1[t] /= Math.sqrt(nonNumIndiv);
            //sum12[t] = Math.sqrt(sum12[t] - mean1[t] * mean1[t]);
        }
        DoubleArrayList groupA = new DoubleArrayList();
        DoubleArrayList groupB = new DoubleArrayList();
        double meanA, meanB, r;
        Long startTime = System.nanoTime();
        for (int t = 0; t < varSize; t++) {
            for (int i = t + 1; i < varSize; i++) {
                groupB.clear();
                groupA.clear();
                for (int k = 0; k < totalPedSubjectNum; k++) {
                    if (!genotypeMissing[t][k] && !genotypeMissing[i][k]) {
                        groupA.add(genotypes[t][k]);
                        groupB.add(genotypes[i][k]);
                    }
                }
                if (groupB.isEmpty()) {
                    continue;
                }
                //meanA = Descriptive.mean(groupA);
                // meanB = Descriptive.mean(groupB);

                r = getCorrelation(groupA, groupB);
                if (r > 0.05) {
                    System.out.print(variantList.get(t).refStartPosition + "\t" + variantList.get(i).refStartPosition + "\t" + r + "\n");
                }
            }
            //System.out.println();
        }
        // System.out.println("formatUnphasedGenotypesCalCorr:" + (System.nanoTime() - startTime) / 1000000000.0);
    }

    private void formatPhasedGenotypesCalLD(List<Variant> variantList, int[] pedEncodeGytIDMap) {
        int varSize = variantList.size();
        int base = 3;
        int alleleNum = 2;
        int gtyID;
        int[] gtys;
        int totalPedSubjectNum = pedEncodeGytIDMap.length;
        int[][] genotypes = new int[totalPedSubjectNum * 2][];
        boolean[][] genotypeMissing = new boolean[totalPedSubjectNum][];

        for (int i = 0; i < totalPedSubjectNum; i++) {
            genotypes[i * 2] = new int[varSize];
            genotypes[i * 2 + 1] = new int[varSize];
            genotypeMissing[i] = new boolean[varSize];
            Arrays.fill(genotypeMissing[i], false);
            Arrays.fill(genotypes[i * 2], 0);
            Arrays.fill(genotypes[i * 2 + 1], 0);
            gtyID = pedEncodeGytIDMap[i];

            for (int t = 0; t < varSize; t++) {
                Variant var = variantList.get(t);
                if (var.getAltAlleles().length > 1) {
                    continue;
                }
                gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, gtyID, totalPedSubjectNum);
                if (gtys == null) {
                    genotypeMissing[i][t] = true;
                } else {
                    genotypes[i * 2][t] = gtys[0];
                    genotypes[i * 2 + 1][t] = gtys[1];
                }
            }
        }

        double meanA, meanB, r;
        Long startTime = System.nanoTime();
        double freqAB = 0;
        double freqA = 0;
        double freqB = 0;
        int nonMissing = 0;
        for (int t = 0; t < varSize; t++) {
            for (int i = t + 1; i < varSize; i++) {
                freqAB = 0;
                freqA = 0;
                freqB = 0;
                nonMissing = 0;
                for (int k = 0; k < totalPedSubjectNum; k++) {
                    if (!genotypeMissing[k][t] && !genotypeMissing[k][i]) {
                        if (genotypes[k * 2][t] == 1) {
                            freqA += 1;
                        }
                        if (genotypes[k * 2][i] == 1) {
                            freqB += 1;
                        }
                        if (genotypes[k * 2][t] == 1 && genotypes[k * 2][i] == 1) {
                            freqAB += 1;
                        }
                        if (genotypes[k * 2 + 1][t] == 1) {
                            freqA += 1;
                        }
                        if (genotypes[k * 2 + 1][i] == 1) {
                            freqB += 1;
                        }
                        if (genotypes[k * 2 + 1][t] == 1 && genotypes[k * 2 + 1][i] == 1) {
                            freqAB += 1;
                        }
                        nonMissing += 2;
                    }
                }
                freqA = freqA / nonMissing;
                freqB = freqB / nonMissing;
                freqAB = freqAB / nonMissing;
                r = freqAB - freqA * freqB;
                r = r / Math.sqrt(freqA * (1 - freqA) * freqB * (1 - freqB));
                if (r > 0.05) {
                    System.out.print(variantList.get(t).refStartPosition + "\t" + variantList.get(i).refStartPosition + "\t" + r + "\n");
                }
            }
            //System.out.println();
        }
        // System.out.println("formatUnphasedGenotypesCalCorr:" + (System.nanoTime() - startTime) / 1000000000.0);
    }

    public void calculateLDVar(Chromosome chromosome, List<Individual> subjectList,
            int[] pedEncodeGytIDMap, boolean isPhased, int threadNum, int windowWidth, double minR2, BufferedWriter bwLD) throws Exception {
        if (chromosome == null) {
            return;
        }
        if (chromosome.variantList == null || chromosome.variantList.isEmpty()) {
            return;
        }
        List<Variant> selectedVariantList = new ArrayList<Variant>();
        for (Variant var : chromosome.variantList) {
            //onle consider biallelic variants
            if (var.getAltAlleles().length > 1) {
                continue;
            }
            selectedVariantList.add(var);
        }
        int[] positions = new int[selectedVariantList.size()];
        for (int i = 0; i < positions.length; i++) {
            positions[i] = selectedVariantList.get(i).refStartPosition;
        }
        String chromName = chromosome.getName();
        int varSize = selectedVariantList.size();
        int[][] bits1 = new int[varSize][];
        int[][] bits2 = new int[varSize][];
        int[][] bits3 = new int[varSize][];
        int[][] bits4 = null;
        double[] mean1 = new double[varSize];
        int[] sum12 = null;
        boolean[] hasMissingGty = new boolean[varSize];
        int totalPedSubjectNum = pedEncodeGytIDMap.length;
        if (isPhased) {
            formatBialleleicPhasedBits(selectedVariantList, totalPedSubjectNum, bits1, bits2, bits3, mean1, hasMissingGty);
            //formatPhasedGenotypesCalLD(selectedVariantList, pedEncodeGytIDMap);
        } else {
            bits4 = new int[varSize][];
            sum12 = new int[varSize];
            formatBialleleicUnphasedBits(selectedVariantList, totalPedSubjectNum, bits1, bits2, bits3, bits4, mean1, sum12, hasMissingGty);
            //formatUnphasedGenotypesCalCorr(selectedVariantList, pedEncodeGytIDMap);
        }

        //cannot be too large
        //double[] corrs = new double[varSize * (varSize - 1) / 2];
        int runningThread = 0;
        // int[] smallBlockIndexes = cc.partitionEvenBlock(size, threadNum > 50 ? threadNum : 50);
        int[] smallBlockIndexes = partitionEvenBlock(varSize, threadNum * 5);
        int smallBlockNum = smallBlockIndexes.length - 1;
        int lastButOneIndex = smallBlockIndexes[smallBlockNum - 1];
        int lastBlockLen = smallBlockIndexes[smallBlockNum] - smallBlockIndexes[smallBlockNum - 1];

        int blockLen = smallBlockIndexes[1] - smallBlockIndexes[0];
        int ti = 0;
        ExecutorService exec = Executors.newFixedThreadPool(threadNum);
        CompletionService serv = new ExecutorCompletionService(exec);
        Long startTime = System.nanoTime();
        for (int i = 0; i < smallBlockNum; i++) {
            final double[][] rArray = new double[smallBlockNum - i][];
            runningThread = 0;
            for (int j = i; j < smallBlockNum; j++) {
                //always prepare a rectangle
                rArray[j - i] = new double[(smallBlockIndexes[i + 1] - smallBlockIndexes[i]) * (smallBlockIndexes[j + 1] - smallBlockIndexes[j])];
                LDCalcTask task = new LDCalcTask(bits1, bits2, bits3, bits4, mean1, sum12, hasMissingGty, totalPedSubjectNum, smallBlockIndexes[i], smallBlockIndexes[i + 1],
                        smallBlockIndexes[j], smallBlockIndexes[j + 1], rArray[j - i], isPhased);
                task.setPositionCutoff(positions, windowWidth);
                serv.submit(task);
                runningThread++;
            }

            for (int index = 0; index < runningThread; index++) {
                Future task1 = serv.take();
                String infor1 = (String) task1.get();
                // System.out.println(infor);
            }

            //unfortunately the last block often have uneven size 
            //correct for dependency
            double r;
            for (int t = smallBlockIndexes[i]; t < smallBlockIndexes[i + 1]; t++) {
                for (int k = t + 1; k < varSize; k++) {
                    if (k < lastButOneIndex) {
                        r = rArray[(k - smallBlockIndexes[i]) / blockLen][(t - smallBlockIndexes[i]) * blockLen + (k - smallBlockIndexes[i]) % blockLen];
                    } else {
                        r = rArray[smallBlockNum - i - 1][(t - smallBlockIndexes[i]) * lastBlockLen + k - lastButOneIndex];
                    }
                    // corrs[ti] = r;
                    r = r * r;
                    if (r > minR2) {
                        bwLD.write(chromName + "\t" + selectedVariantList.get(t).refStartPosition + "\t" + selectedVariantList.get(k).refStartPosition + "\t" + r + "\n");
                    }
                    // tmpFileWriter.write(allData.get(t).id + "\t" + allData.get(k).id + "\t" + r + "\n");
                    //System.out.print("(" + t + "," + k + ":" + r + ")");
                    //  System.out.println(ti + ":" + r);
                    ti++;
                }
                //System.out.println();
            }
        }//end of  for (int t = 0; t < smallBlockNum; t++)
        // System.out.println("LDCalcTask:" + (System.nanoTime() - startTime) / 1000000000.0);
        exec.shutdown();

    }

    public void LDPruneWithinDistance(Chromosome chromosome, List<Individual> subjectList,
            int[] pedEncodeGytIDMap, boolean isPhased, int threadNum, int windowWidth, double r2, AnnotationSummarySet ass) throws Exception {

        if (chromosome == null) {
            return;
        }
        if (chromosome.variantList == null || chromosome.variantList.isEmpty()) {
            return;
        }
        List<Variant> selectedVariantList = new ArrayList<Variant>();

        for (Variant var : chromosome.variantList) {
            //onle consider biallelic variants
            if (var.getAltAlleles().length > 1) {
                continue;
            }
            selectedVariantList.add(var);
        }
        int[] positions = new int[selectedVariantList.size()];
        for (int i = 0; i < positions.length; i++) {
            positions[i] = selectedVariantList.get(i).refStartPosition;
        }
        String chromName = chromosome.getName();
        int varSize = selectedVariantList.size();
        int[][] bits1 = new int[varSize][];
        int[][] bits2 = new int[varSize][];
        int[][] bits3 = new int[varSize][];
        int[][] bits4 = null;
        double[] mean1 = new double[varSize];
        int[] sum12 = null;
        boolean[] hasMissingGty = new boolean[varSize];
        int totalPedSubjectNum = pedEncodeGytIDMap.length;
        if (isPhased) {
            formatBialleleicPhasedBits(selectedVariantList, totalPedSubjectNum, bits1, bits2, bits3, mean1, hasMissingGty);
            //formatPhasedGenotypesCalLD(selectedVariantList, pedEncodeGytIDMap);
        } else {
            bits4 = new int[varSize][];
            sum12 = new int[varSize];
            formatBialleleicUnphasedBits(selectedVariantList, totalPedSubjectNum, bits1, bits2, bits3, bits4, mean1, sum12, hasMissingGty);
            //formatUnphasedGenotypesCalCorr(selectedVariantList, pedEncodeGytIDMap);
        }

        //cannot be too large
        //double[] corrs = new double[varSize * (varSize - 1) / 2];
        int runningThread = 0;
        // int[] smallBlockIndexes = cc.partitionEvenBlock(size, threadNum > 50 ? threadNum : 50);
        int[] smallBlockIndexes = partitionEvenBlock(varSize, threadNum * 5);
        int smallBlockNum = smallBlockIndexes.length - 1;
        int lastButOneIndex = smallBlockIndexes[smallBlockNum - 1];
        int lastBlockLen = smallBlockIndexes[smallBlockNum] - smallBlockIndexes[smallBlockNum - 1];

        int blockLen = smallBlockIndexes[1] - smallBlockIndexes[0];

        ExecutorService exec = Executors.newFixedThreadPool(threadNum);
        CompletionService serv = new ExecutorCompletionService(exec);
        OpenIntIntHashMap redundantVars = new OpenIntIntHashMap();
        for (int i = 0; i < smallBlockNum; i++) {
            final double[][] rArray = new double[smallBlockNum - i][];
            runningThread = 0;
            for (int j = i; j < smallBlockNum; j++) {
                //always prepare a rectangle
                rArray[j - i] = new double[(smallBlockIndexes[i + 1] - smallBlockIndexes[i]) * (smallBlockIndexes[j + 1] - smallBlockIndexes[j])];
                Arrays.fill(rArray[j - i], 0);
                LDCalcTask task = new LDCalcTask(bits1, bits2, bits3, bits4, mean1, sum12, hasMissingGty, totalPedSubjectNum, smallBlockIndexes[i], smallBlockIndexes[i + 1],
                        smallBlockIndexes[j], smallBlockIndexes[j + 1], rArray[j - i], isPhased);
                task.setPositionCutoff(positions, windowWidth);
                serv.submit(task);
                runningThread++;
            }

            for (int index = 0; index < runningThread; index++) {
                Future task1 = serv.take();
                String infor1 = (String) task1.get();
                // System.out.println(infor);
            }

            //unfortunately the last block often have uneven size 
            //correct for dependency
            double r;
            for (int t = smallBlockIndexes[i]; t < smallBlockIndexes[i + 1]; t++) {
                for (int k = t + 1; k < varSize; k++) {
                    if (k < lastButOneIndex) {
                        r = rArray[(k - smallBlockIndexes[i]) / blockLen][(t - smallBlockIndexes[i]) * blockLen + (k - smallBlockIndexes[i]) % blockLen];
                    } else {
                        r = rArray[smallBlockNum - i - 1][(t - smallBlockIndexes[i]) * lastBlockLen + k - lastButOneIndex];
                    }
                    // corrs[ti] = r;
                    r = r * r;
                    if (r >= r2) {
                        if (!redundantVars.containsKey(selectedVariantList.get(t).refStartPosition) && !redundantVars.containsKey(selectedVariantList.get(k).refStartPosition)) {
                            redundantVars.put(selectedVariantList.get(k).refStartPosition, 0);
                        }
                    }
                    // tmpFileWriter.write(allData.get(t).id + "\t" + allData.get(k).id + "\t" + r + "\n");
                    //System.out.print("(" + t + "," + k + ":" + r + ")");
                    //  System.out.println(ti + ":" + r);

                }
                //System.out.println();
            }
        }//end of  for (int t = 0; t < smallBlockNum; t++)
        // System.out.println("LDCalcTask:" + (System.nanoTime() - startTime) / 1000000000.0);
        exec.shutdown();

        if (!redundantVars.isEmpty()) {
            selectedVariantList.clear();
            selectedVariantList.addAll(chromosome.variantList);
            chromosome.variantList.clear();
            for (Variant var : selectedVariantList) {
                //onle consider biallelic variants
                if (!redundantVars.containsKey(var.refStartPosition)) {
                    chromosome.variantList.add(var);
                }
            }
            //System.out.println("On chromosome " + chromName + ", " + redundantVars.size() + " variants in high LD r2>" + r2 + " are removed!");
            redundantVars.clear();
            chromosome.buildVariantIndexMap();
        }
        //  System.out.println(chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size());
    }

    public void calculateLDVarComp(Chromosome chromosome, OpenLongObjectHashMap wahBit, List<Individual> subjectList,
            int[] pedEncodeGytIDMap, boolean isPhased, int threadNum) throws Exception {
        if (chromosome == null) {
            return;
        }
        if (chromosome.variantList == null || chromosome.variantList.isEmpty()) {
            return;
        }
        List<Variant> selectedVariantList = new ArrayList<Variant>();
        for (Variant var : chromosome.variantList) {
            //onle consider biallelic variants
            if (var.getAltAlleles().length > 1) {
                continue;
            }
            selectedVariantList.add(var);
        }
        List<Variant> runningList = new ArrayList<Variant>();
        int testingSize = 1000;
        ExecutorService exec = Executors.newFixedThreadPool(threadNum);
        CompletionService serv = new ExecutorCompletionService(exec);
        do {
            for (int s = 0; s < testingSize; s++) {
                runningList.add(selectedVariantList.get(s));
            }
            System.out.println("Testing size: " + testingSize);
            int varSize = runningList.size();
            int[][] bits1 = new int[varSize][];
            int[][] bits2 = new int[varSize][];
            int[][] bits3 = new int[varSize][];
            int[][] bits4 = new int[varSize][];
            double[] sum1 = new double[varSize];
            int[] sum12 = new int[varSize];
            boolean[] hasMissingGty = new boolean[varSize];
            int totalPedSubjectNum = pedEncodeGytIDMap.length;
            formatBialleleicUnphasedBits(runningList, totalPedSubjectNum, bits1, bits2, bits3, bits4, sum1, sum12, hasMissingGty);
            //  formatUnphasedGenotypesCalCorr(runningList, pedEncodeGytIDMap);
            double[] corrs = new double[varSize * (varSize - 1) / 2];
            int runningThread = 0;
            // int[] smallBlockIndexes = cc.partitionEvenBlock(size, threadNum > 50 ? threadNum : 50);
            int[] smallBlockIndexes = partitionEvenBlock(varSize, threadNum * 5);
            int smallBlockNum = smallBlockIndexes.length - 1;
            int lastButOneIndex = smallBlockIndexes[smallBlockNum - 1];
            int lastBlockLen = smallBlockIndexes[smallBlockNum] - smallBlockIndexes[smallBlockNum - 1];

            int blockLen = smallBlockIndexes[1] - smallBlockIndexes[0];
            int ti = 0;

            Long startTime = System.nanoTime();
            for (int i = 0; i < smallBlockNum; i++) {
                final double[][] rArray = new double[smallBlockNum - i][];
                runningThread = 0;
                for (int j = i; j < smallBlockNum; j++) {
                    //always prepare a rectangle
                    rArray[j - i] = new double[(smallBlockIndexes[i + 1] - smallBlockIndexes[i]) * (smallBlockIndexes[j + 1] - smallBlockIndexes[j])];
                    LDCalcTask task = new LDCalcTask(bits1, bits2, bits3, bits4, sum1, sum12, hasMissingGty, totalPedSubjectNum, smallBlockIndexes[i], smallBlockIndexes[i + 1],
                            smallBlockIndexes[j], smallBlockIndexes[j + 1], rArray[j - i], isPhased);
                    serv.submit(task);
                    runningThread++;
                }

                for (int index = 0; index < runningThread; index++) {
                    Future task1 = serv.take();
                    String infor1 = (String) task1.get();
                    // System.out.println(infor);
                }

                //unfortunately the last block often have uneven size 
                //correct for dependency
                double r;
                for (int t = smallBlockIndexes[i]; t < smallBlockIndexes[i + 1]; t++) {
                    for (int k = t + 1; k < varSize; k++) {
                        if (k < lastButOneIndex) {
                            r = rArray[(k - smallBlockIndexes[i]) / blockLen][(t - smallBlockIndexes[i]) * blockLen + (k - smallBlockIndexes[i]) % blockLen];
                        } else {
                            r = rArray[smallBlockNum - i - 1][(t - smallBlockIndexes[i]) * lastBlockLen + k - lastButOneIndex];
                        }
                        corrs[ti] = r;

                        // tmpFileWriter.write(allData.get(t).id + "\t" + allData.get(k).id + "\t" + r + "\n");
                        //System.out.print("(" + t + "," + k + ":" + r + ")");
                        //  System.out.println(ti + ":" + r);
                        ti++;
                    }
                    //  System.out.println();
                }
            }//end of  for (int t = 0; t < smallBlockNum; t++)
            System.out.println("LDCalcTask：" + (System.nanoTime() - startTime) / 1000000000.0);
            testingSize += 3000;
            runningList.clear();
        } while (testingSize < selectedVariantList.size());

        exec.shutdown();

    }

    public void addPatho(Chromosome chromosome, AnnotationSummarySet ass, HashMap<String, String[]> hmpPatho, int intPatho) {
        int intNum = 0;
        if (chromosome == null) {
            return;
        }
        //exclude the first column
        intPatho = intPatho - 1;
        int varFeatureNum = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            for (int i = 0; i < intPatho; i++) {
                var.setFeatureValue(varFeatureNum + i, null);
            }
        }
        String[] cells = null;
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb == null) {
                continue;
            }
            String strGene = var.geneSymb.toUpperCase();
            cells = hmpPatho.get(strGene);
            if (cells == null) {
                continue;
            }
            for (int i = 0; i < intPatho; i++) {
                var.setFeatureValue(varFeatureNum + i, hmpPatho.get(strGene)[i]);
            }
            intNum++;
        }

        ass.setAnnotNum(ass.getAnnotNum() + intNum);
    }

    public void annotateGenes(Chromosome chromosome, AnnotationSummarySet ass, HashMap<String, String> g2Phe) {
        int intNum = 0;
        if (chromosome == null) {
            return;
        }

        int varFeatureNum = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb == null) {
                var.setFeatureValue(varFeatureNum, null);
                continue;
            }
            String strGene = var.geneSymb.toUpperCase();
            if (g2Phe.containsKey(strGene)) {
                var.setFeatureValue(varFeatureNum, g2Phe.get(strGene));
                intNum++;
            } else {
                var.setFeatureValue(varFeatureNum, ".");
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + intNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - intNum);
    }

    public void annotateGenesArray(Chromosome chromosome, AnnotationSummarySet ass, Map<String, String[]> g2Phe, int valueNum) {
        int intNum = 0;
        if (chromosome == null) {
            return;
        }
        String[] values = null;
        int varFeatureNum = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb == null) {
                var.setFeatureValue(varFeatureNum, null);
                continue;
            }
            String strGene = var.geneSymb.toUpperCase();
            values = g2Phe.get(strGene);
            if (values != null) {
                for (int i = 0; i < valueNum; i++) {
                    var.setFeatureValue(varFeatureNum + i, values[i]);
                }
                intNum++;
            } else {
                for (int i = 0; i < valueNum; i++) {
                    var.setFeatureValue(varFeatureNum + i, ".");
                }
            }

        }

        ass.setAnnotNum(ass.getAnnotNum() + intNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - intNum);
    }

    public void hweTestVar(Chromosome chromosome, double hwdPCase, double hwdPControl, double hwdPAll, AnnotationSummarySet ass) {
        if (chromosome == null) {
            return;
        }

        // AffectedRefHomGtyNum\tAffectedHetGtyNum\tAffectedAltHomGtyNum\tUnaffectedRefHomGtyNum\tUnaffectedHetGtyNum\tUnaffectedAltHomGtyNum
        int[] caseGtys = new int[3];
        int[] controlGtys = new int[3];
        int excludeCase = 0, exludeControl = 0, excludeAll = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();
        for (Variant var : chromosome.variantList) {
            Arrays.fill(caseGtys, 0);
            caseGtys[0] = var.getAffectedRefHomGtyNum();
            caseGtys[1] = var.getAffectedHetGtyNum();
            caseGtys[2] = var.getAffectedAltHomGtyNum();

            Arrays.fill(controlGtys, 0);
            controlGtys[0] = var.getUnaffectedRefHomGtyNum();
            controlGtys[1] = var.getUnaffectedHetGtyNum();
            controlGtys[2] = var.getUnaffectedAltHomGtyNum();

            if (hwdPCase > 0) {
                if (hwdPCase > HardyWeinbergCalculation.hwCalculate(caseGtys[0], caseGtys[1], caseGtys[2])) {
                    excludeCase++;
                    continue;
                }
            }

            if (hwdPControl > 0) {
                if (hwdPControl > HardyWeinbergCalculation.hwCalculate(controlGtys[0], controlGtys[1], controlGtys[2])) {
                    exludeControl++;
                    continue;
                }
            }

            if (hwdPAll > 0) {
                if (hwdPAll > HardyWeinbergCalculation.hwCalculate(caseGtys[0] + controlGtys[0], caseGtys[1] + controlGtys[1], caseGtys[2] + controlGtys[2])) {
                    excludeAll++;
                    continue;
                }
            }
            tmpVarList.add(var);
        }
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarList);
        tmpVarList.clear();
        chromosome.buildVariantIndexMap();

        ass.setAnnotNum(ass.getAnnotNum() + excludeCase);//The annotNum is used to denote the excludeCase temporarily. 
        ass.setLeftNum(ass.getLeftNum() + exludeControl);//The leftNum is used to denote the excludeControl temporarily. 
        ass.setTotalNum(ass.getTotalNum() + excludeAll);//The excludeAll is used to denote the excludeAll temporarily. 

    }

    private class GeneFeatureAnnotTask extends Task implements Callable<String> {

        List<Variant> varList;
        int chrID;
        ReferenceGenome refGenome;
        Set<Byte> featureInSet;
        int featureNum;

        public GeneFeatureAnnotTask(List<Variant> varList, int chrID, ReferenceGenome refGenome, Set<Byte> featureInSet, int featureNum) {
            this.varList = varList;
            this.chrID = chrID;
            this.refGenome = refGenome;
            this.featureInSet = featureInSet;
            this.featureNum = featureNum;
        }

        @Override
        public String call() throws Exception {
            if (refGenome.getName().equals("refgene")) {
                if (varList == null) {
                    return null;
                }
                for (Variant var : varList) {
                    GeneFeature gf = refGenome.getVarFeature(STAND_CHROM_NAMES[chrID], var, true, featureInSet, new RNABoundaryIndex(0));
                    // Feature for UniProtFeature
                    if (gf.getInfor() == null) {
                        var.setFeatureValue(featureNum, ".");
                    } else {
                        var.setFeatureValue(featureNum, gf.getInfor());
                    }
                    var.setRefGeneAnnot(gf.getName());
                }
            } else if (refGenome.getName().equals("gencode")) {
                if (varList == null) {
                    return null;
                }
                for (Variant var : varList) {
                    GeneFeature gf = refGenome.getVarFeature(STAND_CHROM_NAMES[chrID], var, true, featureInSet, new RNABoundaryIndex(0));

                    // Feature for UniProtFeature
                    if (gf.getInfor() == null) {
                        var.setFeatureValue(featureNum, ".");
                    } else {
                        var.setFeatureValue(featureNum, gf.getInfor());
                    }
                    var.setgEncodeAnnot(gf.getName());
                }

            } else if (refGenome.getName().equals("knowngene")) {
                if (varList == null) {
                    return null;
                }
                for (Variant var : varList) {
                    GeneFeature gf = refGenome.getVarFeature(STAND_CHROM_NAMES[chrID], var, true, featureInSet, new RNABoundaryIndex(0));

                    // Feature for UniProtFeature
                    if (gf.getInfor() == null) {
                        var.setFeatureValue(featureNum, ".");
                    } else {
                        var.setFeatureValue(featureNum, gf.getInfor());
                    }
                    var.setKnownGeneAnnot(gf.getName());
                }

            } else if (refGenome.getName().equals("ensembl")) {
                if (varList == null) {
                    return null;
                }
                for (Variant var : varList) {
                    GeneFeature gf = refGenome.getVarFeature(STAND_CHROM_NAMES[chrID], var, true, featureInSet, new RNABoundaryIndex(0));

                    // Feature for UniProtFeature
                    if (gf.getInfor() == null) {
                        var.setFeatureValue(featureNum, ".");
                    } else {
                        var.setFeatureValue(featureNum, gf.getInfor());
                    }
                    var.setEnsemblGeneAnnot(gf.getName());
                }
            } else {
                //customized gene feature 
                if (varList == null) {
                    return null;
                }
                for (Variant var : varList) {
                    GeneFeature gf = refGenome.getVarFeature(STAND_CHROM_NAMES[chrID], var, true, featureInSet, new RNABoundaryIndex(0));

                    // Feature for UniProtFeature
                    if (gf.getInfor() == null) {
                        var.setFeatureValue(featureNum, ".");
                    } else {
                        var.setFeatureValue(featureNum, gf.getInfor());
                    }
                    var.setCustomizedGeneAnnot(gf.getName());
                }
            }
            return "";
        }
    }

    public int[] partitionEvenBlock(int threadNum, int startIndex, int endIndex) throws Exception {
        int totalSnpSize = endIndex - startIndex;
        int intervalLen = totalSnpSize / threadNum;
        int[] bigBlockIndexes = null;

        if (intervalLen == 0) {
            // no need to block
            bigBlockIndexes = new int[2];
            bigBlockIndexes[0] = startIndex;
            bigBlockIndexes[1] = endIndex;

        } else {
            bigBlockIndexes = new int[threadNum + 1];
            Arrays.fill(bigBlockIndexes, startIndex);
            for (int s = 1; s < threadNum; s++) {
                bigBlockIndexes[s] = startIndex + s * intervalLen;
            }
            bigBlockIndexes[threadNum] = endIndex;
        }

        return bigBlockIndexes;

    }

    public class NoncodingRandomForestGFCTask implements Callable<String> {

        Variant[] variantArray;
        MyRandomForest[] myRandomForestList;
        Map genicMap;
        Set<String> geneSymbSet;
        boolean needVerboseNoncode;
        int start;
        int end;
        int featureNum;

        public NoncodingRandomForestGFCTask(MyRandomForest[] myRandomForestList, Map genicMap, int featureNum, int start, int end, Variant[] variantArray, boolean needVerboseNoncode) {
            this.variantArray = variantArray;
            this.genicMap = genicMap;
            this.needVerboseNoncode = needVerboseNoncode;
            this.start = start;
            this.end = end;
            this.featureNum = featureNum;
            this.myRandomForestList = myRandomForestList;

        }

        @Override
        public String call() throws Exception {
            int funseq2Index = -1;
            for (int i = start; i < end; i++) {
                Variant var = variantArray[i];
                String mostImportantGeneFeature = org.cobi.kggseq.Constants.VAR_FEATURE_NAMES[var.smallestFeatureID];
                int index = (Integer) genicMap.get(mostImportantGeneFeature);
                double[] score = new double[var.scores2.length];
                if (funseq2Index < 0) {
                    funseq2Index = score.length - 3;
                }
                for (int j = 0; j < score.length; j++) {
                    if (j == funseq2Index) {
//return back the non-log socre of FunSeq2                        
                        score[j] = Math.log10(var.scores2[j]);
                    } else {
                        score[j] = var.scores2[j] + 0.0;
                    }
                }

                if (!needVerboseNoncode) {
                    //var.scores = null;
                    var.scores2 = new float[10];
                    for (int j = 0; j < var.scores2.length; j++) {
                        var.scores2[j] = (float) score[score.length - var.scores2.length + j];
                    }
                    //System.arraycopy(score, score.length - var.scores.length, var.scores, 0, var.scores.length);
                }
                double tmpP[];
                try {
                    tmpP = myRandomForestList[index].getClassifyDistribution(score);
                    if (tmpP[0] == 1) {
                        var.setFeatureValue(featureNum, "Y");
                        var.setFeatureValue(featureNum + 1, String.valueOf(tmpP[1]));
                    } else {
                        var.setFeatureValue(featureNum, "N");
                        var.setFeatureValue(featureNum + 1, String.valueOf(tmpP[1]));
                    }
                } catch (Exception e) {
                    // TODO Auto-generated catch block
                    System.out.println("Randomforest cannot retrieve the performance!!!");
                }
            }
            return "finish";
        }
    }

    public String[] noncodingRandomForestGFC(Chromosome chroms, boolean needProgressionIndicator,
            boolean filterNonDisMut, int needThreadNumber, Boolean[] isReigonList,
            double[] iniScore, String[] currentLineList, int[] fixedPosition,
            BufferedReader[] lineReaderList, int scoreIndexNum, MyRandomForest[][] myRandomForestList,
            Map genicMap, boolean needVerboseNoncode,
            int featureNum, int chromIndex, FiltrationSummarySet dbNoncodePred91) throws Exception {
        String zeroLabel = ".";
        String oneLabel = "Y";
        String missingLabel = "?";
        ExecutorService executor = Executors.newFixedThreadPool(needThreadNumber);
        final CompletionService<String> serv = new ExecutorCompletionService<String>(executor);
        int varLineBufferCounter;
        int bufferSize = 10000 * needThreadNumber;
        long lineCounter = 0;
        // int bufferSize = 5;
        Variant[] parseVariantArray = new Variant[bufferSize];
        String cells[];
        int localIndex = 0;
        boolean loopAgain = true;
        int fileLineOutNum = 0;
        double[] score = null;
        varLineBufferCounter = 0;
        List<Variant> tempVarList = new ArrayList<Variant>();

        int totalVarNum = 0;
        int counterDis = 0;
        int counterCon = 0;

        boolean continueNextLoci;
        int localPosition;
        boolean hasReadNewRow[] = new boolean[isReigonList.length];
        Arrays.fill(hasReadNewRow, true);

        String[] currentChrom = new String[isReigonList.length];
        int[] currentStartPos = new int[isReigonList.length];
        int[] currentEndPos = new int[isReigonList.length];

        int localListSize = chroms.variantList.size();
        int startPosition, endPosition;

        if (needProgressionIndicator) {
            //   System.out.print("Parsing Chrom" + STAND_CHROM_NAMES[chromIndex] + ", ");
        }

        String[][] cellsTmpStorageList = new String[isReigonList.length][];
        boolean initializeList[] = new boolean[isReigonList.length];
        Arrays.fill(initializeList, false);

        while (loopAgain) {
            lineCounter++;
            if (needProgressionIndicator && lineCounter % 100000 == 0) {
                String prog = String.valueOf(lineCounter);
                System.out.print(prog);
                char[] backSpaces = new char[prog.length()];
                //  Arrays.fill(backSpaces, '\b');
                //  System.out.print(backSpaces);
            }
            if (fileLineOutNum == isReigonList.length) {
                break;
            }
            continueNextLoci = true;
            if (localListSize <= localIndex) {
                break;
            }

            Variant var = chroms.variantList.get(localIndex);
            localPosition = var.refStartPosition;

            for (int j = 0; j < isReigonList.length; j++) {
                if (hasReadNewRow[j]) {
                    if (isReigonList[j]) {
                        if (currentLineList[j] == null) {
                            continue;
                        }
                        cells = Util.tokenize(currentLineList[j], '\t', 2);
                        currentChrom[j] = cells[0];
                        currentStartPos[j] = Util.parseInt(cells[1]);
                        currentEndPos[j] = Util.parseInt(cells[2]);
                    } else {
                        if (currentLineList[j] == null) {
                            continue;
                        }
                        cells = Util.tokenize(currentLineList[j], '\t', 1);
                        currentChrom[j] = cells[0];
                        currentStartPos[j] = Util.parseInt(cells[1]);
                    }
                    if (!initializeList[j]) {
                        cells = Util.tokenize(currentLineList[j], '\t');
                        cellsTmpStorageList[j] = new String[cells.length];
                        initializeList[j] = true;
                    }
                }

                if (!STAND_CHROM_NAMES[chromIndex].equals(currentChrom[j])) {
                    int anchorIndex = 0;
                    for (int j2 = 0; j2 < STAND_CHROM_NAMES.length; j2++) {
                        if (STAND_CHROM_NAMES[j2].equals(currentChrom[j])) {
                            anchorIndex = j2;
                            break;
                        }
                    }

                    if (anchorIndex > chromIndex) {
                        continue;
                    } else if (anchorIndex < chromIndex) {
                        continueNextLoci = false;
                        if ((currentLineList[j] = lineReaderList[j].readLine()) == null) {
                            fileLineOutNum += 1;
                        }
                        hasReadNewRow[j] = true;
                        continue;
                    }
                }
                if (isReigonList[j]) {
                    startPosition = currentStartPos[j];
                    endPosition = currentEndPos[j];
                    if (localPosition >= endPosition) {
                        if ((currentLineList[j] = lineReaderList[j].readLine()) == null) {
                            fileLineOutNum += 1;
                        }
                        continueNextLoci = false;
                        hasReadNewRow[j] = true;
                    } else if (localPosition >= startPosition) {
                        hasReadNewRow[j] = false;
                        if (score == null) {
                            score = new double[scoreIndexNum];
                            System.arraycopy(iniScore, 0, score, 0, iniScore.length);
                        }
//                        cells = Util.tokenize(currentLineList[j], '\t');
                        tokenize(currentLineList[j], '\t', cellsTmpStorageList[j]);
                        cells = cellsTmpStorageList[j];
                        for (int k = 3; k < cells.length; k++) {
                            if (cells[k].equals(oneLabel)) {
                                score[k - 3 + fixedPosition[j]] = 1.0;
                            } else if (cells[k].equals(zeroLabel)) {
                                score[k - 3 + fixedPosition[j]] = 0.0;
                            } else {
                                score[k - 3 + fixedPosition[j]] = Util.parseFloat(cells[k]);
                            }
                        }
                    }
                } else {
                    int position = currentStartPos[j];
                    if (position < localPosition) {
                        continueNextLoci = false;
                        if ((currentLineList[j] = lineReaderList[j].readLine()) == null) {
                            fileLineOutNum += 1;
                        }
                        hasReadNewRow[j] = true;
                    } else if (position == localPosition) {
                        if (score == null) {
                            score = new double[scoreIndexNum];
                            System.arraycopy(iniScore, 0, score, 0, iniScore.length);
                        }
//                        cells = Util.tokenize(currentLineList[j], '\t');

                        tokenize(currentLineList[j], '\t', cellsTmpStorageList[j]);
                        cells = cellsTmpStorageList[j];
                        for (int k = 2; k < cells.length; k++) {
                            if (cells[k].equals(oneLabel)) {
                                score[k - 2 + fixedPosition[j]] = 1.0;
                            } else if (cells[k].equals(zeroLabel)) {
                                score[k - 2 + fixedPosition[j]] = 0.0;
                            } else if (cells[k].equals(missingLabel)) {
                                score[k - 2 + fixedPosition[j]] = Double.NaN;
                            } else {
                                score[k - 2 + fixedPosition[j]] = Util.parseFloat(cells[k]);
                            }
                        }
                        if ((currentLineList[j] = lineReaderList[j].readLine()) == null) {
                            fileLineOutNum += 1;
                        }
                        hasReadNewRow[j] = true;
                    } else {
                        hasReadNewRow[j] = false;
                    }
                }
            }
            if (continueNextLoci || fileLineOutNum == isReigonList.length) {
                localIndex++;
                if (score != null) {
                    String mostImportantGeneFeature = org.cobi.kggseq.Constants.VAR_FEATURE_NAMES[var.smallestFeatureID];
                    if (genicMap.containsKey(mostImportantGeneFeature)) {
                        var.scores2 = new float[score.length];
                        for (int k = 0; k < score.length; k++) {
                            var.scores2[k] = (float) score[k];
                        }
                        parseVariantArray[varLineBufferCounter++] = var;
                        if (varLineBufferCounter >= bufferSize) {
                            int[] blocks = partitionEvenBlock(needThreadNumber, 0, varLineBufferCounter);
                            int blockNum = blocks.length - 1;
                            int runningThread = 0;
                            for (int s = 0; s < blockNum; s++) {
                                NoncodingRandomForestGFCTask task = new NoncodingRandomForestGFCTask(myRandomForestList[s], genicMap, featureNum, blocks[s], blocks[s + 1], parseVariantArray, needVerboseNoncode);
                                serv.submit(task);
                                runningThread++;
                            }
                            for (int s = 0; s < runningThread; s++) {
                                Future<String> task = serv.take();
                                String infor = task.get();
                                // System.out.println(infor);
                            }

                            for (int s = 0; s < varLineBufferCounter; s++) {
                                var = parseVariantArray[s];
                                if (var.featureValues[featureNum] == null) {
                                    if (filterNonDisMut) {
                                        tempVarList.add(var);
                                    }
                                } else if (var.featureValues[featureNum].equals("Y")) {
                                    if (filterNonDisMut) {
                                        tempVarList.add(var);
                                    }
                                    counterDis += 1;
                                } else {
                                    counterCon += 1;
                                }
                            }
                            varLineBufferCounter = 0;
                        }
                    } else {
                        tempVarList.add(var);
                        var.setFeatureValue(featureNum, null);
                        var.setFeatureValue(featureNum + 1, null);
                    }
                    score = null;
                } else {
                    tempVarList.add(var);
                    var.setFeatureValue(featureNum, null);
                    var.setFeatureValue(featureNum + 1, null);
                }
            }
        }

        if (varLineBufferCounter > 0) {
            int[] blocks = partitionEvenBlock(needThreadNumber, 0, varLineBufferCounter);
            int blockNum = blocks.length - 1;
            int runningThread = 0;
            for (int s = 0; s < blockNum; s++) {
                serv.submit(new NoncodingRandomForestGFCTask(myRandomForestList[s], genicMap, featureNum, blocks[s], blocks[s + 1], parseVariantArray, needVerboseNoncode));
                runningThread++;
            }

            for (int s = 0; s < runningThread; s++) {
                Future task = serv.take();
                String infor = (String) task.get();
                // System.out.println(infor);
            }
            for (int s = 0; s < varLineBufferCounter; s++) {
                Variant var = parseVariantArray[s];
                if (var.featureValues[featureNum] == null) {
                    if (filterNonDisMut) {
                        tempVarList.add(var);
                    }
                } else if (var.featureValues[featureNum].equals("Y")) {
                    if (filterNonDisMut) {
                        tempVarList.add(var);
                    }
                    counterDis += 1;
                } else {
                    counterCon += 1;
                }
            }
            varLineBufferCounter = 0;
        }

        if (filterNonDisMut) {
            chroms.variantList.clear();
            chroms.setHasNotOrderVariantList(true);
            chroms.variantList.addAll(tempVarList);
            totalVarNum += tempVarList.size();
            tempVarList.clear();
            chroms.buildVariantIndexMap();
        }

        dbNoncodePred91.increaseCount(0, counterDis);
        dbNoncodePred91.increaseCount(1, counterCon);
        dbNoncodePred91.increaseCount(2, totalVarNum);

        executor.shutdown();

        if (needProgressionIndicator) {
            //  System.out.println("finished");
        }

        return currentLineList;
//			genome.buildVariantIndexMapOnChromosomes();
//			StringBuilder info = new StringBuilder("The number of predicted disease-causal and non-disease-causal variants are " + counterDis + "(#gene " + geneSymbSet.size()
//				+ ") and " + counterCon + " among " + genome.getVarNum() + " variants on the genome according to the Random Forests prediction model");
        // info.append(" trained by ExoVar dataset (http://grass.cgs.hku.hk/limx/kggseq/download/ExoVar.xls).");
        // info.append(" trained by COSMIC dataset (http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/).");
//			LOG.info(info);

//			if (filterNonDisMut) {
//				genome.setVarNum(totalVarNum);
//				info.delete(0, info.length());
//				info.append(totalVarNum).append(" variant(s) are left after filtered by the disease mutation prediction.");
//				LOG.info(info);
//			}
    }

    public String[] noncodingRandomForestGFRegions(Chromosome chroms, boolean needProgressionIndicator,
            boolean filterNonDisMut, int needThreadNumber, Boolean[] isReigonList,
            double[] iniScore, String[] currentLineList, int[] fixedPosition,
            BufferedReader[] lineReaderList, int scoreIndexNum, MyRandomForest[][] myRandomForestList,
            Map genicMap, boolean needVerboseNoncode,
            int featureNum, int chromIndex, FiltrationSummarySet dbNoncodePred91) throws Exception {
        String zeroLabel = ".";
        String oneLabel = "Y";
        String missingLabel = "?";
        ExecutorService executor = Executors.newFixedThreadPool(needThreadNumber);
        final CompletionService<String> serv = new ExecutorCompletionService<String>(executor);
        int varLineBufferCounter;
        int bufferSize = 10000 * needThreadNumber;
        long lineCounter = 0;
        // int bufferSize = 5;
        Variant[] parseVariantArray = new Variant[bufferSize];
        String cells[];
        int localIndex = 0;
        boolean loopAgain = true;
        int fileLineOutNum = 0;
        double[] score = null;
        varLineBufferCounter = 0;
        List<Variant> tempVarList = new ArrayList<Variant>();

        int totalVarNum = 0;
        int counterDis = 0;
        int counterCon = 0;

        boolean continueNextLoci;
        int localPosition;
        boolean hasReadNewRow[] = new boolean[isReigonList.length];
        Arrays.fill(hasReadNewRow, true);

        String[] currentChrom = new String[isReigonList.length];
        int[] currentStartPos = new int[isReigonList.length];
        int[] currentEndPos = new int[isReigonList.length];

        int localListSize = chroms.variantList.size();
        int startPosition, endPosition;

        String[][] cellsTmpStorageList = new String[isReigonList.length][];
        boolean initializeList[] = new boolean[isReigonList.length];
        Arrays.fill(initializeList, false);
        boolean allNA = true;
        while (loopAgain) {
            lineCounter++;
            if (fileLineOutNum == isReigonList.length) {
                break;
            }
            continueNextLoci = true;
            if (localListSize <= localIndex) {
                break;
            }

            Variant var = chroms.variantList.get(localIndex);
            localPosition = var.refStartPosition;

            for (int j = 0; j < isReigonList.length; j++) {
                if (hasReadNewRow[j]) {
                    if (isReigonList[j]) {
                        if (currentLineList[j] == null) {
                            continue;
                        }
                        cells = Util.tokenize(currentLineList[j], '\t', 2);
                        currentChrom[j] = cells[0];
                        currentStartPos[j] = Util.parseInt(cells[1]);
                        currentEndPos[j] = Util.parseInt(cells[2]);
                    } else {
                        if (currentLineList[j] == null) {
                            continue;
                        }
                        cells = Util.tokenize(currentLineList[j], '\t', 2);
                        currentChrom[j] = cells[0];
                        currentStartPos[j] = Util.parseInt(cells[2]);
                    }
                    if (!initializeList[j]) {
                        cells = Util.tokenize(currentLineList[j], '\t');
                        cellsTmpStorageList[j] = new String[cells.length];
                        initializeList[j] = true;
                    }
                }

                if (!STAND_CHROM_NAMES[chromIndex].equals(currentChrom[j])) {
                    int anchorIndex = 0;
                    for (int j2 = 0; j2 < STAND_CHROM_NAMES.length; j2++) {
                        if (STAND_CHROM_NAMES[j2].equals(currentChrom[j])) {
                            anchorIndex = j2;
                            break;
                        }
                    }

                    if (anchorIndex > chromIndex) {
                        continue;
                    } else if (anchorIndex < chromIndex) {
                        continueNextLoci = false;
                        if ((currentLineList[j] = lineReaderList[j].readLine()) == null) {
                            fileLineOutNum += 1;
                        }
                        hasReadNewRow[j] = true;
                        continue;
                    }
                }
                //only consider regions
                if (isReigonList[j]) {
                    startPosition = currentStartPos[j];
                    endPosition = currentEndPos[j];
                    if (localPosition >= endPosition) {
                        if ((currentLineList[j] = lineReaderList[j].readLine()) == null) {
                            fileLineOutNum += 1;
                        }
                        continueNextLoci = false;
                        hasReadNewRow[j] = true;
                    } else if (localPosition >= startPosition) {
                        hasReadNewRow[j] = false;
                        if (score == null) {
                            score = new double[scoreIndexNum];
                            System.arraycopy(iniScore, 0, score, 0, iniScore.length);
                        }
//                        cells = Util.tokenize(currentLineList[j], '\t');
                        tokenize(currentLineList[j], '\t', cellsTmpStorageList[j]);
                        cells = cellsTmpStorageList[j];
                        for (int k = 3; k < cells.length; k++) {
                            if (cells[k].equals(oneLabel)) {
                                score[k - 3 + fixedPosition[j]] = 1.0;
                            } else if (cells[k].equals(zeroLabel)) {
                                score[k - 3 + fixedPosition[j]] = 0.0;
                            } else {
                                score[k - 3 + fixedPosition[j]] = Util.parseFloat(cells[k]);
                            }
                        }
                    }
                }
            }

            String mostImportantGeneFeature = org.cobi.kggseq.Constants.VAR_FEATURE_NAMES[var.smallestFeatureID];
            if (continueNextLoci || fileLineOutNum == isReigonList.length) {
                localIndex++;
                if (genicMap.containsKey(mostImportantGeneFeature)) {
                    if (score != null) {
                        float[] tmpScores = Arrays.copyOf(var.scores2, var.scores2.length);
                        var.scores2 = new float[score.length + tmpScores.length];
                        for (int k = 0; k < score.length; k++) {
                            var.scores2[k] = (float) score[k];
                        }
                        System.arraycopy(tmpScores, 0, var.scores2, score.length, tmpScores.length);
                    } else {
                        allNA = true;
                        for (int i = 0; i < var.scores2.length; i++) {
                            if (!Float.isNaN(var.scores2[i])) {
                                allNA = false;
                                break;
                            }
                        }
                        float[] tmpScores = Arrays.copyOf(var.scores2, var.scores2.length);
                        var.scores2 = new float[scoreIndexNum + tmpScores.length];
                        //Arrays.fill(var.scores2, Float.NaN);
                        Arrays.fill(var.scores2, 0);
                        System.arraycopy(tmpScores, 0, var.scores2, scoreIndexNum, tmpScores.length);

                        if (allNA) {
                            tempVarList.add(var);
                            var.setFeatureValue(featureNum, null);
                            var.setFeatureValue(featureNum + 1, null);

                            /*
                            if (needVerboseNoncode) {
                                for (int i = 0; i < scoreIndexNum; i++) {
                                    var.scores2[i] = 0;
                                }
                                for (int i = scoreIndexNum; i < scoreIndexNum + 10; i++) {
                                    var.scores2[i] = Float.NaN;
                                }
                            } 
                             */
                            continue;
                        }
                    }
                    parseVariantArray[varLineBufferCounter++] = var;
                    if (varLineBufferCounter >= bufferSize) {
                        int[] blocks = partitionEvenBlock(needThreadNumber, 0, varLineBufferCounter);
                        int blockNum = blocks.length - 1;
                        int runningThread = 0;
                        for (int s = 0; s < blockNum; s++) {
                            NoncodingRandomForestGFCTask task = new NoncodingRandomForestGFCTask(myRandomForestList[s], genicMap, featureNum, blocks[s], blocks[s + 1], parseVariantArray, needVerboseNoncode);
                            serv.submit(task);
                            runningThread++;
                        }
                        for (int s = 0; s < runningThread; s++) {
                            Future<String> task = serv.take();
                            String infor = task.get();
                            // System.out.println(infor);
                        }

                        for (int s = 0; s < varLineBufferCounter; s++) {
                            var = parseVariantArray[s];
                            if (var.featureValues[featureNum] == null) {
                                if (filterNonDisMut) {
                                    tempVarList.add(var);
                                }
                            } else if (var.featureValues[featureNum].equals("Y")) {
                                if (filterNonDisMut) {
                                    tempVarList.add(var);
                                }
                                counterDis += 1;
                            } else {
                                counterCon += 1;
                            }

                            //remove the epigenomic scores from scores2
                            if (!needVerboseNoncode && var.scores2.length > 10) {
                                float[] tmpScores = Arrays.copyOf(var.scores2, 10);
                                var.scores2 = new float[scoreIndexNum + tmpScores.length];
                                //Arrays.fill(var.scores2, Float.NaN);
                                Arrays.fill(var.scores2, 0);
                                System.arraycopy(tmpScores, 0, var.scores2, scoreIndexNum, tmpScores.length);
                            }
                        }
                        varLineBufferCounter = 0;
                    }
                } else {
                    tempVarList.add(var);
                    var.setFeatureValue(featureNum, null);
                    var.setFeatureValue(featureNum + 1, null);

                    var.scores2 = new float[scoreIndexNum + var.scores2.length];
                    Arrays.fill(var.scores2, Float.NaN);

                }
                score = null;
            }
        }

        if (varLineBufferCounter > 0) {
            int[] blocks = partitionEvenBlock(needThreadNumber, 0, varLineBufferCounter);
            int blockNum = blocks.length - 1;
            int runningThread = 0;
            for (int s = 0; s < blockNum; s++) {
                serv.submit(new NoncodingRandomForestGFCTask(myRandomForestList[s], genicMap, featureNum, blocks[s], blocks[s + 1], parseVariantArray, needVerboseNoncode));
                runningThread++;
            }

            for (int s = 0; s < runningThread; s++) {
                Future task = serv.take();
                String infor = (String) task.get();
                // System.out.println(infor);
            }
            for (int s = 0; s < varLineBufferCounter; s++) {
                Variant var = parseVariantArray[s];
                if (var.featureValues[featureNum] == null) {
                    if (filterNonDisMut) {
                        tempVarList.add(var);
                    }
                } else if (var.featureValues[featureNum].equals("Y")) {
                    if (filterNonDisMut) {
                        tempVarList.add(var);
                    }
                    counterDis += 1;
                } else {
                    counterCon += 1;
                }

                //remove the epigenomic scores from scores2
                if (!needVerboseNoncode && var.scores2.length > 10) {
                    float[] tmpScores = Arrays.copyOf(var.scores2, 10);
                    var.scores2 = new float[scoreIndexNum + tmpScores.length];
                    //Arrays.fill(var.scores2, Float.NaN);
                    Arrays.fill(var.scores2, 0);
                    System.arraycopy(tmpScores, 0, var.scores2, scoreIndexNum, tmpScores.length);
                }

            }

            varLineBufferCounter = 0;
        }

        if (filterNonDisMut) {
            chroms.variantList.clear();
            chroms.setHasNotOrderVariantList(true);
            chroms.variantList.addAll(tempVarList);
            totalVarNum += tempVarList.size();
            tempVarList.clear();
            chroms.buildVariantIndexMap();
        }

        dbNoncodePred91.increaseCount(0, counterDis);
        dbNoncodePred91.increaseCount(1, counterCon);
        dbNoncodePred91.increaseCount(2, totalVarNum);

        executor.shutdown();

        return currentLineList;
    }

    private class NonCodingPredicLJTask extends Task implements Callable<String> {

        List<Variant> varList;
        Bayes bayesPredictor;
        int startFeatureNum;
        int[] effeIndexes;
        double[] baye1NonCodingPredic;
        double[] baye2NonCodingPredic;
        String chromName;
        boolean doCellType;

        public NonCodingPredicLJTask(List<Variant> varList, Bayes bayesPredictor, int startFeatureNum, int[] effeIndexes,
                double[] baye1NonCodingPredic, double[] baye2NonCodingPredic, String chromName, boolean docellType) {
            this.varList = varList;
            this.bayesPredictor = bayesPredictor;
            this.startFeatureNum = startFeatureNum;
            this.effeIndexes = effeIndexes;
            this.baye1NonCodingPredic = baye1NonCodingPredic;
            this.baye2NonCodingPredic = baye2NonCodingPredic;
            this.chromName = chromName;
            this.doCellType = docellType;
        }

        @Override
        public String call() throws Exception {
            int localListSize = varList.size();
            float cell_p;
            for (int i = 0; i < localListSize; i++) {
                Variant var = varList.get(i);

                double[] scores = bayesPredictor.getBayesScoreCompsit(var.scores2, effeIndexes, baye1NonCodingPredic, baye2NonCodingPredic);
                if (Double.isNaN(scores[0])) {
                    var.setFeatureValue(startFeatureNum, ".");
                } else {
                    var.setFeatureValue(startFeatureNum, String.valueOf(scores[0]));
                }
                //cell-type specific
                if (doCellType) {
                    cell_p = bayesPredictor.getCellSpecificScore(chromName, var.refStartPosition);
                    //sb.append(cell_p + "\t" + composite_p * cell_p / 0.5); 
                    if (Double.isNaN(scores[1])) {
                        var.setFeatureValue(startFeatureNum + 1, ".");
                        var.setFeatureValue(startFeatureNum + 2, ".");
                        var.setFeatureValue(startFeatureNum + 3, ".");
                    } else {
                        var.setFeatureValue(startFeatureNum + 1, String.valueOf(scores[1]));
                        if (Float.isNaN(cell_p)) {
                            var.setFeatureValue(startFeatureNum + 2, ".");
                            var.setFeatureValue(startFeatureNum + 3, ".");
                        } else {
                            var.setFeatureValue(startFeatureNum + 2, String.valueOf(cell_p));
                            var.setFeatureValue(startFeatureNum + 3, String.valueOf(cell_p * scores[1] / 0.5));
                        }
                    }

                } else if (Double.isNaN(scores[1])) {
                    var.setFeatureValue(startFeatureNum + 1, ".");
                } else {
                    var.setFeatureValue(startFeatureNum + 1, String.valueOf(scores[1]));
                }
            }
            return "";
        }
    }

    public void noncodingCompositPredictionMultiThread(Chromosome chrom, int[] effeIndexes,
            Bayes bayesPredictor, int startFeatureNum, int threadNum, double[] baye1NonCodingPredic,
            double[] baye2NonCodingPredic, boolean calCellType) {
        List<Variant> varList = chrom.variantList;

        if (varList.isEmpty()) {
            return;
        }

        String currentLine = null;
        try {
            // System.out.println(" Chromosome " + chromosome.getName() + " " + rsFile.getName());
            // Options.REF_CHROM_NAMES[t]);
            int varSize = varList.size();
            int[] bounderies = partitionEvenBlock(varSize, threadNum);
            int actualThreadNum = bounderies.length - 1;

            ExecutorService exec = Executors.newFixedThreadPool(actualThreadNum);
            final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
            int runningThread = 0;

            for (int f = 1; f < bounderies.length; f++) {
                NonCodingPredicLJTask varTaks = new NonCodingPredicLJTask(varList.subList(bounderies[f - 1], bounderies[f]), bayesPredictor, startFeatureNum, effeIndexes,
                        baye1NonCodingPredic, baye2NonCodingPredic, chrom.getName(), calCellType) {
                    @Override
                    protected void fireTaskComplete() {

                    }
                };
                serv.submit(varTaks);
                runningThread++;
            }

            for (int s = 0; s < runningThread; s++) {
                Future<String> task = serv.take();
                String infor = task.get();
                // System.out.println(infor);
            }
            exec.shutdown();

        } catch (Exception ex) {
            if (currentLine != null) {
                System.err.println("Errors in a row: " + currentLine);
            }
            ex.printStackTrace();
        }
    }

    class RiskPredictionRandomForestTask extends Task implements Callable<String> {

        MyRandomForest myRandomForest;
        List<Variant> varList;
        int featureNum;

        public RiskPredictionRandomForestTask(MyRandomForest myRandomForest, List<Variant> varList, int featureNum) {
            this.myRandomForest = myRandomForest;
            this.varList = varList;
            this.featureNum = featureNum;
        }

        @Override
        public String call() throws Exception {
            double tmpP[];
            String[] predictors = myRandomForest.getName().split("-");
            double[] scores = new double[predictors.length];
            boolean allMissing = true;
            for (Variant var : varList) {
                //ingore intergenic variants
                if (var.smallestFeatureID >= 15) {
                    var.setFeatureValue(featureNum, null);
                    var.setFeatureValue(featureNum + 1, null);
                    continue;
                }
                if (var.scores1 != null) {
                    Arrays.fill(scores, Double.NaN);
                    allMissing = true;
                    for (int j = 0; j < scores.length; j++) {
                        if (!Float.isNaN(var.scores1[j])) {
                            scores[j] = var.scores1[j];
                            allMissing = false;
                        }
                    }
                    if (allMissing) {
                        var.setFeatureValue(featureNum, null);
                        var.setFeatureValue(featureNum + 1, null);
                        continue;
                    }
                    tmpP = myRandomForest.getClassifyDistribution(scores);

                    if (tmpP[0] == 1) {
                        var.setFeatureValue(featureNum, "Y");
                        var.setFeatureValue(featureNum + 1, String.valueOf(tmpP[1]));
                    } else {
                        var.setFeatureValue(featureNum, "N");
                        var.setFeatureValue(featureNum + 1, String.valueOf(tmpP[1]));
                    }
                } else {
                    // keep varaint have no risk scores which may be safer                
                    var.setFeatureValue(featureNum, null);
                    var.setFeatureValue(featureNum + 1, null);
                }
            }
            return "";
        }
    }

    public void riskPredictionRandomForest(Chromosome chromosome, MyRandomForest[] myRandomForest, int maxThreadNum,
            boolean filterNonDisMut, FiltrationSummarySet dbNSFPMendelPred) throws Exception {
        int totalVarNum = 0;
        int counterDis = 0;
        int counterCon = 0;

        // zhicheng
        // System.out.println(myRandomForest.getName());
        List<Variant> keptVarList = new ArrayList<Variant>();
        int featureNum = dbNSFPMendelPred.getAvailableFeatureIndex();
        if (chromosome == null || chromosome.variantList.isEmpty()) {
            return;
        }

        ExecutorService exec = Executors.newFixedThreadPool(maxThreadNum);
        final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
        int runningThread = 0;
        int varNum = chromosome.variantList.size();

        int[] blocks = org.cobi.util.thread.Util.partitionEvenBlock(maxThreadNum, 0, varNum);
        int blockNum = blocks.length - 1;
        runningThread = 0;
        for (int s = 0; s < blockNum; s++) {
            RiskPredictionRandomForestTask task = new RiskPredictionRandomForestTask(myRandomForest[s], chromosome.variantList.subList(blocks[s], blocks[s + 1]), featureNum);
            serv.submit(task);
            runningThread++;
        }

        for (int s = 0; s < runningThread; s++) {
            Future<String> task = serv.take();
            String infor = task.get();
            //  System.out.println(infor);
        }

        exec.shutdown();

        Set<String> geneSymbSet = new HashSet<String>();
        for (Variant var : chromosome.variantList) {
            if (var.featureValues[featureNum] == null) {
                if (filterNonDisMut) {
                    keptVarList.add(var);
                }
                continue;
            } else if (var.featureValues[featureNum].equals("Y")) {
                if (filterNonDisMut) {
                    keptVarList.add(var);
                }
                if (var.geneSymb != null) {
                    geneSymbSet.add(var.geneSymb);
                }
                counterDis += 1;
            } else {
                counterCon += 1;
                if (var.smallestFeatureID != 6) {
                    keptVarList.add(var);
                }
            }
        }

        if (filterNonDisMut) {
            chromosome.variantList.clear();
            chromosome.variantList.addAll(keptVarList);
            totalVarNum += keptVarList.size();
            keptVarList.clear();
        }
        chromosome.buildVariantIndexMap();
        dbNSFPMendelPred.increaseCount(0, counterDis);
        dbNSFPMendelPred.increaseCount(1, geneSymbSet.size());
        dbNSFPMendelPred.increaseCount(2, counterCon);
        dbNSFPMendelPred.increaseCount(3, totalVarNum);
    }

    public void matchTrioSet(List<Individual> subjectIDList, List<int[]> triosIDList) throws Exception {
        for (int i = 0; i < subjectIDList.size(); i++) {
            Individual indiv0 = subjectIDList.get(i);
            if (indiv0.getDadID().equals("0") && indiv0.getMomID().equals("0")) {
                continue;
            }

            int[] setIDs = new int[3];
            setIDs[0] = i;
            setIDs[1] = -9;
            setIDs[2] = -9;

            for (int t = 0; t < subjectIDList.size(); t++) {
                // child father mother ids
                if (!subjectIDList.get(t).getFamilyID().equals(indiv0.getFamilyID())) {
                    continue;
                }
                if (subjectIDList.get(t).getIndividualID().equals(indiv0.getDadID())) {
                    setIDs[1] = t;
                } else if (subjectIDList.get(t).getIndividualID().equals(indiv0.getMomID())) {
                    setIDs[2] = t;
                }
                if (setIDs[0] != -9 && setIDs[1] != -9 && setIDs[2] != -9) {
                    break;
                }
            }

            // at least one parent is needed
            if (setIDs[0] != -9 && (setIDs[1] != -9 || setIDs[2] != -9)) {
                triosIDList.add(setIDs);
            }
        }
    }

    public void summarizeSomatNSVarPerGene(Chromosome chromosome, int somatNumIndex, int readInfoIndex, Set<Byte> dependentGeneFeature, Set<Byte> independentGeneFeature) throws Exception {
        List<Gene> geneList = chromosome.geneList;
        String evenInfor = null;
        int somatEvents = 0;
        Map<String, double[]> geneMutNum = new HashMap<String, double[]>();
        Map<String, StringBuilder> geneReadsTestInfo = new HashMap<String, StringBuilder>();
        Set<String> availableGeneSet = new HashSet<String>();
        String[] somatEventStr = null;
        String geneFeatureAnnot;
        StringBuilder sb = new StringBuilder();
        Set<String> geneSet = new HashSet<String>();
        int index;
        String geneSymb = null, feature = null;
        Byte featureID;
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb == null) {
                continue;
            }
            if (somatNumIndex >= 0) {
                evenInfor = var.getFeatureValues()[somatNumIndex];
                if (evenInfor == null) {
                    somatEvents = 0;
                } else {
                    somatEventStr = Util.tokenize(evenInfor, ' ');
                    if (somatEventStr != null && somatEventStr.length > 0) {
                        somatEvents = (Util.parseInt(somatEventStr[0]));
                    }
                }
            }

            geneFeatureAnnot = var.getRefGeneAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            geneFeatureAnnot = var.getgEncodeAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            geneFeatureAnnot = var.getKnownGeneAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            geneFeatureAnnot = var.getEnsemblGeneAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }

            if (sb.length() > 0) {
                String[] items = Util.tokenize(sb.substring(0, sb.length() - 1), ';');
                for (String item : items) {
                    index = item.indexOf(':');
                    if (index > 0) {
                        geneSymb = item.substring(0, index);
                    }
                    if (geneSet.contains(geneSymb)) {
                        continue;
                    }
                    geneSet.add(geneSymb);

                    index = item.lastIndexOf(':');
                    if (index > 0) {
                        feature = item.substring(index + 1);
                    }
                    featureID = geneFeatureMap.get(feature);
                    double[] mutNums = geneMutNum.get(geneSymb);
                    if (mutNums == null) {
                        mutNums = new double[4];
                        Arrays.fill(mutNums, 0);
                        geneMutNum.put(geneSymb, mutNums);
                    }
                    StringBuilder tesInfo = geneReadsTestInfo.get(geneSymb);
                    if (tesInfo == null) {
                        tesInfo = new StringBuilder();
                        geneReadsTestInfo.put(geneSymb, tesInfo);
                    }

                    // non-synonemous preidcted as Y, non-synonemous preidcted
                    // as N, synonemous
                    if (dependentGeneFeature.contains(featureID)) {
                        // only do it for missense variants; many stop loss
                        // variants have no proection at loss of function
                        // variants
                        mutNums[0] += somatEvents;
                    } else if (independentGeneFeature.contains(featureID)) {
                        mutNums[1] += somatEvents;
                    }

                    if (readInfoIndex >= 0) {
                        String readInfo = var.getFeatureValues()[readInfoIndex];
                        if (readInfo != null) {
                            tesInfo.append(readInfo);
                            tesInfo.append('|');
                            String[] readInfos = readInfo.split("[|]");
                            for (String readS : readInfos) {
                                index = readS.lastIndexOf(':');
                                if (index < 0) {
                                    continue;
                                }
                                // System.out.println(readS);
                                double or = Double.parseDouble(readS.substring(index + 1));
                                // or = Math.log(or);
                                if (featureID <= 6) {
                                    mutNums[2] += or;
                                } else if (featureID == 7) {
                                    mutNums[3] += or;
                                }
                            }
                        }
                    }

                    if (!availableGeneSet.contains(var.geneSymb)) {
                        geneList.add(new Gene(var.geneSymb));
                        availableGeneSet.add(var.geneSymb);
                    }
                }
                sb.delete(0, sb.length());
            }
            geneSet.clear();
        }

        availableGeneSet.clear();
        for (Gene gene : geneList) {
            double[] geneScore = geneMutNum.get(gene.geneSymb);
            if (geneScore == null) {
                gene.addFeatureValue(".");
                gene.addFeatureValue(".");
                gene.addFeatureValue(".");
                gene.addFeatureValue(".");
            } else {
                gene.addFeatureValue(String.valueOf(geneScore[0]));
                gene.addFeatureValue(String.valueOf(geneScore[1]));
                gene.addFeatureValue(String.valueOf(geneScore[2]));
                gene.addFeatureValue(String.valueOf(geneScore[3]));
            }
            sb = geneReadsTestInfo.get(gene.geneSymb);
            if (sb == null) {
                gene.addFeatureValue("");
            } else {
                gene.addFeatureValue(sb.toString());
            }
        }
    }

    public void summarizeAltNSVarPerGene(Chromosome chromosome, Set<Byte> dependentGeneFeature, Set<Byte> independentGeneFeature) throws Exception {
        List<Gene> geneList = chromosome.geneList;

        Map<String, double[]> geneMutNum = new HashMap<String, double[]>();
        Set<String> allAvailableGeneSet = new HashSet<String>();
        int het = 0, hom = 0;
        int total = 0, maxNum = 0;
        String geneFeatureAnnot;
        StringBuilder sb = new StringBuilder();
        Set<String> geneSet = new HashSet<String>();
        int index;
        String geneSymb = null, feature = null;
        Byte featureID;
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb == null) {
                continue;
            }

            geneFeatureAnnot = var.getRefGeneAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            geneFeatureAnnot = var.getgEncodeAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            geneFeatureAnnot = var.getKnownGeneAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            geneFeatureAnnot = var.getEnsemblGeneAnnot();
            if (geneFeatureAnnot != null) {
                sb.append(geneFeatureAnnot).append(';');
            }
            total = var.getAffectedRefHomGtyNum() + var.getAffectedHetGtyNum() + var.getAffectedAltHomGtyNum()
                    + var.getUnaffectedRefHomGtyNum() + var.getUnaffectedHetGtyNum() + var.getUnaffectedAltHomGtyNum();

            maxNum = (int) (0.05 * total);
            het = var.getAffectedHetGtyNum();
            het += var.getUnaffectedHetGtyNum();
            //if it is two large
            if (var.getAffectedRefHomGtyNum() > var.getAffectedAltHomGtyNum()) {
                hom = var.getAffectedAltHomGtyNum();
            } else {
                hom = var.getAffectedRefHomGtyNum();
            }
            if (var.getUnaffectedRefHomGtyNum() > var.getUnaffectedAltHomGtyNum()) {
                hom += var.getUnaffectedAltHomGtyNum();
            } else {
                hom += var.getUnaffectedRefHomGtyNum();
            }

            if (sb.length() > 0) {
                String[] items = Util.tokenize(sb.substring(0, sb.length() - 1), ';');
                for (String item : items) {
                    index = item.indexOf(':');
                    if (index > 0) {
                        geneSymb = item.substring(0, index);
                    }
                    if (geneSet.contains(geneSymb)) {
                        continue;
                    }
                    geneSet.add(geneSymb);

                    index = item.lastIndexOf(':');
                    if (index > 0) {
                        feature = item.substring(index + 1);
                    }
                    featureID = geneFeatureMap.get(feature);
                    double[] mutNums = geneMutNum.get(geneSymb);
                    if (mutNums == null) {
                        mutNums = new double[2];
                        Arrays.fill(mutNums, 0);
                        geneMutNum.put(geneSymb, mutNums);
                    }
                    // non-synonemous preidcted as Y, non-synonemous preidcted
                    // as N, synonemous
                    if (dependentGeneFeature.contains(featureID)) {
                        // only do it for missense variants; many stop loss
                        // variants have no proection at loss of function
                        // variants
                        mutNums[0] += (het + 2 * hom);
                    } else if (independentGeneFeature.contains(featureID)) {
                        mutNums[1] += (het + 2 * hom);
                    }

                    if (!allAvailableGeneSet.contains(geneSymb)) {
                        geneList.add(new Gene(geneSymb));
                        allAvailableGeneSet.add(geneSymb);
                    }
                }
                sb.delete(0, sb.length());
            }
            geneSet.clear();
        }
        allAvailableGeneSet.clear();

        for (Gene gene : geneList) {
            double[] geneScore = geneMutNum.get(gene.geneSymb);
            if (geneScore == null) {
                gene.addFeatureValue(".");
                gene.addFeatureValue(".");
            } else {
                gene.addFeatureValue(String.valueOf(geneScore[0]));
                gene.addFeatureValue(String.valueOf(geneScore[1]));

            }
        }

    }

    class RiskPredictionLogisticTask extends Task implements Callable<String> {

        List<Variant> variantList;
        List<CombOrders> combOrderList;
        List<String> names;
        int geneFeatureNum;
        int bestModelNum;

        public RiskPredictionLogisticTask(List<CombOrders> combOrderList, int bestModelNum, List<String> names, List<Variant> variantList, int geneFeatureNum) {
            this.variantList = variantList;
            this.combOrderList = combOrderList;
            this.names = names;
            this.geneFeatureNum = geneFeatureNum;
            this.bestModelNum = bestModelNum;
        }

        @Override
        public String call() throws Exception {
            RegressionParams tmpRP = null;
            RegressionParams bestRP = null;
            List<Integer> paramIndexes = new ArrayList<Integer>();
            Set<Integer> dataIndexes = new HashSet<Integer>();
            double tmpP, bestP, prior;
            StringBuilder tmpStrB = new StringBuilder();
            boolean isDeleteriousness = false;
            int combListSize = combOrderList.size();
            double sum = 0;
            double[] priors = new double[]{0.05, 0.01, 0.0001};
            double[] mafBins = new double[]{0.01, 0.02, 0.03};

            List<Integer> bestPParamIndexes = new ArrayList<Integer>();
            for (Variant var : variantList) {
                if (var.smallestFeatureID >= 15) {
                    var.setFeatureValue(geneFeatureNum, null);
                    var.setFeatureValue(geneFeatureNum + 1, null);
                    var.setFeatureValue(geneFeatureNum + 2, null);
                    continue;
                }
                if (var.scores1 != null) {
                    dataIndexes.clear();
                    for (int j = 0; j < var.scores1.length; j++) {
                        if (!Float.isNaN(var.scores1[j])) {
                            dataIndexes.add(j);
                        }
                    }
                    if (dataIndexes.isEmpty()) {
                        var.setFeatureValue(geneFeatureNum, null);
                        var.setFeatureValue(geneFeatureNum + 1, null);
                        var.setFeatureValue(geneFeatureNum + 2, null);
                        continue;
                    }
                    tmpRP = null;
                    bestRP = null;
                    bestP = 0;
                    isDeleteriousness = false;
                    tmpStrB.delete(0, tmpStrB.length());
                    for (int t = combListSize - 1; t >= combListSize - bestModelNum; t--) {
                        CombOrders cmbOrder = combOrderList.get(t);
                        /*
                         * if
                         * (!testSet.containsAll(combOrderList.get(t).indexes)||
                         * combOrderList.get(t).indexes.size()!=2) { continue; }
                         */

                        paramIndexes.clear();
                        if (dataIndexes.containsAll(cmbOrder.indexes)) {
                            paramIndexes.addAll(cmbOrder.indexes);
                            tmpRP = cmbOrder.rp;
                            sum = tmpRP.coef[0];
                            for (int j = 1; j < tmpRP.coef.length; j++) {
                                sum += (tmpRP.coef[j] * var.scores1[paramIndexes.get(j - 1)]);
                            }
                            // calculate the conditional probablity with MAF
                            // if (Float.isNaN(var.altAF) || var.altAF <= 0.01)
                            {
                                prior = ((1 - priors[0]) / priors[0]) * tmpRP.sampleCase2CtrRatio;
                                tmpP = 1 + prior * Math.exp(-sum);
                                tmpP = 1 / tmpP;
                            }
                            /*
                             * else if (var.altAF <= 0.03) { prior = ((1 -
                             * priors[1]) / priors[1]) *
                             * tmpRP.sampleCase2CtrRatio; tmpP = 1 + prior *
                             * Math.exp(-sum); tmpP = 1 / tmpP; } else { tmpP =
                             * 0; }
                             */

                            if (bestP <= tmpP) {
                                bestP = tmpP;
                                bestRP = tmpRP;
                                bestPParamIndexes.clear();
                                bestPParamIndexes.addAll(combOrderList.get(t).indexes);
                            }

                            if (tmpP >= tmpRP.optimalCutoff) {
                                var.setFeatureValue(geneFeatureNum, String.valueOf(tmpP));
                                var.setFeatureValue(geneFeatureNum + 1, "Y");
                                for (Integer ind : paramIndexes) {
                                    tmpStrB.append(names.get(ind));
                                    tmpStrB.append(';');
                                }
                                tmpStrB.deleteCharAt(tmpStrB.length() - 1);
                                tmpStrB.append(':');
                                tmpStrB.append(Util.doubleToString(tmpRP.optimalCutoff, 4));
                                tmpStrB.append(':');
                                tmpStrB.append(Util.doubleToString(tmpRP.truePositiveRate, 3));
                                tmpStrB.append(':');
                                tmpStrB.append(Util.doubleToString(tmpRP.trueNegativeRate, 3));
                                var.setFeatureValue(geneFeatureNum + 2, tmpStrB.toString());

                                isDeleteriousness = true;
                                break;
                            }
                        }
                    }

                    if (!isDeleteriousness) {
                        if (tmpRP == null) {
                            // keep varaint have no risk scores which may be
                            // safer
                            var.setFeatureValue(geneFeatureNum, null);
                            var.setFeatureValue(geneFeatureNum + 1, null);
                            var.setFeatureValue(geneFeatureNum + 2, null);
                        } else {
                            var.setFeatureValue(geneFeatureNum, String.valueOf(bestP));
                            var.setFeatureValue(geneFeatureNum + 1, "N");
                            if (!bestPParamIndexes.isEmpty()) {
                                for (Integer ind : bestPParamIndexes) {
                                    tmpStrB.append(names.get(ind));
                                    tmpStrB.append(';');
                                }
                                tmpStrB.deleteCharAt(tmpStrB.length() - 1);
                            }
                            tmpStrB.append(':');
                            tmpStrB.append(Util.doubleToString(bestRP.optimalCutoff, 4));
                            tmpStrB.append(':');
                            tmpStrB.append(Util.doubleToString(bestRP.truePositiveRate, 3));
                            tmpStrB.append(':');
                            tmpStrB.append(Util.doubleToString(bestRP.trueNegativeRate, 3));
                            var.setFeatureValue(geneFeatureNum + 2, tmpStrB.toString());

                        }
                    }
                } else {
                    // keep varaint have no risk scores which may be safer 
                    var.setFeatureValue(geneFeatureNum, null);
                    var.setFeatureValue(geneFeatureNum + 1, null);
                    var.setFeatureValue(geneFeatureNum + 2, null);
                }
            }
            return "";
        }
    }

    public void riskPredictionRareDiseaseAll(Chromosome chromosome, List<CombOrders> combOrderList, int bestModelNum, boolean filterNonDisMut, List<String> names, int maxThreadNum, FiltrationSummarySet dbNSFPMendelPred) throws Exception {
        // humvar predict
        double[] priors = new double[]{0.05, 0.01, 0.0001};
        double[] mafBins = new double[]{0.01, 0.02, 0.03};

        Set<String> geneSymbSet = new HashSet<String>();

        int counterDis = 0;
        int counterCon = 0;

        double sum = 0;
        RegressionParams tmpRP = null;
        RegressionParams bestRP = null;

        double tmpP, bestP;
        StringBuilder tmpStrB = new StringBuilder();
        boolean isDeleteriousness = false;
        List<Integer> bestPParamIndexes = new ArrayList<Integer>();
        int geneFeatureNum = dbNSFPMendelPred.getAvailableFeatureIndex();

        if (chromosome == null) {
            return;
        }
        ExecutorService exec = Executors.newFixedThreadPool(maxThreadNum);
        final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
        int varNum = chromosome.variantList.size();

        int[] blocks = org.cobi.util.thread.Util.partitionEvenBlock(maxThreadNum, 0, varNum);
        int blockNum = blocks.length - 1;
        int runningThread = 0;
        for (int s = 0; s < blockNum; s++) {
            RiskPredictionLogisticTask task = new RiskPredictionLogisticTask(combOrderList, bestModelNum, names, chromosome.variantList.subList(blocks[s], blocks[s + 1]), geneFeatureNum);
            serv.submit(task);
            runningThread++;
        }

        for (int s = 0; s < runningThread; s++) {
            Future<String> task = serv.take();
            String infor = task.get();
            //  System.out.println(infor);
        }

        exec.shutdown();
        List<Variant> keptVarList = new ArrayList<Variant>();
        int totalVarNum = 0;
        for (Variant var : chromosome.variantList) {
            if (var.featureValues[geneFeatureNum + 1] == null) {
                if (filterNonDisMut) {
                    keptVarList.add(var);
                }
            } else if (var.featureValues[geneFeatureNum + 1].equals("Y")) {
                counterDis++;
                if (filterNonDisMut) {
                    keptVarList.add(var);
                }
                if (var.geneSymb != null) {
                    geneSymbSet.add(var.geneSymb);
                }

            } else {
                // noly filter it for missense variants
                if (var.smallestFeatureID != 6) {
                    keptVarList.add(var);
                }
                counterCon++;
            }
        }

        if (filterNonDisMut) {
            chromosome.variantList.clear();
            chromosome.variantList.addAll(keptVarList);
            totalVarNum = keptVarList.size();
            keptVarList.clear();
            chromosome.buildVariantIndexMap();
        }

        dbNSFPMendelPred.increaseCount(0, counterDis);
        dbNSFPMendelPred.increaseCount(1, geneSymbSet.size());
        dbNSFPMendelPred.increaseCount(2, counterCon);
        dbNSFPMendelPred.increaseCount(3, totalVarNum);

    }

    public void dbscSNV(Chromosome chromosome, AnnotationSummarySet dbScSNV, boolean needProgressionIndicator) {
        indexChrom = 0;
        indexPosition = 1;
        indexREF = 2;
        indexALT = 3;
        int indexAda = 4;
        List<Variant> varList = chromosome.variantList;
        int chrID = chromosome.getId();
        BufferedReader br = dbScSNV.getBr();
        StringBuilder newLine = dbScSNV.getLastLine();

        if (varList.isEmpty()) {
            return;
        }
        int feautreNum = dbScSNV.getAvailableFeatureIndex();

        //No matched SNPs means null 
        String missingVal = ".";
        for (Variant var : varList) {
            var.setFeatureValue(feautreNum, missingVal);
        }
        String currentLine = null;
        try {
            String varChrom = chromosome.getName();
            int maxColNum = indexChrom;
            maxColNum = Math.max(maxColNum, indexPosition);
            maxColNum = Math.max(maxColNum, indexREF);
            maxColNum = Math.max(maxColNum, indexALT);
            maxColNum = Math.max(maxColNum, indexAda);

            int lineCounter = 0;

            int filePosition = -1;
            StringBuilder tmpBuffer = new StringBuilder();
            String ref;
            String alt;
            String mafStr = null;
            float score;
            String[] alts;
            String[] mafStrs;
            int[] varIndex = null;

            char[] backSpaces = null;
            int delNum = 0;
            int existVarNum = 0;
            StringBuilder sb = new StringBuilder();
            boolean hitOnce = false;

            String[] cells = null;
            if (newLine.length() == 0) {
                currentLine = br.readLine();
                if (currentLine == null) {
                    return;
                }
            } else {
                currentLine = newLine.toString();
                newLine.delete(0, newLine.length());
            }
            int fileChrID;
            do {
                /*
                 if (currentLine.indexOf("BP") >= 0 || currentLine.indexOf("bp") >= 0) {
                 continue;
                 }*/
                lineCounter++;
                if (needProgressionIndicator && lineCounter % 50000 == 0) {
                    String prog = String.valueOf(lineCounter);
                    System.out.print(prog);
                    backSpaces = new char[prog.length()];
                    Arrays.fill(backSpaces, '\b');
                    System.out.print(backSpaces);
                }

                //StringTokenizer st = new StringTokenizer(currentLine.trim());
                cells = Util.tokenize(currentLine, '\t');
                if (cells.length < 2) {
                    cells = Util.tokenizeIngoreConsec(currentLine, ' ');
                }
                //initialize varaibles

                mafStrs = null;

                fileChrID = chromNameIndexMap.get(cells[indexChrom]);
                if (chrID < fileChrID) {
                    newLine.append(currentLine);
                    break;
                } else if (chrID > fileChrID) {
                    continue;
                }
                filePosition = Util.parseInt(cells[indexPosition]);
                ref = cells[indexREF];
                alt = cells[indexALT];
                if (cells.length > indexAda) {
                    mafStr = cells[indexAda];
                }

                alts = alt.split("N");
                if (mafStr != null) {
                    mafStrs = mafStr.split(",", -1);
                }

                //  System.err.println(currentLine);
                hitOnce = false;
                int tmpPos = 0;
                //once the variant is in db, it at least has a zero freq
                for (int s = 0; s < alts.length; s++) {
                    if (alts[s] == null || alts[s].isEmpty()) {
                        continue;
                    }
                    score = Float.NaN;
                    if (mafStrs != null && s < mafStrs.length) {
                        //this missing score is denoted by .
                        if (mafStrs[s] != null && !mafStrs[s].isEmpty() && !mafStrs[s].equals(".")) {
                            score = Util.parseFloat(mafStrs[s]);
                        }
                    }
                    tmpPos = filePosition;

                    alt = alts[s];

                    varIndex = chromosome.lookupVariantIndexes(tmpPos);
                    if (varIndex == null) {
                        continue;
                    }

                    // System.out.println(fileChr);
                    for (int index : varIndex) {
                        Variant var = varList.get(index);
                        if (var.isIndel) {
                            continue;
                        } else {
                            String[] altAlleles = var.getAltAlleles();
                            for (String str : altAlleles) {
                                if (str.charAt(0) == alt.charAt(0)) {
                                    hitOnce = true;
                                    if (Float.isNaN(score)) {
                                        var.setFeatureValue(feautreNum, ".");
                                    } else {
                                        var.setFeatureValue(feautreNum, String.valueOf(score));
                                    }
                                    break;
                                }
                            }
                        }

                    }
                }
                if (hitOnce) {
                    existVarNum++;
                }
            } while ((currentLine = br.readLine()) != null);

            dbScSNV.setLeftNum(existVarNum + dbScSNV.getLeftNum());
            dbScSNV.setAnnotNum(existVarNum + dbScSNV.getAnnotNum());
            dbScSNV.setTotalNum(lineCounter + dbScSNV.getTotalNum());
            if (needProgressionIndicator) {
                backSpaces = new char[7];
                Arrays.fill(backSpaces, '\b');
                System.out.print(backSpaces);
            }
        } catch (Exception ex) {
            if (currentLine != null) {
                System.err.println("Errors in a row: " + currentLine);
            }
            ex.printStackTrace();
        }
    }

    public void geneFeatureAnnot(Chromosome chromosome, int chrID, ReferenceGenome refGenome, Set<Byte> featureInSet, int featureNum, int maxThreadNum) throws Exception {
        ExecutorService exec = Executors.newFixedThreadPool(maxThreadNum);
        final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
        int runningThread = 0;
        int varNum = chromosome.variantList.size();

        int[] blocks = org.cobi.util.thread.Util.partitionEvenBlock(maxThreadNum, 0, varNum);
        int blockNum = blocks.length - 1;
        for (int s = 0; s < blockNum; s++) {
            GeneFeatureAnnotTask task = new GeneFeatureAnnotTask(chromosome.variantList.subList(blocks[s], blocks[s + 1]), chrID, refGenome, featureInSet, featureNum);
            serv.submit(task);
            runningThread++;
        }

        for (int s = 0; s < runningThread; s++) {
            Future<String> task = serv.take();
            String infor = task.get();
            //  System.out.println(infor);
        }

        exec.shutdown();
    }

    public void geneFeatureAnnot(Chromosome chromosome, int chrID, ReferenceGenome refGenome, Set<Byte> featureInSet, int feautreNum) throws Exception {

        if (refGenome.getName().equals("refgene")) {
            if (chromosome == null) {
                return;
            }
            for (Variant var : chromosome.variantList) {
                GeneFeature gf = refGenome.getVarFeature(STAND_CHROM_NAMES[chrID], var, true, featureInSet, new RNABoundaryIndex(0));

                // Feature for UniProtFeature
                if (gf.getInfor() == null) {
                    var.setFeatureValue(feautreNum, ".");
                } else {
                    var.setFeatureValue(feautreNum, gf.getInfor());
                }
                var.setRefGeneAnnot(gf.getName());
            }

        } else if (refGenome.getName().equals("gencode")) {
            if (chromosome == null) {
                return;
            }
            for (Variant var : chromosome.variantList) {
                GeneFeature gf = refGenome.getVarFeature(STAND_CHROM_NAMES[chrID], var, true, featureInSet, new RNABoundaryIndex(0));

                // Feature for UniProtFeature
                if (gf.getInfor() == null) {
                    var.setFeatureValue(feautreNum, ".");
                } else {
                    var.setFeatureValue(feautreNum, gf.getInfor());
                }
                var.setgEncodeAnnot(gf.getName());
            }

        } else if (refGenome.getName().equals("knowngene")) {
            if (chromosome == null) {
                return;
            }
            for (Variant var : chromosome.variantList) {
                GeneFeature gf = refGenome.getVarFeature(STAND_CHROM_NAMES[chrID], var, true, featureInSet, new RNABoundaryIndex(0));

                // Feature for UniProtFeature
                if (gf.getInfor() == null) {
                    var.setFeatureValue(feautreNum, ".");
                } else {
                    var.setFeatureValue(feautreNum, gf.getInfor());
                }
                var.setKnownGeneAnnot(gf.getName());
            }

        } else if (refGenome.getName().equals("ensembl")) {
            if (chromosome == null) {
                return;
            }
            for (Variant var : chromosome.variantList) {
                GeneFeature gf = refGenome.getVarFeature(STAND_CHROM_NAMES[chrID], var, true, featureInSet, new RNABoundaryIndex(0));

                // Feature for UniProtFeature
                if (gf.getInfor() == null) {
                    var.setFeatureValue(feautreNum, ".");
                } else {
                    var.setFeatureValue(feautreNum, gf.getInfor());
                }
                var.setEnsemblGeneAnnot(gf.getName());

            }

        }
    }

    class PathwayP {

        double p;
        String pathway;

        public PathwayP(double p, String pathway) {
            this.p = p;
            this.pathway = pathway;
        }
    }

    class PathwayPComparator implements Comparator<PathwayP> {

        @Override
        public int compare(PathwayP arg0, PathwayP arg1) {
            return Double.compare(arg0.p, arg1.p);
        }
    }

    private String[] searchIbdRegions(List<String[]> regionItems, String chrom, int poss) {
        // to do
        String[] region = null;
        for (String[] cells : regionItems) {
            if (cells[0].equals(chrom)) {
                int start = Util.parseInt(cells[1]);
                int end = Util.parseInt(cells[2]);
                if (poss >= start && poss <= end) {
                    region = new String[2];
                    region[0] = "chr" + chrom + ":" + cells[1] + "-" + cells[2];
                    region[1] = String.valueOf(Util.parseInt(cells[2]) - Util.parseInt(cells[1]));
                    break;
                }

            }
        }
        return region;
    }

    public void pubMedIDIdeogramExploreGene(Chromosome chromosome, AnnotationSummarySet ass, List<String> pubmedMeshList, List<String[]> ideogramItems) throws Exception {
        if (!pubmedMeshList.isEmpty() && pubmedMeshList.get(0).equals("ANY")) {
            return;
        }

        NCBIRetriever ncbiRetriever = new NCBIRetriever();
        int account = 1;
        List<String> ideoGrams = null;
        Map<String, String> historyResultMap = new HashMap<String, String>();
        StringBuilder result = new StringBuilder();

        if (chromosome == null) {
            return;
        }

        for (mRNA mrna : chromosome.mRNAList) {
            result.delete(0, result.length());
            ideoGrams = searchIdeogramRegions(ideogramItems, chromosome.getName(), mrna.getStart(), mrna.getEnd());
            for (String reg : ideoGrams) {
                // ignore regions which are two broad
                /*
                 * if (reg.indexOf(".") < 0) { continue; }
                 */
                String pubIDs = historyResultMap.get(reg);
                if (pubIDs == null) {
                    LOG.info(account + ": Searching NCBI PubMed for " + pubmedMeshList.toString() + " and " + reg);
                    while ((pubIDs = ncbiRetriever.pubMedIDESearch(pubmedMeshList, reg)) == null) {
                        // System.out.print("reconnecting...");
                    }
                    // System.out.println(pubIDs);
                    if (pubIDs == null) {
                        pubIDs = "";
                    }
                    historyResultMap.put(reg, pubIDs);
                    account++;
                }
                if (pubIDs.trim().length() > 0) {
                    result.append(reg).append(":[").append(pubIDs).append("] ");
                }
            }
            mrna.addFeatureValue(result.toString());
        }

        ass.setAnnotNum(ass.getAnnotNum() + account);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - account);
    }

    private List<String> searchIdeogramRegions(List<String[]> ideogramItems, String chrom, int startPoss, int endPoss) {
        // to do
        List<String> ideoRegions = new ArrayList<String>();
        for (String[] cells : ideogramItems) {
            if (cells[0].equals(chrom)) {
                if (cells[1] != null && cells[1].trim().length() > 0) {
                    int start = Util.parseInt(cells[2]);
                    int end = Util.parseInt(cells[3]);
                    if (startPoss >= start && startPoss <= end) {
                        ideoRegions.add(cells[0] + cells[1]);
                    } else if (endPoss >= start && endPoss <= end) {
                        ideoRegions.add(cells[0] + cells[1]);
                    }
                }
            }
        }
        // normall the detailed region in the latter
        Collections.reverse(ideoRegions);
        return ideoRegions;

    }

    private List<String> searchIdeogramRegions(List<String[]> ideogramItems, String chrom, int poss) {
        // to do
        List<String> ideoRegions = new ArrayList<String>();
        for (String[] cells : ideogramItems) {
            if (cells[0].equals(chrom)) {
                if (cells[1] != null && cells[1].trim().length() > 0) {
                    int start = Util.parseInt(cells[2]);
                    int end = Util.parseInt(cells[3]);
                    if (poss >= start && poss <= end) {
                        ideoRegions.add(cells[0] + cells[1]);
                    }
                }
            }
        }
        // normall the detailed region in the latter
        Collections.reverse(ideoRegions);
        return ideoRegions;

    }

}
