/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.jet.stat.Probability;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import jsc.independentsamples.MannWhitneyTest;
import jsc.tests.H1;
import org.apache.log4j.Logger;
import org.cobi.kggseq.GlobalManager;
import org.cobi.kggseq.entity.ArrayListComparatorDouble;
import org.cobi.kggseq.entity.Chromosome;
import org.cobi.kggseq.entity.Gene;
import org.cobi.kggseq.entity.Genome;
import org.cobi.kggseq.entity.GeneSet;
import org.cobi.kggseq.entity.GeneSetWilcoxPValueComparator;
import org.cobi.kggseq.entity.RefGene;
import org.cobi.kggseq.entity.SeqSegment;
import org.cobi.kggseq.entity.Variant;
import org.cobi.util.file.LocalFileFunc;
import org.cobi.util.plot.PValuePainter;
import org.cobi.util.stat.MultipleTestingMethod;
import org.cobi.util.text.BGZFInputStream;
import org.cobi.util.text.LocalExcelFile;
import org.cobi.util.text.StringArrayDoubleComparator;
import org.cobi.util.text.Util;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;

/**
 *
 * @author mxli
 */
public class GeneAnnotator {

    private static final Logger LOG = Logger.getLogger(GeneAnnotator.class);

    public void readZebrafishPhenotype(Map<String, String> xebraP2Phe) throws Exception {
        File geneSymbols1 = new File(GlobalManager.RESOURCE_PATH + "/HgncGene.txt");
        BufferedReader br = new BufferedReader(new FileReader(geneSymbols1));
        String strLine;
        Map<String, String> geneEtrIDSymbMap = new HashMap<String, String>();
//                    int i=1;
        while ((strLine = br.readLine()) != null) {
//                        System.out.println(i++);
            String[] strItems = strLine.split("\t", -1);
            if (strItems[9] != null) {
                strItems[9] = strItems[9].trim();
                if (strItems[9].length() > 0) {
                    geneEtrIDSymbMap.put(strItems[9], strItems[1].toUpperCase());
                }
            }
        }
        br.close();

        File fisher = new File(GlobalManager.RESOURCE_PATH + "/pheno_fish.txt.gz");
        br = LocalFileFunc.getBufferedReader(fisher.getCanonicalPath());
        //skip the headline
        br.readLine();
        br.readLine();
        Map<String, Set> xebraP2Phe0 = new HashMap<String, Set>();

        while ((strLine = br.readLine()) != null) {
            String[] strItems = strLine.split("\t", -1);
            String geneSymb = geneEtrIDSymbMap.get(strItems[2]);
            if (geneSymb == null) {
                continue;
            }
            Set pheno1 = xebraP2Phe0.get(geneSymb);
            if (pheno1 == null) {
                pheno1 = new HashSet<String>();
                xebraP2Phe0.put(geneSymb, pheno1);
            }

            if (strItems[9] != null && strItems[9].length() > 0) {
                pheno1.add(strItems[9]);
            }
            if (strItems[18] != null && strItems[18].length() > 0) {
                pheno1.add(strItems[18]);
            }
        }
        for (Map.Entry<String, Set> mPath : xebraP2Phe0.entrySet()) {
            xebraP2Phe.put(mPath.getKey(), mPath.getValue().toString());
        }
        xebraP2Phe0.clear();
        geneEtrIDSymbMap = null;
        br.close();
    }

    public void retrieveVariants(RefGene mgeneExons, String scorePath, BGZFInputStream bfScoreIndexes, String mafPath, BGZFInputStream bfMAFIndexes, int actualThreadNum, int[] scoreIndexes) throws Exception {
        IntArrayList exonBounders = new IntArrayList();
        for (SeqSegment mgs : mgeneExons.getMergedSegments()) {
            exonBounders.add(mgs.getStart());
            exonBounders.add(mgs.getEnd());
        }
        int boundNum = exonBounders.size();
        if (exonBounders.getQuick(0) > exonBounders.getQuick(boundNum - 1)) {
            exonBounders.reverse();
        }

        long[] pos = bfScoreIndexes.findIndex(exonBounders.getQuick(0), exonBounders.getQuick(boundNum - 1));
        if (pos[0] == pos[1]) {
            return;
        }

        File rsFile = new File(scorePath);
        final BGZFInputStream bfScore = new BGZFInputStream(rsFile.getCanonicalPath(), 1, true);
        bfScore.adjustPos(pos[0], pos[1]);
        bfScore.creatSpider(pos[0] != 0);
        int indexCHROM = -1, indexPOS = 0;
        VarAnnotTask varAnnoTask = new VarAnnotTask(bfScore.spider[0], indexCHROM, indexPOS, scoreIndexes);
        IntArrayList positions = new IntArrayList();
        List<char[]> alleleList = new ArrayList<char[]>();
        List<float[]> scoreList = new ArrayList<float[]>();
        int[] regions = new int[exonBounders.size()];
        for (int i = 0; i < regions.length; i++) {
            regions[i] = exonBounders.getQuick(i);
        }
        List<double[]> frequnceList = new ArrayList<double[]>(positions.size());
        for (int i = positions.size() - 1; i >= 0; i++) {
            frequnceList.add(new double[]{Double.NaN});
        }
        varAnnoTask.dbNSFPAnnotSimple(positions, alleleList, scoreList, regions);
        varAnnoTask.markByANNOVARefFormatSimpleVar(positions, alleleList, frequnceList);
//should define the format of alleles and frequencies
        File mapfFile = new File(mafPath);
        final BGZFInputStream bfMAF = new BGZFInputStream(mapfFile.getCanonicalPath(), 1, true);
        bfMAF.adjustPos(pos[0], pos[1]);
        bfMAF.creatSpider(pos[0] != 0);

        int sss = 0;

    }

    public void geneScoreSimulation(Map<String, double[]> geneScores, String chromName, Map<String, RefGene> mergedGeneExonMap, String dbNSFPFilePath, String alleFreqFilePath, int actualThreadNum, int[] scoreIndexes) {
        try {
            File rsFile = new File(dbNSFPFilePath);
            BGZFInputStream bfScoreIndexes = new BGZFInputStream(rsFile.getCanonicalPath(), actualThreadNum, true);
            if (!bfScoreIndexes.checkIndex()) {
                bfScoreIndexes.adjustPos();
                bfScoreIndexes.buildIndex(rsFile.getCanonicalPath(), -1, 0, true);
            }
            bfScoreIndexes.readIndex(false, null);

            File mafFile = new File(alleFreqFilePath);
            BGZFInputStream bfMAFIndexes = new BGZFInputStream(mafFile.getCanonicalPath(), actualThreadNum, true);
            if (!bfMAFIndexes.checkIndex()) {
                bfMAFIndexes.adjustPos();
                bfMAFIndexes.buildIndex(rsFile.getCanonicalPath(), 0, 1, false);
            }
            bfMAFIndexes.readIndex(false, chromName);

            for (Map.Entry<String, double[]> item : geneScores.entrySet()) {
                String geneSymbDis = item.getKey();
                RefGene mgeneExons = mergedGeneExonMap.get(geneSymbDis);
                
                retrieveVariants(mgeneExons, dbNSFPFilePath, bfScoreIndexes, alleFreqFilePath, bfMAFIndexes, actualThreadNum, scoreIndexes);
//                sss
            }
        } catch (Exception ex) {

        }

    }

    public void readMousePhenotype(Map<String, String[]> g2MP) throws Exception {
        File fleMouse1 = new File(GlobalManager.RESOURCE_PATH + "/" + "HMD_HumanPhenotype.rpt.txt");
        File fleMouse2 = new File(GlobalManager.RESOURCE_PATH + "/" + "VOC_MammalianPhenotype.rpt.txt");
        File fleMouse3 = new File(GlobalManager.RESOURCE_PATH + "/" + "ALL_genotype_phenotype.mouse.gz");

        BufferedReader br = new BufferedReader(new FileReader(fleMouse2));
        String strLine;
        Map<String, String> mp2Phe1 = new HashMap<String, String>();
//                    int i=1;
        while ((strLine = br.readLine()) != null) {
//                        System.out.println(i++);
            String[] strItems = strLine.split("\t", -1);
            if (strItems.length > 2) {
                if (strItems[2] != null && strItems[2].length() > 0) {
                    mp2Phe1.put(strItems[0], strItems[1] + "->" + strItems[2]);
                } else {
                    mp2Phe1.put(strItems[0], strItems[1]);
                }
            } else {
                mp2Phe1.put(strItems[0], strItems[1]);
            }
        }
        br.close();

        String strPhe1 = null;

        br = new BufferedReader(new FileReader(fleMouse1));
        while ((strLine = br.readLine()) != null) {
            String[] strItems = strLine.split("\t", -1);
            if (strItems.length < 7 || strItems[6] == null || strItems[6].trim().length() == 0) {
                continue;
            }

            strPhe1 = ".";
            String[] strMP = strItems[6].trim().split(" ");
            for (String temp : strMP) {
                if (mp2Phe1.containsKey(temp)) {
                    if (strPhe1.length() == 1) {
                        strPhe1 = mp2Phe1.get(temp);
                    } else {
                        strPhe1 += ("||" + mp2Phe1.get(temp));
                    }
                }
            }
            g2MP.put(strItems[0], new String[]{strPhe1, "."});
        }
        br.close();

        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fleMouse3))));
        String geneSymb = null;
        while ((strLine = br.readLine()) != null) {
            // String[] strItems = strLine.split(",(?=([^\"]*\"[^\"]*\")*[^\"]*$)", -1);
            String[] strItems = strLine.split("\t", -1);
            geneSymb = strItems[1].toUpperCase().trim();
            String[] values = g2MP.get(geneSymb);

            if (values == null) {
                g2MP.put(geneSymb, new String[]{".", strItems[20].replace("\"", "") + "," + strItems[22].replace("\"", "")});
            } else {
                values[1] = strItems[20].replace("\"", "") + "," + strItems[22].replace("\"", "");
            }
        }
        br.close();
    }

    public Map<String, double[]> readMutationGeneScore(String geneCoVarFilePath, List<String> heads) throws Exception {
        Map<String, double[]> driverGeneScores = new HashMap<String, double[]>();
        String line;
        String[] cells;
        BufferedReader brScore = new BufferedReader(new FileReader(geneCoVarFilePath));
        // "expr",  , "GIH", "JPT", "LWK", "MEX", "MKK", "YRI"
        String[] names = {"gene", "expr", "reptime", "hic", "constraint_score"};
        //String[] names = {"gene", "reptime", "hic", "constraint_score", "YRI"};
        //   names = null;
        if (names == null || names.length == 0) {
            return driverGeneScores;
        }
        int[] effecIndex = new int[names.length];
        Arrays.fill(effecIndex, -1);

        line = brScore.readLine();
        cells = line.split("\t");

        for (int i = 0; i < effecIndex.length; i++) {
            for (int j = 0; j < cells.length; j++) {
                if (names[i].equals(cells[j])) {
                    effecIndex[i] = j;
                    break;
                }
            }
        }

        for (int i = 1; i < effecIndex.length; i++) {
            heads.add(cells[effecIndex[i]]);
        }

        while ((line = brScore.readLine()) != null) {
            if (line.trim().length() == 0) {
                continue;
            }
            cells = line.split("\t");
            double[] itemCell = new double[effecIndex.length - 1];
            Arrays.fill(itemCell, Double.NaN);
            boolean hasError = false;
            for (int i = 1; i < effecIndex.length; i++) {
                hasError = !Util.isNumeric(cells[effecIndex[i]]) || cells[effecIndex[i]].equals("NaN");
                if (hasError) {
                    break;
                }
                itemCell[i - 1] = Double.parseDouble(cells[effecIndex[i]]);
            }
            if (hasError) {
                continue;
            }
            //convert the expression 
            if (itemCell[0] > 0) {
                itemCell[0] = Math.log10(itemCell[0]);
            }
            driverGeneScores.put(cells[effecIndex[0]], itemCell);
        }
        brScore.close();
        return driverGeneScores;
    }

    public Map<String, Map<String, Integer>> readCosmicGeneAnnotation(String dbPath) throws Exception {
        int indexChrom = 0;
        int indexref = 3;
        int indexalt = 4;
        int indexaaref = 3;
        int indexaaalt = 4;
        int indexhg19pos = 1;
        int indexCancerInfo = 7;
        int indexGene = 4;

        int maxColNum = indexChrom;

        maxColNum = Math.max(maxColNum, indexCancerInfo);
        maxColNum = Math.max(maxColNum, indexGene);

        String currentLine = null;

        boolean incomplete;

        String geneStr = null;
        String cancerInfo = null;

        StringBuilder tmpBuffer = new StringBuilder();
        int cosmicVarNum = 0;

        File rsFile = new File(dbPath);
        if (!rsFile.exists()) {
            LOG.error(rsFile.getCanonicalPath() + " does not exist!");
            return null;
        }
        // System.out.print(" Chromosome " + Options.REF_CHROM_NAMES[chromID]);

        BufferedReader br = LocalFileFunc.getBufferedReader(rsFile.getCanonicalPath());
        // skip to the head line
        br.readLine();

        Map<String, Map<String, Integer>> cosmicGeneMut = new HashMap<String, Map<String, Integer>>();

        int num = 0;
        Integer num1;
        while ((currentLine = br.readLine()) != null) {
            // System.out.println(currentLine);
            String[] cells = Util.tokenize(currentLine, '\t', maxColNum);
            // initialize varaibles
            incomplete = true;

            geneStr = null;
            cancerInfo = null;
            geneStr = cells[indexGene];
            cancerInfo = cells[indexCancerInfo];;

            Map<String, Integer> canNum = cosmicGeneMut.get(geneStr);
            if (canNum == null) {
                canNum = new HashMap<String, Integer>();
                cosmicGeneMut.put(geneStr, canNum);
            }
            cancerInfo = cancerInfo.substring(1, cancerInfo.length() - 1);
            String[] cells1 = cancerInfo.split(",");
            for (String cell : cells1) {
                cells = Util.tokenize(cell, '=');
                num = Integer.parseInt(cells[1].trim());
                cells[0] = cells[0].trim();
                num1 = canNum.get(cells[0]);
                if (num1 != null) {
                    num += num1;
                }
                canNum.put(cells[0], num);
            }
            cosmicVarNum++;
        }
        br.close();
        StringBuilder info = new StringBuilder();
        info.append(cosmicVarNum).append(" variants in the COSMIC database are read.");
        LOG.info(info);
        return cosmicGeneMut;
    }

    public void enrichmentTestGeneSet(Map<String, GeneSet> dbPathwaySet, Map<String, double[]> genePValueMap, int minSetSize, double genePValueCutoff,
            String geneSumOutFile) throws Exception {

        DoubleArrayList zValueListOutPathway = new DoubleArrayList();
        DoubleArrayList zValueListInGeneSet = new DoubleArrayList();

        List<String[]> availableGWASGeneList = new ArrayList<String[]>();
        double[] p1 = null;
        int subPopulationSize;
        DoubleArrayList geneSetPValuesForQQPlot = new DoubleArrayList();

        List<GeneSet> effectiveGeneSetList = new ArrayList<GeneSet>();
        double p;
        for (Map.Entry<String, GeneSet> mPath : dbPathwaySet.entrySet()) {
            GeneSet curPath = mPath.getValue();
            String pathID = mPath.getKey();

            long start1 = System.currentTimeMillis();
            Set<String> pGenes = curPath.getGeneSymbols();
            zValueListInGeneSet.clear();
            zValueListOutPathway.clear();
            availableGWASGeneList.clear();

            Set<String> geneSymbolWeightMap = curPath.getGeneSymbols();

            for (String pGeneSymb : geneSymbolWeightMap) {
                //may introduce biase
                if (!genePValueMap.containsKey(pGeneSymb)) {
                    continue;
                }

                p1 = genePValueMap.get(pGeneSymb);
                if (p1 == null) {
                    continue;
                }
                if (Double.isNaN(p1[0])) {
                    continue;
                }
                availableGWASGeneList.add(new String[]{pGeneSymb, String.valueOf(p1[0])});
                zValueListInGeneSet.add(p1[0]);
            }
            if (availableGWASGeneList.isEmpty()) {
                continue;
            }
            GeneSet newPathway = new GeneSet(pathID, pathID, curPath.getURL());
            newPathway.getGeneValues().addAll(availableGWASGeneList);
            for (Map.Entry<String, double[]> mGene : genePValueMap.entrySet()) {
                if (!pGenes.contains(mGene.getKey())) {
                    if (Double.isNaN(mGene.getValue()[0])) {
                        continue;
                    }
                    zValueListOutPathway.add(mGene.getValue()[0]);
                }
            }

            subPopulationSize = zValueListInGeneSet.size();
            if (subPopulationSize < minSetSize) {
                continue;
            }
            double[] zInGeneSet = new double[subPopulationSize];
            for (int t = 0; t < subPopulationSize; t++) {
                zInGeneSet[t] = zValueListInGeneSet.getQuick(t);
            }

            double[] zOutGeneSet = new double[zValueListOutPathway.size()];
            for (int t = 0; t < zOutGeneSet.length; t++) {
                zOutGeneSet[t] = zValueListOutPathway.getQuick(t);
            }
            /*
             if (pathID.equals("ADENYL_NUCLEOTIDE_BINDING")) {
             BufferedWriter testWriter = new BufferedWriter(new FileWriter("test.txt"));
             testWriter.write("group\tvalue\n");
             for (int t = 0; t < zInGeneSet.length; t++) {
             testWriter.write("1\t" + zInGeneSet[t] + "\n");
             }
             for (int t = 0; t < zOutGeneSet.length; t++) {
             testWriter.write("2\t" + zOutGeneSet[t] + "\n");
             }
             testWriter.close();
             }
             */

            MannWhitneyTest mt = new MannWhitneyTest(zInGeneSet, zOutGeneSet, H1.LESS_THAN);
            // WilcoxonTest wt = new WilcoxonTest(zInGeneSet, median, H1.LESS_THAN);
            // newPathway.setWilcoxonPValue(mt.getSP());
            // allPValues.add(newPathway.getWilcoxonPValue());
            //outPathwayList.add(newPathway);                
            //  runningResultTopComp.insertText(pathID + 2 + " Time " + (end - start));

            p = mt.getSP();
            newPathway.setWilcoxonPValue(p);

            geneSetPValuesForQQPlot.add(p);
            effectiveGeneSetList.add(newPathway);

        }

        ArrayListComparatorDouble acd = new ArrayListComparatorDouble(1);
        List<String[]> outGeneSetList = new ArrayList<String[]>();
        outGeneSetList.add(new String[]{"GeneSetID", "WilcoxonP", "Gene", "GeneP"});
        Collections.sort(effectiveGeneSetList, new GeneSetWilcoxPValueComparator());
        for (GeneSet pathway : effectiveGeneSetList) {
            List<String[]> geneList = pathway.getGeneValues();

            String[] items = new String[4];
            items[0] = pathway.getID();
            items[1] = String.valueOf(pathway.getWilcoxonPValue());

            Collections.sort(geneList, acd);
            items[2] = geneList.get(0)[0];
            items[3] = geneList.get(0)[1];
            outGeneSetList.add(items);
            for (int i = 1; i < geneList.size(); i++) {
                items = new String[4];
                items[2] = geneList.get(i)[0];
                items[3] = geneList.get(i)[1];
                outGeneSetList.add(items);
            }
        }

        LocalExcelFile.writeArray2XLSXFile(geneSumOutFile, outGeneSetList, true, 0, 2);
        String info = "The results of geneset enrichment test are saved in " + (new File(geneSumOutFile)).getCanonicalPath() + "!";

        LOG.info(info);
        geneSetPValuesForQQPlot.quickSort();
        List<DoubleArrayList> pvalueLists = new ArrayList<DoubleArrayList>();
        List<String> nameList = new ArrayList<String>();
        pvalueLists.add(geneSetPValuesForQQPlot);
        PValuePainter pvPainter = new PValuePainter(450, 450);
        File plotFile = new File(geneSumOutFile.substring(0, geneSumOutFile.length() - 5) + "qq.png");
        nameList.add("WilcoxonP");
        pvPainter.drawMultipleQQPlot(pvalueLists, nameList, null, plotFile.getCanonicalPath(), 1E-10);
        info = "The QQ plot saved in " + plotFile.getCanonicalPath();
        LOG.info(info);

        StringBuilder sumSb = new StringBuilder();
        double adjGenePValueCutoff = MultipleTestingMethod.benjaminiHochbergFDR(genePValueCutoff, geneSetPValuesForQQPlot);
        Set<String> wholeSigCandiGeneSet = new HashSet<String>();
        for (GeneSet pathway : effectiveGeneSetList) {
            if (pathway.getWilcoxonPValue() <= adjGenePValueCutoff) {
                wholeSigCandiGeneSet.add(pathway.getID());
            }
        }

        sumSb.append(wholeSigCandiGeneSet.size()).append(" genesets with p-value <= ").append(adjGenePValueCutoff).append(" pass the Benjamini Hochberg FDR q value cutoff ")
                .append(genePValueCutoff).append(".\n");
        LOG.info(sumSb.toString());
        wholeSigCandiGeneSet.clear();
    }

    public void geneMutationRateTest(Genome genome, Map<String, Double> geneLengths, Map<String, Map<String, Integer>> cosmicGeneMut, Map<String, double[]> geneScores,
            List<String> scoreHeads, double genePValueCutoff, Map<String, double[]> geneMutRegPValueAllVar, String geneSumOutFile) throws Exception {
        Chromosome[] chroms = genome.getChromosomes();
        int somaticNSVarIndex = -1;
        int somaticSVarIndex = -1;
        int somaticNSRatioIndex = -1;
        int somaticSRatioIndex = -1;

        boolean needPubMedAnno = false;
        boolean needCosmic = false;
        int shiftCol = 0;

        List<String> featureLabels = genome.getGeneFeatureLabels();
        for (int i = 0; i < featureLabels.size(); i++) {
            if (featureLabels.get(i).equals("#DependentVar")) {
                somaticNSVarIndex = i;
            } else if (featureLabels.get(i).equals("#IndependentVar")) {
                somaticSVarIndex = i;
            } else if (featureLabels.get(i).equals("NonsynonymousReadsRatio")) {
                somaticNSRatioIndex = i;
            } else if (featureLabels.get(i).equals("SynonymousReadsRatio")) {
                somaticSRatioIndex = i;
            }
            if (somaticSVarIndex >= 0 && somaticNSVarIndex >= 0 && somaticNSRatioIndex >= 0 && somaticSRatioIndex >= 0) {
                break;
            }
        }

        if (cosmicGeneMut != null) {
            needCosmic = true;
        }

        String[] basicHeadRow = new String[]{"GeneSymbol", "#DependentVar", "#IndependentVar", "AvgCodingLen"};
        List<String> geneTableHeadRow = new ArrayList<String>();

        geneTableHeadRow.addAll(Arrays.asList(basicHeadRow));

        if (needPubMedAnno) {
            geneTableHeadRow.add("PubMedID");
            shiftCol++;
        }
        if (needCosmic) {
            geneTableHeadRow.add("COSMICCancerInfo");
            shiftCol++;
        }

        for (String name : scoreHeads) {
            geneTableHeadRow.add(name);
        }

        geneTableHeadRow.add("RegResidue");
        geneTableHeadRow.add("RegP");
        boolean regReads = false;
        if (regReads) {
            geneTableHeadRow.add("RegReadsResidue");
            geneTableHeadRow.add("RegReadsP");
        }

        RConnection rcon = new RConnection();
        try {
            rcon.eval("pack=\"MASS\";  if (!require(pack,character.only = TRUE))      {        install.packages(pack,dep=TRUE,repos='http://cran.us.r-project.org');        if(!require(pack,character.only = TRUE)) stop(\"Package not found\")    }");
            rcon.eval("library(MASS)");
        } catch (RserveException ex) {
            LOG.error("KGGSeq failed to install R package MASS! Please report this problem to your administrator.");
        }
        try {
            rcon.eval("pack=\"countreg\";  if (!require(pack,character.only = TRUE))      {        install.packages(pack,dep=TRUE,repos='http://R-Forge.R-project.org');        if(!require(pack,character.only = TRUE)) stop(\"Package not found\")    }");
            rcon.eval("library(countreg)");
        } catch (RserveException ex) {
            LOG.error("KGGSeq failed to install R package countreg! Please report this problem to your administrator.");
        }
        if (regReads) {
            rcon.eval("pack=\"mvtnorm\";  if (!require(pack,character.only = TRUE))      {        install.packages(pack,dep=TRUE,repos='http://cran.us.r-project.org');        if(!require(pack,character.only = TRUE)) stop(\"Package not found\")    }");
            rcon.eval("library(mvtnorm)");
            rcon.eval("pack=\"tmvtnorm\";  if (!require(pack,character.only = TRUE))      {        install.packages(pack,dep=TRUE,repos='http://cran.us.r-project.org');        if(!require(pack,character.only = TRUE)) stop(\"Package not found\")    }");
            rcon.eval("library(tmvtnorm)");
        }

        List<String[]> geneMutRateSheet = new ArrayList<String[]>();
        List<double[]> countsRegList = new ArrayList<double[]>();
        List<double[]> readsRegList = new ArrayList<double[]>();
        DoubleArrayList genePValuesForQQPlot = new DoubleArrayList();
        int colNum = geneTableHeadRow.size();
        String nsVar = null, sVar = null;
        int scoreNum = 0;
        //the gene scores as covariables
        boolean hasGeneScore = false;
        if (!geneScores.isEmpty()) {
            scoreNum = geneScores.get((new ArrayList(geneScores.keySet())).get(0)).length;
            hasGeneScore = true;
        }
        scoreNum += 3;
        int nonZeroSynNum = 0;
        double nSN = 0;
        String COSMICCancerInfo;
        double[] geneScore = null;
        List<String> geneOrder = new ArrayList<String>();
        DoubleArrayList obsCounts = new DoubleArrayList();

        List<String[]> geneZeroMutRateSheet = new ArrayList<String[]>();
        boolean toTest = false;
        for (int i = 0; i < chroms.length; i++) {
            if (chroms[i] == null) {
                continue;
            }

            for (Gene gene : chroms[i].geneList) {
                if (gene.geneSymb == null) {
                    continue;
                }

                Double geneLen = geneLengths.get(gene.geneSymb);
                if (geneLen == null) {
                    continue;
                }
                if (hasGeneScore) {
                    geneScore = geneScores.get(gene.geneSymb);
                    if (geneScore == null) {
                        continue;
                    }
                }

                List<String> features = gene.featureValues;

                nsVar = features.get(somaticNSVarIndex);
                sVar = features.get(somaticSVarIndex);

                if (nsVar.equals(".")) {
                    continue;
                }
                nSN = Double.parseDouble(nsVar);
                if (nSN <= 0) {
                    if (toTest) {
                        String[] row = new String[colNum];
                        row[0] = gene.geneSymb;
                        row[1] = nsVar;
                        row[2] = sVar;
                        row[3] = geneLen.toString();
                        if (hasGeneScore) {
                            for (int t = 0; t < geneScore.length; t++) {
                                row[t + shiftCol + 4] = String.valueOf(geneScore[t]);
                            }
                        }
                        geneZeroMutRateSheet.add(row);
                    }
                    continue;
                }

                if (sVar.equals(".")) {
                    sVar = "0";
                }

                String[] row = new String[colNum];
                row[0] = gene.geneSymb;
                row[1] = nsVar;
                row[2] = sVar;
                row[3] = geneLen.toString();
                if (needPubMedAnno) {
                    //  geneTableRow.add(info[2]);
                }
                if (needCosmic) {
                    Map<String, Integer> sb = cosmicGeneMut.get(row[0]);
                    if (sb == null) {
                        COSMICCancerInfo = ".";
                    } else {
                        COSMICCancerInfo = sb.toString();
                    }
                    row[4] = COSMICCancerInfo;
                }

                double[] scoresM = new double[scoreNum];
                scoresM[0] = nSN;
                scoresM[1] = Double.parseDouble(sVar);
                scoresM[2] = geneLen;
                if (hasGeneScore) {
                    for (int t = 0; t < geneScore.length; t++) {
                        scoresM[3 + t] = (geneScore[t]);
                        row[t + shiftCol + 4] = String.valueOf(geneScore[t]);
                    }
                }
                obsCounts.add(nSN);
                countsRegList.add(scoresM);
                geneMutRateSheet.add(row);
                nonZeroSynNum += scoresM[1];

                if (regReads) {
                    double[] scoresS = new double[scoreNum];
                    scoresS[0] = Double.parseDouble(features.get(somaticNSRatioIndex));
                    scoresS[1] = Double.parseDouble(features.get(somaticSRatioIndex));
                    scoresS[2] = geneLen;
                    for (int t = 0; t < geneScore.length; t++) {
                        scoresS[3 + t] = (geneScore[t]);
                    }
                    readsRegList.add(scoresS);
                }
                geneOrder.add(gene.geneSymb);
            }
        }
        if (nonZeroSynNum < 5) {
            String info = "There seems no synomemous variants in your input data! The cancer gene analysis may have some problem!";
            LOG.error(info);
            return;
        }

        StringBuilder sb = new StringBuilder();
        double[] countsMatrix = new double[countsRegList.size() * countsRegList.get(0).length];
        for (int j = 0; j < countsRegList.size(); j++) {
            double[] v = countsRegList.get(j);
            System.arraycopy(v, 0, countsMatrix, j * v.length, v.length);
        }
        double[] readsMatrix = null;
        if (regReads) {
            readsMatrix = new double[readsRegList.size() * readsRegList.get(0).length];
            for (int j = 0; j < readsRegList.size(); j++) {
                double[] v = readsRegList.get(j);
                System.arraycopy(v, 0, readsMatrix, j * v.length, v.length);
            }
        }

        /*
        rcon.voidEval("qres.ztnb <- function(glm.obj)\n" + "{\n" + "#	Quantile residuals for Negative Binomial glm\n" + "#	GKS  22 Jun 97\n" + "#\n" + "	y <- glm.obj$y\n"
                + "	size <- glm.obj$theta\n" + "	mu <- fitted(glm.obj)\n" + "	p <- size/(mu + size)\n" + "	a <- ifelse(y > 1, pbeta(1-p, y+1, size)/(1-p^size), 0)\n"
                + "	b <- pbeta(1-p, y, size)/(1-p^size)\n" + "	u <- runif(n = length(y), min = a, max = b)\n" + "	qnorm(u)\n" + "}");
         */
        // update the codes in the orginal package with bugs
        rcon.voidEval("summary2.zerotrunc <- function(object,...)\n" + "{\n" + "  ## pearson residuals\n" + "  object$residuals <- residuals(object, type = \"pearson\")\n"
                + "\n" + "  ## compute z statistics\n" + "  cf <- object$coefficients\n" + "  se <- sqrt(diag(object$vcov))\n" + "  k <- length(cf)\n" + "  \n"
                + "  if(object$dist == \"negbin\") {\n" + "    cf <- c(cf, \"Log(theta)\" = as.vector(log(object$theta)))\n" + "    se <- c(se, object$SE.logtheta)\n" + "  }\n"
                + "  zstat <- cf/se\n" + "  pval <- 2*pnorm(-abs(zstat))\n" + "  cf <- cbind(cf, se, zstat, pval)\n"
                + "  colnames(cf) <- c(\"Estimate\", \"Std. Error\", \"z value\", \"Pr(>|z|)\")\n" + "  object$coefficients <- cf\n" + "\n" + "  ## number of iterations\n"
                + "  object$iterations <- tail(na.omit(object$optim$count), 1)\n" + "  \n" + "  ## delete some slots\n"
                + "  object$fitted.values <- object$terms <- object$model <- object$y <-\n" + "    object$x <- object$levels <- object$contrasts <- object$start <- NULL\n" + "\n"
                + "  ## return\n" + "  class(object) <- \"summary.zerotrunc\"\n" + "  object\n" + "}");

        String strPrefix = "genemutetest" + (int) (Math.random() * 10000);
        rcon.assign("valMat", countsMatrix);
        rcon.voidEval("valMat<-matrix(valMat,nrow=" + countsRegList.size() + ",ncol=" + countsRegList.get(0).length + ", byrow = TRUE)");
        rcon.voidEval(strPrefix + "<-data.frame(valMat);");
        // rcon.voidEval("m1 <- zerotrunc(valMat[,1] ~ valMat[,2] + valMat[,3] + valMat[,4] +valMat[,5] +valMat[,6], dist=\"negbin\")");
        sb.delete(0, sb.length());
        sb.append("colnames(").append(strPrefix).append(") <- c(\"DepVar\",\"IndepVar\",\"AvgCodingLen\"");
        for (String name : scoreHeads) {
            sb.append(",\"").append(name).append("\"");
        }
        sb.append(")");

        rcon.voidEval(sb.toString());

        sb.delete(0, sb.length());
        sb.append("m1 <- zerotrunc(DepVar ~ IndepVar + AvgCodingLen");
        for (String name : scoreHeads) {
            sb.append("+").append(name);
        }
        sb.append(", data = ").append(strPrefix).append(", dist=\"negbin\")");
        rcon.voidEval(sb.toString());

        rcon.voidEval("sumR<-summary2.zerotrunc(m1)");

        /*
        System.out.println("Running mutation rate test");
        File file = new File(geneSumOutFile + "summary.txt");
        String infor = "sink('" + file.getCanonicalPath() + "')";
        System.out.println(infor);
        // rcon.voidEval("sink('" + file.getCanonicalPath() + "')");
        System.out.println("Running mutation rate test");
        rcon.voidEval("print(sumR)");
        rcon.voidEval("sink()");
        //rcon.voidEval("save(textSR " + ", file=" + file.getCanonicalPath() + ")");
        String summary = rcon.eval("readChar(file('" + file.getCanonicalPath() + "', 'r'), 200000)").asString();
        LOG.info(summary);
         */
        double[] residueCounter = rcon.eval("resid(m1,\"pearson\")").asDoubles();

        // double[] residueCounter = rcon.eval("m1$fitted.values").asDoubles();

        /*
        
         double[] obsCountArray = new double[obsCounts.size()];
         for (int t = 0; t < obsCountArray.length; t++) {
         obsCountArray[t] = obsCounts.getQuick(t);
         System.out.println(fitCounter[t] + "\t" + obsCountArray[t]);
         }
         rcon.assign("obs", obsCountArray);
         double r = rcon.eval("cor(obs, m1$fitted.value, method = \"spearman\"").asDouble();
         r = rcon.eval("cor(mydf$DepVar, m1$fitted.values, method = \"spearman\"").asDouble();
         String info = "Correlation between observed and fitted values: " + r;
         LOG.info(summary + "\n" + info);
         */
        double[] residueReads = null;
        if (regReads) {
            rcon.assign("valMat", readsMatrix);
            rcon.voidEval("valMat<-matrix(valMat,nrow=" + readsRegList.size() + ",ncol=" + readsRegList.get(0).length + ", byrow = TRUE)");
            rcon.voidEval(strPrefix + "<-data.frame(valMat);");
            rcon.voidEval("colnames(" + strPrefix + ") <- c(\"DepVarReads\",\"IndepVarReads\",\"AvgCodingLen\",\"expr\",\"reptime\",\"hic\",\"constraint_score\")");
            rcon.voidEval("fit <- lm(DepVarReads ~ IndepVarReads +AvgCodingLen + expr +reptime +hic+constraint_score,data = " + strPrefix + ")");
            rcon.voidEval("res <- residuals(fit)");
            rcon.voidEval("re1<-sort(res)");
            rcon.voidEval("lower <- c(re1[length(re1)/20])");
            rcon.voidEval("upper <- c(re1[length(re1)*19/20])");
            rcon.voidEval("re1<-re1[re1>=lower&re1<=upper]");
            rcon.voidEval("re1<-matrix(re1,length(re1),1)");
            rcon.voidEval("fit1 <- mle.tmvnorm(re1, lower=lower, upper=upper)");
            rcon.voidEval("tt<-summary(fit1)");
            rcon.voidEval("mu <- tt@coef[1,1]");
            rcon.voidEval("var <- tt@coef[2,1]");
            // rcon.voidEval("res<-(res-mu)/sqrt(gene)");
            residueReads = rcon.eval("(res-mu)/sqrt(var)").asDoubles();
        }

        if (rcon != null) {
            rcon.close();
        }

        Map<String, Double> countsRegGeneP = new HashMap<String, Double>();
        Map<String, Double> lmRegGeneP = new HashMap<String, Double>();
        int geneNum = geneOrder.size();
        int popuSize = 19061;
        if (geneNum == geneMutRateSheet.size()) {
            for (int t = 0; t < geneNum; t++) {
                countsRegGeneP.put(geneOrder.get(t), residueCounter[t]);
                if (regReads) {
                    lmRegGeneP.put(geneOrder.get(t), residueReads[t]);
                }
                // System.out.println(res[t]);
            }

            List<String[]> geneMutRateSheet1 = new ArrayList<String[]>();
            List<String[]> geneMutRateSheet2 = new ArrayList<String[]>();
            boolean hasP = false;

            for (String[] v : geneMutRateSheet) {
                Double p = countsRegGeneP.get(v[0]);
                hasP = false;
                if (p != null) {
                    // System.out.println(p);
                    v[colNum - 2] = p.toString();
                    // p = (p - mean1) / sd1;
                    double[] pvalues = new double[1];
                    if (p > 0) {
                        p = Probability.normal(-p);
                    } else {
                        p = 1 - Probability.normal(p);
                    }
                    v[colNum - 1] = p.toString();
                    //  geneMutRegPValuePosVar.put(v[0], p);
                    genePValuesForQQPlot.add(p);
                    hasP = true;
                    pvalues[0] = p;
                    geneMutRegPValueAllVar.put(v[0], pvalues);
                }

                if (regReads) {
                    p = lmRegGeneP.get(v[0]);
                    if (p != null) {
                        // System.out.println(p);
                        v[colNum - 2] = p.toString();
                        // p = (p - mean1) / sd1;
                        if (p > 0) {
                            p = Probability.normal(-p);
                        } else {
                            p = 1 - Probability.normal(p);
                        }
                        v[colNum - 1] = p.toString();
                        hasP = true;
                        // genePValues.add(Double.parseDouble(v[5]));

                        // System.out.println(p);
                    }
                }
                if (hasP) {
                    geneMutRateSheet1.add(v);
                } else {
                    geneMutRateSheet2.add(v);
                }
            }

            Collections.sort(geneMutRateSheet1, new StringArrayDoubleComparator(colNum - 1));

            geneMutRateSheet.clear();
            geneMutRateSheet.addAll(geneMutRateSheet1);
            geneMutRateSheet.addAll(geneMutRateSheet2);
            if (toTest) {
                geneMutRateSheet.addAll(geneZeroMutRateSheet);
                Set<String> geneSymbs = new HashSet<String>();
                for (String[] items : geneMutRateSheet) {
                    geneSymbs.add(items[0]);
                }
                for (Map.Entry<String, double[]> geneP : geneScores.entrySet()) {
                    if (!geneSymbs.contains(geneP.getKey())) {
                        Double geneLen = geneLengths.get(geneP.getKey());
                        if (geneLen == null) {
                            continue;
                        }
                        geneScore = geneP.getValue();
                        String[] row = new String[colNum];
                        row[0] = geneP.getKey();
                        row[1] = "0";
                        row[2] = "0";
                        row[3] = geneLen.toString();
                        if (hasGeneScore) {
                            for (int t = 0; t < geneScore.length; t++) {
                                row[t + shiftCol + 4] = String.valueOf(geneScore[t]);
                            }
                        }
                        geneMutRateSheet.add(row);
                    }
                }

            }
            genePValuesForQQPlot.quickSort();
            List<DoubleArrayList> pvalueLists = new ArrayList<DoubleArrayList>();
            List<String> nameList = new ArrayList<String>();
            pvalueLists.add(genePValuesForQQPlot);
            PValuePainter pvPainter = new PValuePainter(450, 450);
            File plotFile = new File(geneSumOutFile.substring(0, geneSumOutFile.length() - 5) + "gene.qq.png");
            nameList.add("Somat");
            pvPainter.drawMultipleQQPlot(pvalueLists, nameList, null, plotFile.getCanonicalPath(), 1E-10);
            String info = "The QQ plot saved in " + plotFile.getCanonicalPath();
            LOG.info(info);
        } else {
            String msg = "Error, the hurdle regression p-values are not correctly estimate due to missing data";
            LOG.error(msg);
        }

        double adjGenePValueCutoff = MultipleTestingMethod.benjaminiHochbergFDR(genePValueCutoff, genePValuesForQQPlot);

        genePValuesForQQPlot.clear();
        StringBuilder sumSb = new StringBuilder();

        Set<String> wholeSigCandiGeneSet = new HashSet<String>();
        for (Map.Entry<String, double[]> geneP : geneMutRegPValueAllVar.entrySet()) {
            if (geneP.getValue()[0] <= adjGenePValueCutoff) {
                wholeSigCandiGeneSet.add(geneP.getKey());
            }
        }

        sumSb.append(wholeSigCandiGeneSet.size()).append(" genes with p-value <= ").append(adjGenePValueCutoff).append(" pass the Benjamini Hochberg FDR q value cutoff ")
                .append(genePValueCutoff).append(".\n");
        LOG.info(sumSb.toString());
        wholeSigCandiGeneSet.clear();

        geneMutRateSheet.add(0, geneTableHeadRow.toArray(new String[0]));
        LocalExcelFile.writeArray2XLSXFile(geneSumOutFile, geneMutRateSheet, true, -1, 0);
        String info = "The results of gene-based mutation rate test are saved in " + (new File(geneSumOutFile)).getCanonicalPath() + "!";

        LOG.info(info);
    }

}
