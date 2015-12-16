// (c) 2009-2011 Miaoxin Li
// This file is distributed as part of the KGG source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
// Permission is granted for you to use this file to compile IGG.
// All computer programs have bugs. Use this file at your own risk.
// Tuesday, March 01, 2011
package org.cobi.kggseq.controller;

import cern.colt.function.DoubleProcedure;
import cern.colt.list.DoubleArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.HashSet;
import java.util.Iterator;

import org.apache.log4j.Logger;
import org.cobi.kggseq.Constants;
import org.cobi.kggseq.entity.GeneSet;
import org.cobi.util.file.LocalFileFunc;
import org.cobi.util.stat.MultipleTestingMethod;

/**
 *
 * @author Miaoxin Li
 */
public class GeneSetExplorer implements Constants {

    // private static final Log LOG = Log.getInstance(GeneSetExplorer.class);
    private static final Logger LOG = Logger.getLogger(GeneSetExplorer.class);
    private Map<String, GeneSet> geneSetSet;
    private int uniqueTotalPathwayGeneNum = 0;
    String genesetFileName = null;

    public int getUniqueTotalPathwayGeneNum() {
        return uniqueTotalPathwayGeneNum;
    }

    public void setUniqueTotalPathwayGeneNum(int uniqueTotalPathwayGeneNum) {
        this.uniqueTotalPathwayGeneNum = uniqueTotalPathwayGeneNum;
    }

    /**
     *
     * @throws Exception
     */
    public GeneSetExplorer(String pathFile) throws Exception {
        genesetFileName = pathFile;
        geneSetSet = new HashMap<String, GeneSet>();
    }

    public Map<String, GeneSet> getGeneSetSet() {
        return geneSetSet;
    }

    public void setGeneSetSet(Map<String, GeneSet> geneSetSet) {
        this.geneSetSet = geneSetSet;
    }

    /**
     *
     * @param kggPathwayFileName
     * @param biocartaPathwayFileName
     * @throws Exception
     */
    public void loadGSEAGeneSets(int minGene, int maxGene) throws Exception {
        String line = null;

        geneSetSet.clear();
        HashSet<String> smalllPathIDToRemove = new HashSet<String>();
        HashSet<String> largePathIDToRemove = new HashSet<String>();
        HashSet<String> allPathwayGenes = new HashSet<String>();


        File dataFile = new File(genesetFileName);
//        LineReader br = null;
        BufferedReader br=null;
        if (dataFile.exists() && dataFile.getName().endsWith(".zip")) {
//            br = new CompressedFileReader(dataFile);
            br=LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
        } else {
            if (dataFile.exists() && dataFile.getName().endsWith(".gz")) {
//                br = new CompressedFileReader(dataFile);
                br=LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
            } else {
                if (dataFile.exists()) {
//                    br = new AsciiLineReader(dataFile);
                    br=LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
                } else {
                    throw new Exception("No input file: " + dataFile.getCanonicalPath());
                }
            }
        }

        while ((line = br.readLine()) != null) {
            StringTokenizer st = new StringTokenizer(line);
            //skip gene ID
            String pathID = st.nextToken().trim();
            String url = st.nextToken().trim();

            GeneSet geneset = geneSetSet.get(pathID);
            if (geneset == null) {
                geneset = new GeneSet(pathID, pathID, url);
            }
            while (st.hasMoreTokens()) {
                String geneSymb = st.nextToken().trim();
                geneset.addGeneSymbol(geneSymb);
            }

            if (geneset.getGeneSymbols().size() < minGene) {
                smalllPathIDToRemove.add(geneset.getID());
            } else if (geneset.getGeneSymbols().size() > maxGene) {
                largePathIDToRemove.add(geneset.getID());
            } else {
                geneSetSet.put(pathID, geneset);
                allPathwayGenes.addAll(geneset.getGeneSymbols());
            }
        }
        br.close();
        allPathwayGenes.clear();
        uniqueTotalPathwayGeneNum = allPathwayGenes.size();

        /*
        StringBuilder info = new StringBuilder();
        Iterator<String> iter;
        if (largePathIDToRemove.size() > 0) {
        info.append("The following pathways are excluded due to their gene number is over ");
        info.append(maxGene);
        info.append("\n: [");
        iter = largePathIDToRemove.iterator();
        while (iter.hasNext()) {
        info.append(iter.next());
        info.append("; ");
        }
        largePathIDToRemove.clear();
        info.append("]\n\n");
        }
        
        
        if (smalllPathIDToRemove.size() > 0) {
        info.append("The following pathways are excluded due to their gene number is less than ");
        info.append(minGene);
        info.append("\n: [");
        iter = smalllPathIDToRemove.iterator();
        while (iter.hasNext()) {
        String id = iter.next();
        info.append(id);
        info.append("; ");
        }
        largePathIDToRemove.clear();
        info.append("]\n");
        }
        
        LOG.info(info);
         */

    }

    /**
     *
     * @param totalGWASGenes
     * @param seedGenes
     * @param sigGenes
     * @param searchedPathwayList
     * @throws Exception
     */
    public void searchEnrichedPathwaysby2Sets(HashSet<String> totalGWASGenes, Map<String, String> seedGenes, Map<String, String> sigGenes, List<GeneSet> searchedPathwayList) throws Exception {
        int seedGeneInPathNum = 0;
        int sigGeneInPathNum = 0;
        int subPopulationSize = 0;

        Map<String, GeneSet> searchedPathways = new HashMap<String, GeneSet>();
        for (Map.Entry<String, GeneSet> mPath : geneSetSet.entrySet()) {
            HashSet<String> pGenes = mPath.getValue().getGeneSymbols();
            String pathID = mPath.getKey();
            GeneSet curPath = searchedPathways.get(pathID);
            if (curPath == null) {
                curPath = new GeneSet(pathID, mPath.getValue().getName(), mPath.getValue().getURL());
                curPath.setURL(mPath.getValue().getURL());
            }
            seedGeneInPathNum = 0;
            sigGeneInPathNum = 0;
            subPopulationSize = 0;
            Iterator<String> it = pGenes.iterator();
            while (it.hasNext()) {
                String geneSym = it.next();
                if (totalGWASGenes.contains(geneSym)) {
                    subPopulationSize++;
                }
                if (seedGenes.containsKey(geneSym)) {
                    curPath.addGeneSymbol(geneSym);
                    seedGeneInPathNum++;
                }
                if (sigGenes.containsKey(geneSym)) {
                    curPath.addGeneSymbol(geneSym);
                    sigGeneInPathNum++;
                }
            }

            /*
            if (seedGeneInPathNum > 0 && sigGeneInPathNum > 0) {
            curPath.setEnrichedPValue(MultipleTestingMethod.hypergeometricEnrichmentTest(GENOME_GENE_NUM, pGenes.size(), seedGenes.size(), seedGeneInPathNum)
             * MultipleTestingMethod.hypergeometricEnrichmentTest(totalGWASGenes.size(), subPopulationSize, sigGenes.size(), sigGeneInPathNum));
            searchedPathways.put(pathID, curPath);
            }
             * 
             */
        }
        searchedPathwayList.addAll(searchedPathways.values());
        //Collections.sort(searchedPathwayList, new PathwayPValueComparator());
    }

    /**
     *
     * @param totalGWASGenes
     * @param sigGenes
     * @param pathwayList
     * @throws Exception
     */
    public void searchEnrichedPathways(HashSet<String> totalGWASGenes, Map<String, String> sigGenes, List<GeneSet> pathwayList) throws Exception {
        Map<String, GeneSet> searchedPathways = new HashMap<String, GeneSet>();

        int sigGeneInPathNum = 0;
        int subPopulationSize = 0;
        for (Map.Entry<String, GeneSet> mPath : geneSetSet.entrySet()) {
            HashSet<String> pGenes = mPath.getValue().getGeneSymbols();
            String pathID = mPath.getKey();
            GeneSet curPath = searchedPathways.get(pathID);
            if (curPath == null) {
                curPath = new GeneSet(pathID, mPath.getValue().getName(), mPath.getValue().getURL());
                curPath.setURL(mPath.getValue().getURL());
            }

            sigGeneInPathNum = 0;
            subPopulationSize = 0;
            Iterator<String> it = pGenes.iterator();
            while (it.hasNext()) {
                String geneSym = it.next();
                if (totalGWASGenes.contains(geneSym)) {
                    subPopulationSize++;
                }
                if (sigGenes.containsKey(geneSym)) {
                    curPath.addGeneSymbol(geneSym);
                    sigGeneInPathNum++;
                }
            }

            if (sigGeneInPathNum > 0) {
                curPath.setEnrichedPValue(MultipleTestingMethod.hypergeometricEnrichmentTest(totalGWASGenes.size(),
                        sigGenes.size(), subPopulationSize, sigGeneInPathNum));
                searchedPathways.put(pathID, curPath);
            }
        }
        pathwayList.addAll(searchedPathways.values());
        //  Collections.sort(pathwayList, new PathwayPValueComparator());
    }

    /**
     *
     * @param allGWASGenes
     * @param selectedGWASGenes
     * @param pathwayList
     * @throws Exception
     */
    public void searchEnrichedPathwaysby2Sets(HashSet<String> allGWASGenes,
            HashSet<String> selectedGWASGenes, HashSet<String> seedGenes, List<GeneSet> pathwayList) throws Exception {
        //note the selectedGWASGenes belongs to the allGWASGenes
        Map<String, GeneSet> searchedPathways = new HashMap<String, GeneSet>();
        int sigGeneInPathNum = 0;
        int subPopulationSize = 0;
        int allGWASGeneNum = allGWASGenes.size();
        int selectedGWASGeneNum = selectedGWASGenes.size();

        for (Map.Entry<String, GeneSet> mPath : geneSetSet.entrySet()) {
            HashSet<String> pGenes = mPath.getValue().getGeneSymbols();
            String pathID = mPath.getKey();
            GeneSet curPath = searchedPathways.get(pathID);
            if (curPath == null) {
                curPath = new GeneSet(pathID, mPath.getValue().getName(), mPath.getValue().getURL());
                curPath.setURL(mPath.getValue().getURL());
            }

            sigGeneInPathNum = 0;
            subPopulationSize = 0;
            Iterator<String> pGeneIter = pGenes.iterator();
            if (seedGenes != null && seedGenes.size() > 0) {
                while (pGeneIter.hasNext()) {
                    String pGeneSymb = pGeneIter.next();
                    if (allGWASGenes.contains(pGeneSymb)) {
                        if (selectedGWASGenes.contains(pGeneSymb)) {
                            curPath.addGeneSymbol(pGeneSymb);
                            sigGeneInPathNum++;
                        }
                        //include the seed gene into the geneset anyway
                        if (seedGenes.contains(pGeneSymb)) {
                            curPath.addGeneSymbol(pGeneSymb);
                            sigGeneInPathNum++;
                        }
                        subPopulationSize++;
                    }
                }
            } else {
                while (pGeneIter.hasNext()) {
                    String pGeneSymb = pGeneIter.next();
                    if (allGWASGenes.contains(pGeneSymb)) {
                        if (selectedGWASGenes.contains(pGeneSymb)) {
                            curPath.addGeneSymbol(pGeneSymb);
                            sigGeneInPathNum++;
                        }
                        subPopulationSize++;
                    }
                }
            }
            curPath.setTotalGeneNum(subPopulationSize);
            if (sigGeneInPathNum > 0) {
                curPath.setEnrichedPValue(MultipleTestingMethod.hypergeometricEnrichmentTest(allGWASGeneNum,
                        selectedGWASGeneNum, subPopulationSize, sigGeneInPathNum));
                searchedPathways.put(pathID, curPath);
            }
        }
        pathwayList.addAll(searchedPathways.values());
    }

    /**
     *
     * @param genePValues
     * @param pathwayPThreshold
     * @param genePThreshold
     * @param pathwayList
     * @return
     * @throws Exception
     */
    public DoubleArrayList searchSignificantPathwaysSimesTest(Map<String, Double> genePValues,
            double presentationPathwayPThreshold, double presetationGenePThreshold, List<GeneSet> pathwayList) throws Exception {
        Map<String, GeneSet> searchedPathways = new HashMap<String, GeneSet>();
        DoubleArrayList pValues = new DoubleArrayList();
        int geneNum;
        double pMin;
        DoubleArrayList pathwayPValues = new DoubleArrayList();
        HashSet<String> genesKept = new HashSet<String>();
        for (Map.Entry<String, GeneSet> mPath : geneSetSet.entrySet()) {
            GeneSet curPath = mPath.getValue();
            String pathID = mPath.getKey();
            pValues.clear();
            genesKept.clear();
            HashSet<String> pGenes = curPath.getGeneSymbols();
            Iterator<String> pGeneIter = pGenes.iterator();
            while (pGeneIter.hasNext()) {
                String pGeneSymb = pGeneIter.next();
                Double pV = genePValues.get(pGeneSymb);
                if (pV != null) {
                    //note all avaible genes are involved in the test
                    pValues.add(pV);
                    if (pV <= presetationGenePThreshold) {
                        genesKept.add(pGeneSymb);
                    }
                }
            }
            pValues.quickSort();
            //use Simes test to combine p-values
            geneNum = pValues.size();
            if (geneNum == 0) {
                continue;
            }

            GeneSet newPathway = new GeneSet(pathID, pathID, curPath.getURL());
            newPathway.getGeneSymbols().addAll(genesKept);

            pMin = pValues.get(0) * geneNum;
            if (geneNum > 1) {
                for (int i = 1; i < geneNum; i++) {
                    if ((pValues.getQuick(i) * geneNum / (i + 1)) < pMin) {
                        pMin = pValues.getQuick(i) * geneNum / (i + 1);
                    }
                }
            }
            newPathway.setTotalGeneNum(geneNum);
            pathwayPValues.add(pMin);
            if (pMin <= presentationPathwayPThreshold) {
                newPathway.setEnrichedPValue(pMin);
                searchedPathways.put(pathID, newPathway);
            }
        }
        pathwayList.addAll(searchedPathways.values());
        return pathwayPValues;
    }

    /**
     *
     * @param args
     */
    public static void main(String[] args) {
        int size = 100;
        int[] buff = new int[size];
        long start = System.currentTimeMillis();
        DoubleArrayList pValues = new DoubleArrayList();
        for (int i = 0; i < size; i++) {
            buff[i] = 1;
            pValues.add(1);
        }

        System.out.println(System.currentTimeMillis() - start);

        start = System.currentTimeMillis();
        for (int i = size - 1; i <= 0; i--) {
            buff[i] = 1;
        }
        System.out.println(System.currentTimeMillis() - start);
        try {
            //calcluate pvalues
            int allGeneSize = 15300;
            DoubleProcedure dp = new DoubleProcedure() {

                @Override
                public boolean apply(double a) {
                    System.out.print(a);
                    return true;
                }
            };


        } catch (Exception ex) {
        }
    }
}
