// (c) 2009-2011 Miaoxin Li
// This file is distributed as part of the KGG source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
// Permission is granted for you to use this file to compile IGG.
// All computer programs have bugs. Use this file at your own risk.
// Tuesday, March 01, 2011
/**
 *
 * @author mxli
 */
package org.cobi.kggseq.controller;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.log4j.Logger;
import org.cobi.kggseq.entity.PPIGraph;
import org.cobi.kggseq.entity.GeneSet;

import org.cobi.util.text.Util;

/**
 * 
 * @author mxli
 */
public class CandidateGeneExtender {

    private Set<String> seedGeneSet;
    private Set<String> otherInterestingGeneSymbols;
    private Set<String> ppiExtendedGeneSet;
    private Set<String> pathwayExtendedGeneSet;
    private PPIGraph ppiNetwork;
    private GeneSetExplorer genesetExplorer;
    // private static final Log LOG = Log.getInstance(CandidateGeneExtender.class);
    private static final Logger LOG = Logger.getLogger(CandidateGeneExtender.class);

    public String countCoverages(List<String> gwasGeneList, int ppiLevel) {
        List<String[]> forwardInteractionItems = ppiNetwork.getForwardInteractionItems();
        HashSet<String> ppiGeneSet = new HashSet<String>();
        int geneSize = forwardInteractionItems.size();
        for (int i = 0; i < geneSize; i++) {
            ppiGeneSet.add(forwardInteractionItems.get(i)[0]);
            ppiGeneSet.add(forwardInteractionItems.get(i)[1]);
        }

        Set<String> geneSet = new HashSet<String>();
        Map<String, GeneSet> genesets = genesetExplorer.getGeneSetSet();

        for (Map.Entry<String, GeneSet> mPathway : genesets.entrySet()) {
            HashSet<String> pathwayGenes = mPathway.getValue().getGeneSymbols();
            geneSet.addAll(pathwayGenes);
        }

        geneSize = gwasGeneList.size();
        double geneNumPPI = 0;
        double geneNumPathway = 0;
        double geneNumAll = 0;
        double geneNumEither = 0;
        for (int i = 0; i < geneSize; i++) {
            String geneSymb = gwasGeneList.get(i);
            if (ppiGeneSet.contains(geneSymb)) {
                geneNumPPI += 1;
            }
            if (geneSet.contains(geneSymb)) {
                geneNumPathway += 1;
            }
            if (ppiGeneSet.contains(geneSymb) || geneSet.contains(geneSymb)) {
                geneNumEither += 1;
            } else if (ppiGeneSet.contains(geneSymb) && geneSet.contains(geneSymb)) {
                geneNumAll += 1;
            }
        }
        StringBuilder infor = new StringBuilder();
        infor.append("There are ").append(geneSize).append(" genes in your GWAS dataset; ");
        infor.append((int) geneNumPPI);
        infor.append("(");
        double percentage = 100 * geneNumPPI / geneSize;
        infor.append(Util.doubleToString(percentage, 3));
        infor.append("%) genes have registry in our PPI dataset; ");
        infor.append((int) geneNumPathway);
        infor.append("(");
        percentage = 100 * geneNumPathway / geneSize;
        infor.append(Util.doubleToString(percentage, 3));
        infor.append("%) genes have registry in our pathway dataset.\n    ");
        infor.append((int) geneNumEither);
        infor.append("(");
        percentage = 100 * geneNumEither / geneSize;
        infor.append(Util.doubleToString(percentage, 3));
        infor.append("%) genes have registry in either of datasets; ");
        infor.append((int) geneNumAll);
        infor.append("(");
        percentage = 100 * geneNumAll / geneSize;
        infor.append(Util.doubleToString(percentage, 3));
        infor.append("%) genes have registry in both datasets simultaneously.\n");

        geneNumPPI = 0;
        geneNumAll = 0;
        geneNumPathway = 0;
        geneNumEither = 0;
        for (int i = 0; i < geneSize; i++) {
            String geneSymb = gwasGeneList.get(i);
            if (ppiExtendedGeneSet.contains(geneSymb)) {
                geneNumPPI += 1;
            }
            if (pathwayExtendedGeneSet.contains(geneSymb)) {
                geneNumPathway += 1;
            }
            if (ppiExtendedGeneSet.contains(geneSymb) || pathwayExtendedGeneSet.contains(geneSymb)) {
                geneNumEither += 1;
            } else if (ppiExtendedGeneSet.contains(geneSymb) && pathwayExtendedGeneSet.contains(geneSymb)) {
                geneNumAll += 1;
            }
        }

        infor.append((int) geneNumPPI);
        infor.append("(");
        percentage = 100 * geneNumPPI / geneSize;
        infor.append(Util.doubleToString(percentage, 3));
        infor.append("%) genes have ").append(ppiLevel).append(" - (or less-) level PPIs with seed candidate genes;");
        infor.append((int) geneNumPathway);
        infor.append("(");
        percentage = 100 * geneNumPathway / geneSize;
        infor.append(Util.doubleToString(percentage, 3));
        infor.append("%) genes share same pathways with seed candidate genes.\n    ");
        infor.append((int) geneNumEither);
        infor.append("(");
        percentage = 100 * geneNumEither / geneSize;
        infor.append(Util.doubleToString(percentage, 3));
        infor.append("%) genes meet either of conditions; ");
        infor.append((int) geneNumAll);
        infor.append("(");
        percentage = 100 * geneNumAll / geneSize;
        infor.append(Util.doubleToString(percentage, 3));
        infor.append("%) genes meet both of conditions.");

        return infor.toString();
    }

    public Set<String> getSeedGeneSet() {
        return seedGeneSet;
    }

    public void setSeedGeneSet(Set<String> seedGeneSet) {
        this.seedGeneSet = seedGeneSet;
    }

    /**
     *
     * @throws Exception
     */
    public CandidateGeneExtender() throws Exception {
        seedGeneSet = new HashSet<String>();
        otherInterestingGeneSymbols = new HashSet<String>();
        ppiExtendedGeneSet = new HashSet<String>();
        pathwayExtendedGeneSet = new HashSet<String>();
    }

    /**
     *
     * @return
     * @throws Exception
     */
    public final HashSet<String> getAllGeneSymbs() throws Exception {
        HashSet<String> geneSet = new HashSet<String>();
        geneSet.addAll(seedGeneSet);
        geneSet.addAll(ppiExtendedGeneSet);
        geneSet.addAll(pathwayExtendedGeneSet);
        geneSet.addAll(otherInterestingGeneSymbols);
        return geneSet;
    }

    /**
     *
     * @param fileName
     * @throws Exception
     */
    public final void readSeedGenesFromFile(final String fileName)
            throws Exception {
        String line = null;
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        StringBuilder tmpBuf = new StringBuilder();
        try {
            while ((line = br.readLine()) != null) {
                StringTokenizer st = new StringTokenizer(line, "\t");
                // System.out.println(line);
                // source
                st.nextToken();
                tmpBuf.append(st.nextToken().trim());
                String symbol = tmpBuf.toString();
                tmpBuf.delete(0, tmpBuf.length());
                st.nextToken();
                st.nextToken();
                st.nextToken();
                String asSeed = st.nextToken().trim();
                if (asSeed.equals("true")) {
                    seedGeneSet.add(symbol);
                } else {
                    otherInterestingGeneSymbols.add(symbol);
                }
            }
            StringBuilder inforString = new StringBuilder();
            inforString.append(seedGeneSet.size());
            inforString.append(" seed candidate genes, ");
            inforString.append(otherInterestingGeneSymbols.size());
            inforString.append(" non-seed candidate genes are read.");
            LOG.info(inforString);
        } finally {
            br.close();
        }
    }

    /**
     *
     * @throws Exception
     */
    public final void extendPathwayGenes() throws Exception {
        Map<String, GeneSet> allPathways = genesetExplorer.getGeneSetSet();
        pathwayExtendedGeneSet.clear();
        for (Map.Entry<String, GeneSet> mPathway : allPathways.entrySet()) {
            HashSet<String> pathwayGenes = mPathway.getValue().getGeneSymbols();
            Iterator<String> itSeed = seedGeneSet.iterator();
            while (itSeed.hasNext()) {
                if (pathwayGenes.contains(itSeed.next())) {
                    pathwayExtendedGeneSet.addAll(pathwayGenes);
                    break;
                }
            }
        }

    }

    /**
     *
     * @param level
     * @throws Exception
     */
    public final void extendPPIGenes(final int level) throws Exception {
        ppiExtendedGeneSet.clear();
        Iterator<String> itSeed = seedGeneSet.iterator();
        Graph<String, Integer> ppiGraph = new UndirectedSparseGraph<String, Integer>();

        while (itSeed.hasNext()) {
            String geneSymb = itSeed.next();
            ppiGraph.addVertex(geneSymb);
            ppiNetwork.buildGraph(ppiGraph, geneSymb, 0, level);
        }

        Collection<String> genes = ppiGraph.getVertices();
        Iterator<String> ite = genes.iterator();
        String geneSyb = null;
        while (ite.hasNext()) {
            geneSyb = ite.next();
            if (!seedGeneSet.contains(geneSyb)) {
                ppiExtendedGeneSet.add(geneSyb);
            }
        }
    }

    /**
     *
     * @param level
     * @throws Exception
     */
    public final void loadPPIDB(String ppidbFile) throws Exception {
        // build ppi network
        ppiExtendedGeneSet = new HashSet<String>();
        ppiNetwork = new PPIGraph(ppidbFile);
        ppiNetwork.readPPIItems();
    }

    public final Graph<String, Integer> pickRelevantPPIGrapph(final int level) throws Exception {
        Iterator<String> itSeed = seedGeneSet.iterator();
        Graph<String, Integer> ppiGraph = new UndirectedSparseGraph<String, Integer>();

        while (itSeed.hasNext()) {
            String geneSymb = itSeed.next();
            ppiGraph.addVertex(geneSymb);
            ppiNetwork.buildGraph(ppiGraph, geneSymb, 0, level);
        }
        return ppiGraph;
    }

    public final void loadGeneSetDB(String pathwayFile, int minPathwayGene, int maxPathwayGene) throws Exception {
        pathwayExtendedGeneSet = new HashSet<String>();
        genesetExplorer = new GeneSetExplorer(pathwayFile);
        genesetExplorer.loadGSEAGeneSets(minPathwayGene, maxPathwayGene);
    }

    public final Map<String, GeneSet> pickRelevantGeneSet() throws Exception {
        Map<String, GeneSet> selPathways = new HashMap<String, GeneSet>();
        Map<String, GeneSet> allPathways = genesetExplorer.getGeneSetSet();

        for (Map.Entry<String, GeneSet> mPathway : allPathways.entrySet()) {
            HashSet<String> pathwayGenes = mPathway.getValue().getGeneSymbols();
            Iterator<String> itSeed = seedGeneSet.iterator();
            while (itSeed.hasNext()) {
                if (pathwayGenes.contains(itSeed.next())) {
                    selPathways.put(mPathway.getKey(), mPathway.getValue());
                    break;
                }
            }
        }
        return selPathways;
    }

    public final Map<String, GeneSet> pickAllPathways() throws Exception {
        Map<String, GeneSet> allPathways = genesetExplorer.getGeneSetSet();
        return allPathways;
    }

    /**
     *
     * @param seedSet
     * @throws Exception
     */
    public final void replaceSeedGenes(final HashSet<String> seedSet)
            throws Exception {
        seedGeneSet.clear();
        seedGeneSet.addAll(seedSet);
    }
}
