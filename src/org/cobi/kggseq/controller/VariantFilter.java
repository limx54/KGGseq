/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import cern.colt.list.FloatArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.map.OpenLongObjectHashMap;
import cern.jet.stat.Probability;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.log4j.Logger;
import org.cobi.kggseq.GlobalManager;
import org.cobi.kggseq.entity.AnnotationSummarySet;
import org.cobi.kggseq.entity.Chromosome;
import org.cobi.kggseq.entity.FiltrationSummarySet;
import org.cobi.kggseq.entity.Genome;
import org.cobi.kggseq.entity.Individual;
import org.cobi.kggseq.entity.RNABoundaryIndex;
import org.cobi.kggseq.entity.RefDup;
import org.cobi.kggseq.entity.ReferenceGenome;
import org.cobi.kggseq.entity.Variant;
import org.cobi.util.net.NCBIRetriever;
import org.cobi.util.stat.ContingencyTable;
import org.cobi.util.text.Util;

/**
 *
 * @author mxli
 */
public class VariantFilter {

    // private static final Log LOG = Log.getInstance(VariantAnnotator.class);
    private static final Logger LOG = Logger.getLogger(VariantFilter.class);

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

    public boolean checkCompoundHeteroGeneSudo(IntArrayList effectiveIndivIDs, List<Individual> sortedSubjectList, int[] pedEncodeGytIDMap, List<int[]> triosIDList, List<Variant> geneDisVars,
            int chromID, boolean isPhasedGty, boolean noNeedHomo, String lastGeneSymb, Set<String> caseGenes, Set<String> controlGenes,
            List<Variant> tmpVarListGene, OpenLongObjectHashMap wahBit, int[] hitIndivCount1, Set<String> effectiveVarSet, String[] couts1, String[] couts2,
            boolean allPsudoControl, int fixedColNum) {

        int effectiveIndivSize = effectiveIndivIDs.size();
        StringBuilder varTransmitInfo = new StringBuilder();
        StringBuilder varUntransmitInfo = new StringBuilder();
        boolean isEffective;
        boolean isEffectiveUntransm;
        boolean hasADouble = false;
        StringBuilder varLabel = new StringBuilder();

        int caseGeneCounts = 0;
        int sudoControlGeneCounts = 0;
        List<Variant> tmpVarListIndiv = new ArrayList<Variant>();
        int[] childGty = null;
        int[] fatherGty = null;
        int[] motherGty = null;
        int base = 0;
        int alleleNum = 0;
        boolean[] bits = new boolean[32];
        int startIndex;
        for (int j = 0; j < effectiveIndivSize; j++) {
            int index = effectiveIndivIDs.getQuick(j);
            Individual mIndivChild = sortedSubjectList.get(triosIDList.get(index)[0]);

            int countN = 0;
            boolean transmitAltFromFa = false;
            boolean transmitAltFromMo = false;

            boolean transmitAltFromFaSudo = false;
            boolean transmitAltFromMoSudo = false;
            varTransmitInfo.delete(0, varTransmitInfo.length());
            varUntransmitInfo.delete(0, varUntransmitInfo.length());
            tmpVarListIndiv.clear();

            for (Variant tVar : geneDisVars) {
                isEffective = false;
                isEffectiveUntransm = false;

                alleleNum = tVar.getAltAlleles().length + 1;

                if (isPhasedGty) {
                    base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
                    if (tVar.compressedGty) {
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[0]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        childGty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[0]]);
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[1]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        fatherGty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[1]]);
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[2]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        motherGty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[2]]);
                    } else {
                        childGty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[0]], pedEncodeGytIDMap.length);
                        fatherGty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[1]], pedEncodeGytIDMap.length);
                        motherGty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[2]], pedEncodeGytIDMap.length);
                    }

                } else {
                    base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
                    if (tVar.compressedGty) {
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[0]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        childGty = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[0]]);
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[1]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        fatherGty = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[1]]);
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[2]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        motherGty = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[2]]);
                    } else {
                        childGty = BinaryGtyProcessor.getUnphasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[0]], pedEncodeGytIDMap.length);
                        fatherGty = BinaryGtyProcessor.getUnphasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[1]], pedEncodeGytIDMap.length);
                        motherGty = BinaryGtyProcessor.getUnphasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[2]], pedEncodeGytIDMap.length);
                    }

                }
                if (childGty == null || fatherGty == null || motherGty == null) {
                    continue;
                }
                if (noNeedHomo && childGty[0] == childGty[1]) {
                    continue;
                }

                // assume the reference allele has not causal effect.
                if (childGty[0] == 0 && childGty[1] != 0) {
                    // Warning: this will exclude the F:0/1 & M:0/1->O:0/1
                    // genotypes of parents
                    if ((childGty[1] == fatherGty[0] || childGty[1] == fatherGty[1]) && childGty[1] != motherGty[0] && childGty[1] != motherGty[1]) {
                        if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                            continue;
                        }
                        countN++;
                        transmitAltFromFa = true;
                        isEffective = true;
                    } else if ((childGty[1] == motherGty[0] || childGty[1] == motherGty[1]) && childGty[1] != fatherGty[0] && childGty[1] != fatherGty[1]) {
                        if (motherGty[0] != 0 && motherGty[1] != 0) {
                            continue;
                        }
                        countN++;
                        transmitAltFromMo = true;
                        isEffective = true;
                    }
                } else if (childGty[0] != 0 && childGty[1] != 0) {
                    // ignore variants with more than 2 alleles
                    if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                        continue;
                    }
                    if (motherGty[0] != 0 && motherGty[1] != 0) {
                        continue;
                    }

                    // then offspreing is homozygous
                    if ((childGty[0] == fatherGty[0] || childGty[0] == fatherGty[1]) && (childGty[1] == motherGty[0] || childGty[1] == motherGty[1])) {
                        countN++;
                        countN++;
                        transmitAltFromFa = true;
                        transmitAltFromMo = true;
                        isEffective = true;
                    } else if ((childGty[1] == fatherGty[0] || childGty[1] == fatherGty[1]) && (childGty[0] == motherGty[0] || childGty[0] == motherGty[1])) {
                        countN++;
                        countN++;
                        transmitAltFromFa = true;
                        transmitAltFromMo = true;
                        isEffective = true;
                    }
                } else if (childGty[0] != 0) {
                    // Warning: this will exclude the F:0/1 & M:0/1->O:0/1
                    // genotypes of parents
                    if ((childGty[0] == fatherGty[0] || childGty[0] == fatherGty[1]) && childGty[0] != motherGty[0] && childGty[0] != motherGty[1]) {
                        if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                            continue;
                        }
                        countN++;
                        transmitAltFromFa = true;
                        isEffective = true;
                    } else if ((childGty[0] == motherGty[0] || childGty[0] == motherGty[1]) && childGty[0] != fatherGty[0] && childGty[0] != fatherGty[1]) {
                        if (motherGty[0] != 0 && motherGty[1] != 0) {
                            continue;
                        }
                        countN++;
                        transmitAltFromMo = true;
                        isEffective = true;
                    }
                }

                int[] uchildGty = unTransmittedGty(fatherGty, motherGty, childGty);
                if (uchildGty != null) {
                    if (uchildGty[0] == 0 && uchildGty[1] != 0) {
                        // Warning: this will exclude the F:0/1 & M:0/1->O:0/1
                        // genotypes of parents
                        if ((uchildGty[1] == fatherGty[0] || uchildGty[1] == fatherGty[1]) && uchildGty[1] != motherGty[0] && uchildGty[1] != motherGty[1]) {
                            if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                                continue;
                            }
                            transmitAltFromFaSudo = true;
                            isEffectiveUntransm = true;
                        } else if ((uchildGty[1] == motherGty[0] || uchildGty[1] == motherGty[1]) && uchildGty[1] != fatherGty[0] && uchildGty[1] != fatherGty[1]) {
                            if (motherGty[0] != 0 && motherGty[1] != 0) {
                                continue;
                            }
                            transmitAltFromMoSudo = true;
                            isEffectiveUntransm = true;
                        }
                    } else if (uchildGty[0] != 0 && uchildGty[1] != 0) {
                        // ignore variants with more than 2 alleles
                        if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                            continue;
                        }
                        if (motherGty[0] != 0 && motherGty[1] != 0) {
                            continue;
                        }

                        // then offspreing is homozygous
                        if ((uchildGty[0] == fatherGty[0] || uchildGty[0] == fatherGty[1]) && (uchildGty[1] == motherGty[0] || uchildGty[1] == motherGty[1])) {
                            transmitAltFromFaSudo = true;
                            transmitAltFromMoSudo = true;
                            isEffectiveUntransm = true;
                        } else if ((uchildGty[1] == fatherGty[0] || uchildGty[1] == fatherGty[1]) && (uchildGty[0] == motherGty[0] || uchildGty[0] == motherGty[1])) {
                            transmitAltFromFaSudo = true;
                            transmitAltFromMoSudo = true;
                            isEffectiveUntransm = true;
                        }
                    } else if (uchildGty[0] != 0) {
                        // Warning: this will exclude the F:0/1 & M:0/1->O:0/1
                        // genotypes of parents
                        if ((uchildGty[0] == fatherGty[0] || uchildGty[0] == fatherGty[1]) && uchildGty[0] != motherGty[0] && uchildGty[0] != motherGty[1]) {
                            if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                                continue;
                            }
                            transmitAltFromFaSudo = true;
                            isEffectiveUntransm = true;
                        } else if ((uchildGty[0] == motherGty[0] || uchildGty[0] == motherGty[1]) && uchildGty[0] != fatherGty[0] && uchildGty[0] != fatherGty[1]) {
                            if (motherGty[0] != 0 && motherGty[1] != 0) {
                                continue;
                            }
                            transmitAltFromMoSudo = true;
                            isEffectiveUntransm = true;
                        }
                    }
                }

                if (isEffectiveUntransm) {
                    varUntransmitInfo.append(tVar.refStartPosition).append(":");
                    varUntransmitInfo.append(fatherGty[0]);
                    varUntransmitInfo.append('/');
                    varUntransmitInfo.append(fatherGty[1]);
                    varUntransmitInfo.append(',');
                    varUntransmitInfo.append(motherGty[0]);
                    varUntransmitInfo.append('/');
                    varUntransmitInfo.append(motherGty[1]);
                    varUntransmitInfo.append(',');
                    varUntransmitInfo.append(uchildGty[0]);
                    varUntransmitInfo.append('/');
                    varUntransmitInfo.append(uchildGty[1]);
                    varUntransmitInfo.append(";");
                }

                if (isEffective) {
                    varTransmitInfo.append(tVar.refStartPosition).append(":");
                    varTransmitInfo.append(fatherGty[0]);
                    varTransmitInfo.append('/');
                    varTransmitInfo.append(fatherGty[1]);
                    varTransmitInfo.append(',');
                    varTransmitInfo.append(motherGty[0]);
                    varTransmitInfo.append('/');
                    varTransmitInfo.append(motherGty[1]);
                    varTransmitInfo.append(',');
                    varTransmitInfo.append(childGty[0]);
                    varTransmitInfo.append('/');
                    varTransmitInfo.append(childGty[1]);
                    varTransmitInfo.append(";");
                    tmpVarListIndiv.add(tVar);
                }
                /*
                 * if (fatherGty[0] != fatherGty[1] && motherGty[0] !=
                 * motherGty[1] && childGty[0] != 0 && childGty[1] != 0) {
                 * System.out.println(lastGeneSymb + "\t" +
                 * mIndivChild.getLabelInChip() + "\t" + varInfo.toString()); }
                 */
            }

            if (transmitAltFromFaSudo && transmitAltFromMoSudo) {
                // controlDoubleHitGenes.add(lastGeneSymb);
                sudoControlGeneCounts++;
                if (allPsudoControl) {
                    hasADouble = true;
                    couts2[j + fixedColNum] = "Sudo:" + varUntransmitInfo.substring(0, varUntransmitInfo.length() - 1);
                }

                // System.out.println(mIndivFather.getLabelInChip() + "," +
                // mIndivMother.getLabelInChip() + "," +
                // mIndivChild.getLabelInChip() + ";" + varUntransmitInfo);
            }

            if (transmitAltFromFa && transmitAltFromMo) {
                couts1[j + fixedColNum] = String.valueOf(1);
                if (couts2[j + fixedColNum] != null) {
                    couts2[j + fixedColNum] = couts2[j + fixedColNum] + varTransmitInfo.substring(0, varTransmitInfo.length() - 1);
                } else {
                    couts2[j + fixedColNum] = varTransmitInfo.substring(0, varTransmitInfo.length() - 1);
                }

                if (mIndivChild.getAffectedStatus() == 2) {
                    caseGenes.add(lastGeneSymb);
                    caseGeneCounts++;
                }
                hasADouble = true;
                for (Variant tVar : tmpVarListIndiv) {
                    varLabel.append(tVar.refStartPosition).append(":").append(tVar.getRefAllele());
                    for (String altA : tVar.getAltAlleles()) {
                        varLabel.append(":");
                        varLabel.append(altA);
                    }
                    if (!effectiveVarSet.contains(varLabel.toString())) {
                        tmpVarListGene.add(tVar);
                        effectiveVarSet.add(varLabel.toString());
                    }
                    varLabel.delete(0, varLabel.length());
                }
            } else {
                couts1[j + fixedColNum] = "0";
                if (couts2[j + fixedColNum] == null) {
                    couts2[j + fixedColNum] = ".";
                }
            }

            if (countN >= 2) {
                hitIndivCount1[j]++;
            }
        }// end of scan at all individuals
        couts1[3] = String.valueOf(caseGeneCounts);
        couts2[3] = String.valueOf(caseGeneCounts);
        if (allPsudoControl) {
            couts1[4] = String.valueOf(sudoControlGeneCounts);
            couts2[4] = String.valueOf(sudoControlGeneCounts);
        }
        return hasADouble;
    }

    public void doubleHitGeneExploreVarPhasedGty(Chromosome chromosome, OpenLongObjectHashMap wahBit, int[] pedEncodeGytIDMap, int[] caseIDs, int[] controlIDs, List<String> pubmedMeshList,
            Map<String, String[]> geneNamesMap, Map<String, String> genePubMedID, boolean noNeedHomo, int pathogenicPredicIndex,
            List<String[]> hitDisCountsGenes, List<String[]> hitDisCounReads, Set<String> caseDoubleHitGenes, Set<String> controlDoubleHitGenes, FiltrationSummarySet doubleHitModelFilter, boolean needSearchPubMed)
            throws Exception {
        NCBIRetriever ncbiRetriever = new NCBIRetriever();

        if (pubmedMeshList == null || pubmedMeshList.isEmpty()) {
            needSearchPubMed = false;
        }

        Set<String> caseGenes = new HashSet<String>();
        Set<String> controlGenes = new HashSet<String>();

        String predicType = null;

        List<String> geneNames = new ArrayList<String>();

        Map<String, int[]> controlGeneCounts = new HashMap<String, int[]>();
        Map<String, int[]> caseGeneCounts = new HashMap<String, int[]>();

        boolean hasADouble = false;
        int account = 0;
        String lastGeneFeature = null;
        StringBuilder varInfo = new StringBuilder();
        Set<String> effectiveVarSet = new HashSet<String>();
        List<Variant> tmpVarListIndiv = new ArrayList<Variant>();
        List<Variant> tmpVarListGene = new ArrayList<Variant>();
        List<Variant> tmpVarListChrom = new ArrayList<Variant>();

        StringBuilder varLabel = new StringBuilder();
        int leftVarNum = 0;
        int controlNum = controlIDs == null ? 0 : controlIDs.length;
        int caseNum = caseIDs == null ? 0 : caseIDs.length;

        int effectiveIndivSize = caseNum + controlNum;

        if (chromosome.variantList == null || chromosome.variantList.isEmpty()) {
            return;
        }

        List<Variant> geneDisVars = null;
        List<Variant> geneNeuVars = null;

        boolean ignoreNonPathogenic = false;
        int[] caseGeneDoublehitSynoCounts = new int[1];
        int[] controlGeneDoublehitSynoCounts = new int[1];

        Map<String, List<Variant>> neuGeneVarMap = new HashMap<String, List<Variant>>();
        Map<String, List<Variant>> disGeneVarMap = new HashMap<String, List<Variant>>();

        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                if (ignoreNonPathogenic && pathogenicPredicIndex >= 0) {
                    predicType = var.featureValues[pathogenicPredicIndex];
                    if (predicType != null && predicType.equals("N")) {
                        continue;
                    }
                }
                //assume the synonymouse variants are neutural
                if (var.smallestFeatureID == 7) {
                    geneNeuVars = neuGeneVarMap.get(var.geneSymb);
                    if (geneNeuVars == null) {
                        geneNeuVars = new ArrayList<Variant>();
                        neuGeneVarMap.put(var.geneSymb, geneNeuVars);
                    }
                    geneNeuVars.add(var);
                } else {
                    geneDisVars = disGeneVarMap.get(var.geneSymb);
                    if (geneDisVars == null) {
                        geneDisVars = new ArrayList<Variant>();
                        disGeneVarMap.put(var.geneSymb, geneDisVars);
                    }
                    geneDisVars.add(var);
                }
            }
        }

        int fixedColNum = 5;

        long[][] counts1 = {{7, 12}, {18, 13}};

        boolean hasLess5 = false;
        double prob = 0;
        OpenIntIntHashMap fatherTransmitted = new OpenIntIntHashMap();
        OpenIntIntHashMap motherTransmitted = new OpenIntIntHashMap();
        int base = 0;
        int alleleNum = 0;
        boolean[] bits = new boolean[32];
        int startIndex;

        for (Map.Entry<String, List<Variant>> item : disGeneVarMap.entrySet()) {
            String geneSymbDis = item.getKey();

            geneDisVars = item.getValue();
            geneNeuVars = neuGeneVarMap.get(geneSymbDis);
            tmpVarListGene.clear();

            hasADouble = false;
            effectiveVarSet.clear();
            if (!geneDisVars.isEmpty()) {
                String[] couts1 = new String[effectiveIndivSize + fixedColNum];
                String[] couts2 = new String[effectiveIndivSize + fixedColNum];
                tmpVarListGene.clear();
                hasADouble = false;

                int[] childGty = null;
                int countH = 0;
                int countU = 0;
                for (int j = 0; j < controlNum; j++) {
                    int index = controlIDs[j];
                    varInfo.delete(0, varInfo.length());
                    tmpVarListIndiv.clear();
                    fatherTransmitted.clear();
                    motherTransmitted.clear();

                    for (Variant tVar : geneDisVars) {
                        alleleNum = tVar.getAltAlleles().length + 1;
                        base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
                        if (tVar.compressedGty) {
                            startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[index];
                            for (int i = 0; i < base; i++) {
                                bits[i] = wahBit.containsKey(startIndex);
                                startIndex += pedEncodeGytIDMap.length;
                            }
                            childGty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[index]);
                        } else {
                            childGty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[index], pedEncodeGytIDMap.length);
                        }

                        if (childGty == null) {
                            continue;
                        }
                        if (noNeedHomo && childGty[0] == childGty[1]) {
                            continue;
                        }

                        // assume the reference allele has not hitting
                        // effect.
                        if (childGty[0] != 0 || childGty[1] != 0) {
                            fatherTransmitted.put(childGty[0], 0);
                            motherTransmitted.put(childGty[1], 0);
                            varInfo.append(tVar.refStartPosition).append(":");
                            varInfo.append(childGty[0]);
                            varInfo.append('|');
                            varInfo.append(childGty[1]);
                            varInfo.append(";");
                            tmpVarListIndiv.add(tVar);
                            //to enable  if (fatherTransmitted.size() > 1 && motherTransmitted.size() > 1) work
                            if (childGty[0] != 0 && childGty[1] != 0) {
                                fatherTransmitted.put(0, 0);
                                motherTransmitted.put(0, 0);
                            }
                        }
                    }

                    if (fatherTransmitted.size() > 1 && motherTransmitted.size() > 1) {
                        couts1[j + fixedColNum] = String.valueOf(1);
                        couts2[j + fixedColNum] = varInfo.substring(0, varInfo.length() - 1);
                        countH++;
                        controlGenes.add(geneSymbDis);
                        hasADouble = true;
                        for (Variant tVar : tmpVarListIndiv) {
                            varLabel.append(tVar.refStartPosition).append(":").append(tVar.getRefAllele());
                            for (String altA : tVar.getAltAlleles()) {
                                varLabel.append(":");
                                varLabel.append(altA);
                            }
                            if (!effectiveVarSet.contains(varLabel.toString())) {
                                tmpVarListGene.add(tVar);
                                effectiveVarSet.add(varLabel.toString());
                            }
                            varLabel.delete(0, varLabel.length());
                        }

                        int[] geneCount = controlGeneCounts.get(geneSymbDis);
                        if (geneCount == null) {
                            controlGeneCounts.put(geneSymbDis, new int[]{1});
                        } else {
                            geneCount[0]++;
                        }
                    } else {
                        couts1[j + fixedColNum] = "0";
                        couts2[j + fixedColNum] = ".";
                        countU++;
                    }
                }// end of scan at all controls

                if (controlNum > 0) {
                    couts1[3] = String.valueOf(countH);
                    couts2[3] = String.valueOf(countH);
                    counts1[1][0] = countH;
                    counts1[1][1] = countU;
                }

                countH = 0;
                countU = 0;
                for (int j = 0; j < caseNum; j++) {
                    int index = caseIDs[j];

                    varInfo.delete(0, varInfo.length());
                    tmpVarListIndiv.clear();
                    fatherTransmitted.clear();
                    motherTransmitted.clear();
                    for (Variant tVar : geneDisVars) {
                        alleleNum = tVar.getAltAlleles().length + 1;
                        base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
                        if (tVar.compressedGty) {
                            startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[index];
                            for (int i = 0; i < base; i++) {
                                bits[i] = wahBit.containsKey(startIndex);
                                startIndex += pedEncodeGytIDMap.length;
                            }
                            childGty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[index]);
                        } else {
                            childGty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[index], pedEncodeGytIDMap.length);
                        }

                        if (childGty == null) {
                            continue;
                        }
                        if (noNeedHomo && childGty[0] == childGty[1]) {
                            continue;
                        }

                        // assume the reference allele has not hitting
                        // effect.
                        if (childGty[0] != 0 || childGty[1] != 0) {
                            fatherTransmitted.put(childGty[0], 0);
                            motherTransmitted.put(childGty[1], 0);
                            varInfo.append(tVar.refStartPosition).append(":");
                            varInfo.append(childGty[0]);
                            varInfo.append('|');
                            varInfo.append(childGty[1]);
                            varInfo.append(";");
                            tmpVarListIndiv.add(tVar);
                            //to enable  if (fatherTransmitted.size() > 1 && motherTransmitted.size() > 1) work
                            if (childGty[0] != 0 && childGty[1] != 0) {
                                fatherTransmitted.put(0, 0);
                                motherTransmitted.put(0, 0);
                            }
                        }
                    }

                    if (fatherTransmitted.size() > 1 && motherTransmitted.size() > 1) {
                        couts1[j + controlNum + fixedColNum] = String.valueOf(1);
                        couts2[j + controlNum + fixedColNum] = varInfo.substring(0, varInfo.length() - 1);
                        caseGenes.add(geneSymbDis);
                        hasADouble = true;
                        for (Variant tVar : tmpVarListIndiv) {
                            varLabel.append(tVar.refStartPosition).append(":").append(tVar.getRefAllele());
                            for (String altA : tVar.getAltAlleles()) {
                                varLabel.append(":");
                                varLabel.append(altA);
                            }
                            if (!effectiveVarSet.contains(varLabel.toString())) {
                                tmpVarListGene.add(tVar);
                                effectiveVarSet.add(varLabel.toString());
                            }
                            varLabel.delete(0, varLabel.length());
                        }
                        countH++;
                        int[] geneCount = caseGeneCounts.get(geneSymbDis);
                        if (geneCount == null) {
                            caseGeneCounts.put(geneSymbDis, new int[]{1});
                        } else {
                            geneCount[0]++;
                        }

                    } else {
                        couts1[j + controlNum + fixedColNum] = "0";
                        couts2[j + controlNum + fixedColNum] = ".";
                        countU++;
                    }
                }// end of scan at all cases

                if (caseNum > 0) {
                    couts1[4] = String.valueOf(countH);
                    couts2[4] = String.valueOf(countH);

                    counts1[0][0] = countH;
                    counts1[0][1] = countU;
                }
                if (controlNum > 0 && caseNum > 0) {
                    if (counts1[0][0] <= 5) {
                        hasLess5 = true;
                    } else if (counts1[0][1] <= 5) {
                        hasLess5 = true;
                    } else if (counts1[1][0] <= 5) {
                        hasLess5 = true;
                    } else if (counts1[1][1] <= 5) {
                        hasLess5 = true;
                    }

                    if (hasLess5) {
                        prob = ContingencyTable.fisherExact22(counts1, 2, 2, 2);
                    } else {
                        prob = ContingencyTable.pearsonChiSquared22(counts1);
                        prob = Probability.chiSquareComplemented(1, prob);
                    }
                    couts1[2] = Util.formatPValue(prob);
                    couts2[2] = couts1[2];
                }
                if (hasADouble) {
                    if (needSearchPubMed) {
                        account++;
                        geneNames.clear();

                        // unforturnately after I add the alais, there
                        // are too many false postive hits
                        // so I finally withraw this function
                        String[] names = geneNamesMap.get(geneSymbDis);
                        if (names != null && names.length > 0) {
                            geneNames.addAll(Arrays.asList(names));
                        }

                        geneNames.add(geneSymbDis);
                        LOG.info(account + ": Searching NCBI PubMed for " + pubmedMeshList.toString() + " and " + geneNames.toString());
                        while ((lastGeneFeature = ncbiRetriever.pubMedIDESearch(pubmedMeshList, geneNames, GlobalManager.pubMedFilter)) == null) {
                            // System.out.print("reconnecting...");
                        }
                        // System.out.println();
                        genePubMedID.put(geneSymbDis, lastGeneFeature);
                    }

                    couts1[0] = geneSymbDis;
                    couts1[1] = lastGeneFeature;

                    couts2[0] = geneSymbDis;
                    couts2[1] = lastGeneFeature;

                    hitDisCountsGenes.add(couts1);
                    hitDisCounReads.add(couts2);
                    tmpVarListChrom.addAll(tmpVarListGene);
                }
            }
        }
        caseDoubleHitGenes.addAll(caseGeneCounts.keySet());
        controlDoubleHitGenes.addAll(controlGeneCounts.keySet());
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarListChrom);
        tmpVarListChrom.clear();
        chromosome.buildVariantIndexMap();
        doubleHitModelFilter.increaseCount(0, chromosome.variantList.size());

    }

    public void doubleHitGeneExploreVarTriosSudoControl(Chromosome chromosome, OpenLongObjectHashMap wahBit, boolean isPhasedGty, List<Individual> sortedSubjectList, int[] pedEncodeGytIDMap, List<int[]> triosIDList, IntArrayList effectiveIndivIDs, List<String> pubmedMeshList,
            Map<String, String[]> geneNamesMap, Map<String, String> genePubMedID, boolean noNeedHomo, List<String[]> hitDisCountsGenes, List<String[]> hitDisCounReads, Set<String> caseDoubleHitGenes, Set<String> controlDoubleHitGenes,
            boolean allPsudoControl, FiltrationSummarySet doubleHitModelFilter, boolean needSearchPubMed) throws Exception {
        NCBIRetriever ncbiRetriever = new NCBIRetriever();

        if (pubmedMeshList == null || pubmedMeshList.isEmpty()) {
            needSearchPubMed = false;
        }
        int indivSize = sortedSubjectList.size();

        indivSize = triosIDList.size();
        List<String> geneNames = new ArrayList<String>();

        int effectiveIndivSize = effectiveIndivIDs.size();

        int[] hitIndivCount1 = new int[indivSize];
        int[] hitIndivCount2 = new int[indivSize];
        Arrays.fill(hitIndivCount1, 0);

        boolean hasADouble = false;
        int account = 0;
        String lastGeneFeature = null;

        Set<String> effectiveVarSet = new HashSet<String>();
        List<Variant> tmpVarListGene = new ArrayList<Variant>();
        List<Variant> tmpVarListChrom = new ArrayList<Variant>();

        int fixedColNum = 4;

        if (allPsudoControl) {
            fixedColNum = 5;
        }
        if (chromosome.variantList == null || chromosome.variantList.isEmpty()) {
            return;
        }

        List<Variant> geneDisVars = null;
        List<Variant> geneNeuVars = null;
        effectiveVarSet.clear();

        int[] caseGeneDoublehitSynoCounts = new int[1];
        int[] controlGeneDoublehitSynoCounts = new int[1];

        Map<String, List<Variant>> neuGeneVarMap = new HashMap<String, List<Variant>>();
        Map<String, List<Variant>> disGeneVarMap = new HashMap<String, List<Variant>>();
        int chromID = chromosome.getId();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                //assume the synonymouse variants are neutural
                if (var.smallestFeatureID == 7) {
                    geneNeuVars = neuGeneVarMap.get(var.geneSymb);
                    if (geneNeuVars == null) {
                        geneNeuVars = new ArrayList<Variant>();
                        neuGeneVarMap.put(var.geneSymb, geneNeuVars);
                    }
                    geneNeuVars.add(var);
                } else {
                    geneDisVars = disGeneVarMap.get(var.geneSymb);
                    if (geneDisVars == null) {
                        geneDisVars = new ArrayList<Variant>();
                        disGeneVarMap.put(var.geneSymb, geneDisVars);
                    }
                    geneDisVars.add(var);
                }
            }
        }

        long[][] counts1 = {{7, 12}, {18, 13}};

        boolean hasLess5 = false;
        double prob = 0;
        for (Map.Entry<String, List<Variant>> item : disGeneVarMap.entrySet()) {
            String geneSymbDis = item.getKey();

            geneDisVars = item.getValue();
            geneNeuVars = neuGeneVarMap.get(geneSymbDis);
            tmpVarListGene.clear();
            String[] couts1 = new String[effectiveIndivSize + fixedColNum];
            String[] couts2 = new String[effectiveIndivSize + fixedColNum];
            hasADouble = false;

            if (!geneDisVars.isEmpty()) {
                Arrays.fill(caseGeneDoublehitSynoCounts, 0);
                Arrays.fill(controlGeneDoublehitSynoCounts, 0);
                hasADouble = checkCompoundHeteroGeneSudo(effectiveIndivIDs, sortedSubjectList, pedEncodeGytIDMap, triosIDList, geneDisVars, chromID, isPhasedGty, noNeedHomo, geneSymbDis,
                        caseDoubleHitGenes, controlDoubleHitGenes, tmpVarListGene, wahBit, hitIndivCount1, effectiveVarSet,
                        couts1, couts2, allPsudoControl, fixedColNum);
                if (hasADouble) {
                    if (needSearchPubMed) {
                        account++;
                        geneNames.clear();

                        // unforturnately after I add the alais,
                        // there are too many false postive hits
                        // so I finally withraw this function
                        String[] names = geneNamesMap.get(geneSymbDis);
                        if (names != null && names.length > 0) {
                            geneNames.addAll(Arrays.asList(names));
                        }

                        geneNames.add(geneSymbDis);
                        LOG.info(account + ": Searching NCBI PubMed for " + pubmedMeshList.toString() + " and " + geneNames.toString());
                        while ((lastGeneFeature = ncbiRetriever.pubMedIDESearch(pubmedMeshList, geneNames, GlobalManager.pubMedFilter)) == null) {
                            // System.out.print("reconnecting...");
                        }
                        // System.out.println();
                        genePubMedID.put(geneSymbDis, lastGeneFeature);
                    }

                    couts1[0] = geneSymbDis;
                    couts1[1] = lastGeneFeature;

                    couts2[0] = geneSymbDis;
                    couts2[1] = lastGeneFeature;

                    tmpVarListChrom.addAll(tmpVarListGene);

                    //Do not do association analysis
                    /*
                     if (geneNeuVars != null && !geneNeuVars.isEmpty()) {
                     checkCompoundHeteroGeneSudoSimple(effectiveIndivIDs, pedEncodeGytIDMap, triosIDList, geneNeuVars, chromID, isPhasedGty, noNeedHomo, geneSymbDis,
                     caseGeneDoublehitSynoCounts, controlGeneDoublehitSynoCounts);
                     }
  
                     couts1[3] = String.valueOf(caseGeneDoublehitSynoCounts[0]);
                     couts2[3] = String.valueOf(caseGeneDoublehitSynoCounts[0]);
                     if (allPsudoControl) {
                     couts1[5] = String.valueOf(controlGeneDoublehitSynoCounts[0]);
                     couts2[5] = String.valueOf(controlGeneDoublehitSynoCounts[0]);
                     counts1[0][0] = Integer.parseInt(couts1[3]);
                     counts1[0][1] = Integer.parseInt(couts1[4]);
                     counts1[1][0] = Integer.parseInt(couts1[5]);
                     counts1[1][1] = Integer.parseInt(couts1[6]);

                     hasLess5 = false;
                     for (int t = 0; t < 2; t++) {
                     for (int k = 0; k < 2; k++) {
                     if (counts1[t][k] <= 5) {
                     hasLess5 = true;
                     break;
                     }
                     }
                     }
                     if (hasLess5) {
                     prob = ContingencyTable.fisherExact22(counts1, 2, 2, 2);
                     } else {
                     prob = ContingencyTable.pearsonChiSquared22(counts1);
                     prob = Probability.chiSquareComplemented(1, prob);
                     }

                     couts1[7] = String.valueOf(prob);
                     couts2[7] = String.valueOf(prob);
                     }
                     */
                    hitDisCountsGenes.add(couts1);
                    hitDisCounReads.add(couts2);
                }
                //do statistical test
            }

        }

        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarListChrom);
        tmpVarListChrom.clear();
        chromosome.buildVariantIndexMap();
        doubleHitModelFilter.increaseCount(0, chromosome.variantList.size());

    }

    public boolean checkCompoundHeteroGeneSudoSimple(IntArrayList effectiveIndivIDs, int[] pedEncodeGytIDMap, List<int[]> triosIDList, List<Variant> geneVars,
            OpenLongObjectHashMap wahBit, int chromID, boolean isPhasedGty, boolean noNeedHomo, String lastGeneSymb, int[] caseGeneCounts, int[] sudoControlGeneCounts) {
        int effectiveIndivSize = effectiveIndivIDs.size();
        boolean hasADouble = false;
        int base = 0;
        int alleleNum = 0;
        boolean isEffective;
        boolean[] bits = new boolean[32];
        int startIndex;
        for (int j = 0; j < effectiveIndivSize; j++) {
            int index = effectiveIndivIDs.getQuick(j);

            int countN = 0;
            boolean transmitAltFromFa = false;
            boolean transmitAltFromMo = false;
            boolean transmitAltFromFaSudo = false;
            boolean transmitAltFromMoSudo = false;

            for (Variant tVar : geneVars) {
                isEffective = false;
                int[] childGty = null;
                int[] fatherGty = null;
                int[] motherGty = null;
                alleleNum = tVar.getAltAlleles().length + 1;
                if (isPhasedGty) {
                    base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
                    if (tVar.compressedGty) {
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[0]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        childGty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[0]]);
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[1]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        fatherGty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[1]]);
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[2]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        motherGty = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[2]]);
                    } else {
                        childGty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[0]], pedEncodeGytIDMap.length);
                        fatherGty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[1]], pedEncodeGytIDMap.length);
                        motherGty = BinaryGtyProcessor.getPhasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[2]], pedEncodeGytIDMap.length);
                    }

                } else {
                    base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
                    if (tVar.compressedGty) {
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[0]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        childGty = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[0]]);
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[1]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        fatherGty = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[1]]);
                        startIndex = tVar.encodedGty[0] + pedEncodeGytIDMap[triosIDList.get(index)[2]];
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        motherGty = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[2]]);
                    } else {
                        childGty = BinaryGtyProcessor.getUnphasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[0]], pedEncodeGytIDMap.length);
                        fatherGty = BinaryGtyProcessor.getUnphasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[1]], pedEncodeGytIDMap.length);
                        motherGty = BinaryGtyProcessor.getUnphasedGtyAt(tVar.encodedGty, alleleNum, base, pedEncodeGytIDMap[triosIDList.get(index)[2]], pedEncodeGytIDMap.length);
                    }
                }
                if (childGty == null || fatherGty == null || motherGty == null) {
                    continue;
                }
                if (noNeedHomo && childGty[0] == childGty[1]) {
                    continue;
                }
                // assume the reference allele has not hitting effect.
                if (childGty[0] == 0 && childGty[1] != 0) {
                    // Warning: this will exclude the 0/1 0/1 genotypes of
                    // parents
                    if ((childGty[1] == fatherGty[0] || childGty[1] == fatherGty[1]) && childGty[1] != motherGty[0] && childGty[1] != motherGty[1]) {
                        if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                            continue;
                        }
                        countN++;
                        transmitAltFromFa = true;
                        isEffective = true;
                    } else if ((childGty[1] == motherGty[0] || childGty[1] == motherGty[1]) && childGty[1] != fatherGty[0] && childGty[1] != fatherGty[1]) {
                        if (motherGty[0] != 0 && motherGty[1] != 0) {
                            continue;
                        }
                        countN++;
                        transmitAltFromMo = true;
                        isEffective = true;
                    }
                } else if (childGty[0] != 0 && childGty[1] != 0) {
                    if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                        continue;
                    }
                    if (motherGty[0] != 0 && motherGty[1] != 0) {
                        continue;
                    }

                    if ((childGty[0] == fatherGty[0] || childGty[0] == fatherGty[1]) && (childGty[1] == motherGty[0] || childGty[1] == motherGty[1])) {
                        countN++;
                        countN++;
                        transmitAltFromFa = true;
                        transmitAltFromMo = true;
                        isEffective = true;
                    } else if ((childGty[1] == fatherGty[0] || childGty[1] == fatherGty[1]) && (childGty[0] == motherGty[0] || childGty[0] == motherGty[1])) {
                        countN++;
                        countN++;
                        transmitAltFromFa = true;
                        transmitAltFromMo = true;
                        isEffective = true;
                    }
                } else if (childGty[0] != 0) {
                    // Warning: this will exclude the 0/1 0/1 genotypes of
                    // parents
                    if ((childGty[0] == fatherGty[0] || childGty[0] == fatherGty[1]) && childGty[0] != motherGty[0] && childGty[0] != motherGty[1]) {
                        if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                            continue;
                        }
                        countN++;
                        transmitAltFromFa = true;
                        isEffective = true;
                    } else if ((childGty[0] == motherGty[0] || childGty[0] == motherGty[1]) && childGty[0] != fatherGty[0] && childGty[0] != fatherGty[1]) {
                        if (motherGty[0] != 0 && motherGty[1] != 0) {
                            continue;
                        }
                        countN++;
                        transmitAltFromMo = true;
                        isEffective = true;
                    }
                }

                int[] uchildGty = unTransmittedGty(fatherGty, motherGty, childGty);

                if (uchildGty != null) {
                    if (uchildGty[0] == 0 && uchildGty[1] != 0) {
                        // Warning: this will exclude the F:0/1 & M:0/1->O:0/1
                        // genotypes of parents
                        if ((uchildGty[1] == fatherGty[0] || uchildGty[1] == fatherGty[1]) && uchildGty[1] != motherGty[0] && uchildGty[1] != motherGty[1]) {
                            if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                                continue;
                            }
                            transmitAltFromFaSudo = true;
                        } else if ((uchildGty[1] == motherGty[0] || uchildGty[1] == motherGty[1]) && uchildGty[1] != fatherGty[0] && uchildGty[1] != fatherGty[1]) {
                            if (motherGty[0] != 0 && motherGty[1] != 0) {
                                continue;
                            }
                            transmitAltFromMoSudo = true;
                        }
                    } else if (uchildGty[0] != 0 && uchildGty[1] != 0) {
                        // ignore variants with more than 2 alleles
                        if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                            continue;
                        }
                        if (motherGty[0] != 0 && motherGty[1] != 0) {
                            continue;
                        }

                        // then offspreing is homozygous
                        if ((uchildGty[0] == fatherGty[0] || uchildGty[0] == fatherGty[1]) && (uchildGty[1] == motherGty[0] || uchildGty[1] == motherGty[1])) {
                            transmitAltFromFaSudo = true;
                            transmitAltFromMoSudo = true;
                            isEffective = true;
                        } else if ((uchildGty[1] == fatherGty[0] || uchildGty[1] == fatherGty[1]) && (uchildGty[0] == motherGty[0] || uchildGty[0] == motherGty[1])) {
                            transmitAltFromFaSudo = true;
                            transmitAltFromMoSudo = true;
                        }
                    } else if (uchildGty[0] != 0) {
                        // Warning: this will exclude the F:0/1 & M:0/1->O:0/1
                        // genotypes of parents
                        if ((uchildGty[0] == fatherGty[0] || uchildGty[0] == fatherGty[1]) && uchildGty[0] != motherGty[0] && uchildGty[0] != motherGty[1]) {
                            if (fatherGty[0] != 0 && fatherGty[1] != 0) {
                                continue;
                            }
                            transmitAltFromFaSudo = true;

                        } else if ((uchildGty[0] == motherGty[0] || uchildGty[0] == motherGty[1]) && uchildGty[0] != fatherGty[0] && uchildGty[0] != fatherGty[1]) {
                            if (motherGty[0] != 0 && motherGty[1] != 0) {
                                continue;
                            }
                            transmitAltFromMoSudo = true;
                        }
                    }
                }

                /*
                 * if (fatherGty[0] != fatherGty[1] && motherGty[0] !=
                 * motherGty[1] && childGty[0] != 0 && childGty[1] != 0) {
                 * System.out.println(lastGeneSymb + "\t" +
                 * mIndivChild.getLabelInChip() + "\t" + varInfo.toString()); }
                 */
            }

            if (transmitAltFromFa && transmitAltFromMo) {
                caseGeneCounts[0]++;
                hasADouble = true;
            }

            if (transmitAltFromFaSudo && transmitAltFromMoSudo) {
                sudoControlGeneCounts[0]++;
            }
        }// end of scan at all individuals
        return hasADouble;
    }

    public int[] unTransmittedGty(int[] fatherGty, int[] motherGty, int[] childGty) {
        int[] unTransGty = new int[2];
        if (childGty[0] == fatherGty[0] && childGty[1] == motherGty[0]) {
            unTransGty[0] = fatherGty[1];
            unTransGty[1] = motherGty[1];
        } else if (childGty[0] == fatherGty[0] && childGty[1] == motherGty[1]) {
            unTransGty[0] = fatherGty[1];
            unTransGty[1] = motherGty[0];
        } else if (childGty[0] == fatherGty[1] && childGty[1] == motherGty[0]) {
            unTransGty[0] = fatherGty[0];
            unTransGty[1] = motherGty[1];
        } else if (childGty[0] == fatherGty[1] && childGty[1] == motherGty[1]) {
            unTransGty[0] = fatherGty[0];
            unTransGty[1] = motherGty[0];
        } else if (childGty[0] == motherGty[0] && childGty[1] == fatherGty[0]) {
            unTransGty[0] = motherGty[1];
            unTransGty[1] = fatherGty[1];
        } else if (childGty[0] == motherGty[0] && childGty[1] == fatherGty[1]) {
            unTransGty[0] = motherGty[1];
            unTransGty[1] = fatherGty[0];
        } else if (childGty[0] == motherGty[1] && childGty[1] == fatherGty[0]) {
            unTransGty[0] = motherGty[0];
            unTransGty[1] = fatherGty[1];
        } else if (childGty[0] == motherGty[1] && childGty[1] == fatherGty[1]) {
            unTransGty[0] = motherGty[0];
            unTransGty[1] = fatherGty[0];
        } else {
            return null;
        }

        return unTransGty;
    }

    public void geneFeatureFilter(Chromosome chromsome, int[] variantsCounters, Set<Byte> featureInSet, FiltrationSummarySet geneModelFilter) throws Exception {
        int leftVariantNum = 0;
        List<Variant> temVariants = new ArrayList<Variant>();

        if (chromsome == null) {
            return;
        }
        temVariants.clear();
        temVariants.addAll(chromsome.variantList);
        chromsome.variantList.clear();

        for (Variant var : temVariants) {
            if (featureInSet.contains(var.smallestFeatureID)) {
                chromsome.variantList.add(var);
                leftVariantNum++;
            }
            variantsCounters[var.smallestFeatureID]++;
        }

        chromsome.buildVariantIndexMap();
        geneModelFilter.increaseCount(0, leftVariantNum);

    }

    public void matchTumorNontumorPair(String[] indivPairs, List<Individual> subjectIDList, List<int[]> pairedCancerSampleList, List<String> pairedCancerSampleLableList)
            throws Exception {
        int[] caseControSatus = new int[2];
        StringBuilder sb = new StringBuilder();
        for (String pair : indivPairs) {
            String[] invdivs = pair.split(":");
            int[] pariedIDs = new int[2];
            pariedIDs[0] = -9;
            pariedIDs[1] = -9;
            Arrays.fill(caseControSatus, 0);

            for (int t = 0; t < subjectIDList.size(); t++) {
                if (subjectIDList.get(t).getLabelInChip().equals(invdivs[0])) {
                    caseControSatus[0] = subjectIDList.get(t).getAffectedStatus();
                    pariedIDs[0] = t;
                } else if (subjectIDList.get(t).getLabelInChip().equals(invdivs[1])) {
                    caseControSatus[1] = subjectIDList.get(t).getAffectedStatus();
                    pariedIDs[1] = t;
                }
                if (pariedIDs[0] != -9 && pariedIDs[1] != -9) {
                    break;
                }

            }
            if (caseControSatus[0] == caseControSatus[1] && caseControSatus[0] == 1) {
                String infor = "Both " + invdivs[0] + " and " + invdivs[1] + " have unaffected status (1)!";
                throw new Exception(infor);
            } else if (caseControSatus[0] == caseControSatus[1] && caseControSatus[0] == 2) {
                String infor = "Both " + invdivs[0] + " and " + invdivs[1] + " have affected status (1)!";
                throw new Exception(infor);
            }
            if (caseControSatus[0] == 2 && caseControSatus[1] == 1) {
                int tmp = pariedIDs[0];
                pariedIDs[0] = pariedIDs[1];
                pariedIDs[1] = tmp;

            }

            if (pariedIDs[0] != -9 && pariedIDs[1] != -9) {
                sb.append("Matched Pair (Nontumor<->Tumor): ").append(subjectIDList.get(pariedIDs[0]).getLabelInChip()).append(" <-> ").append(subjectIDList.get(pariedIDs[1]).getLabelInChip());
                sb.append("\n");
                pairedCancerSampleList.add(pariedIDs);
                pairedCancerSampleLableList.add(subjectIDList.get(pariedIDs[0]).getLabelInChip() + " <-> " + subjectIDList.get(pariedIDs[1]).getLabelInChip());
            }
        }

        LOG.info(sb.toString());
    }

    public void somaticMutationFilterVar(Chromosome chromosome, OpenLongObjectHashMap wahBit, boolean isPhased, List<Individual> subjectList, int[] pedEncodeGytIDMap, int[] controlSetID, List<int[]> setSampleIDList,
            List<String> setSampleLabelList, FiltrationSummarySet somaticModelFilter, double somatP) throws Exception {
        int pairedSampleLen = setSampleIDList.size();
        double p = 1;
        double odds = 1;
        StringBuilder readEnrichTestInfo = new StringBuilder();

        boolean isNotSomaticVar = true;
        int geneFeatureNum = somaticModelFilter.getAvailableFeatureIndex();

        int controlNum = controlSetID.length;

        OpenIntIntHashMap mutatedHomoAllele = new OpenIntIntHashMap();
        OpenIntIntHashMap mutatedAllele = new OpenIntIntHashMap();

        OpenIntIntHashMap controlHomoAllele = new OpenIntIntHashMap();
        OpenIntIntHashMap controlAllele = new OpenIntIntHashMap();

        int[] gtys = null;
        int[] gtys0 = null;
        int[] gtys1 = null;

        int subID = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();
        int hitNum = 0;
        boolean hasLess5;
        int somatEvenNum = 0;
        double a00, a01, a10, a11;
        long[][] readCountsInt = new long[2][2];
        boolean[] bits = new boolean[32];
        int base = 0;
        int alleleNum = 0;
        String currChr = chromosome.getName();

        int chromID = chromosome.getId();
        tmpVarList.clear();
        int[] countsT = new int[2];
        int[] countsNT = new int[2];
        int id0, id1, id2;
        long startIndex;
        for (Variant var : chromosome.variantList) {
            mutatedHomoAllele.clear();
            mutatedAllele.clear();
            controlHomoAllele.clear();
            controlAllele.clear();
            /*
             if (var.refStartPosition == 142728364) {
             int sss = 0;
             }
             */
            alleleNum = var.getAltAlleles().length + 1;
            if (isPhased) {
                base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
            } else {
                base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
            }
            for (int j = controlNum - 1; j >= 0; j--) {
                subID = controlSetID[j];
                subID = pedEncodeGytIDMap[subID];
                if (subID < 0) {
                    continue;
                }
                if (var.compressedGty) {
                    startIndex = var.encodedGtyIndex[0] + subID;
                    for (int i = 0; i < base; i++) {
                        bits[i] = wahBit.containsKey(startIndex);
                        startIndex += pedEncodeGytIDMap.length;
                    }
                    if (isPhased) {
                        gtys = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                    } else {
                        gtys = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                    }
                } else if (isPhased) {
                    gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, pedEncodeGytIDMap.length);
                } else {
                    gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, pedEncodeGytIDMap.length);
                }

                if (gtys == null) {
                    continue;
                }

                if (gtys[0] == gtys[1]) {
                    controlHomoAllele.put(gtys[0], 0);
                    controlAllele.put(gtys[0], 0);
                } else {
                    controlAllele.put(gtys[0], 0);
                    controlAllele.put(gtys[1], 0);
                }
            }

            p = 1;
            odds = 1;
            readEnrichTestInfo.delete(0, readEnrichTestInfo.length());
            // pList.clear();
            // orList.clear();
            isNotSomaticVar = true;
            somatEvenNum = 0;
            for (int s = 0; s < pairedSampleLen; s++) {
                int[] pairID = setSampleIDList.get(s);
                id0 = pairID[0];
                id0 = pedEncodeGytIDMap[id0];
                if (id0 < 0) {
                    continue;
                }
                if (var.compressedGty) {
                    startIndex = var.encodedGtyIndex[0] + id0;
                    for (int i = 0; i < base; i++) {
                        bits[i] = wahBit.containsKey(startIndex);
                        startIndex += pedEncodeGytIDMap.length;
                    }
                    if (isPhased) {
                        gtys0 = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, id0);
                    } else {
                        gtys0 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, id0);
                    }
                } else if (isPhased) {
                    gtys0 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, id0, pedEncodeGytIDMap.length);
                } else {
                    gtys0 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, id0, pedEncodeGytIDMap.length);
                }

                if (gtys0 == null) {
                    continue;
                }
                id1 = pairID[1];
                id1 = pedEncodeGytIDMap[id1];
                if (id1 < 0) {
                    continue;
                }
                if (var.compressedGty) {
                    startIndex = var.encodedGtyIndex[0] + id1;
                    for (int i = 0; i < base; i++) {
                        bits[i] = wahBit.containsKey(startIndex);
                        startIndex += pedEncodeGytIDMap.length;
                    }
                    if (isPhased) {
                        gtys1 = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, id1);
                    } else {
                        gtys1 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, id1);
                    }
                } else if (isPhased) {
                    gtys1 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, id1, pedEncodeGytIDMap.length);
                } else {
                    gtys1 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, id1, pedEncodeGytIDMap.length);
                }

                if (gtys1 == null) {
                    continue;
                }

                // assumption: non-tumor tissue is homozygous and tumor
                // tissue is heterogygous
                if (gtys0[0] != gtys1[0] || gtys0[1] != gtys1[1]) {
                    //the KED file has not readInfor
                    if (var.readInfor != null) {
                        countsT[0] = var.readInfor[id1 + id1];
                        countsT[1] = var.readInfor[id1 + id1 + 1];

                        countsNT[0] = var.readInfor[id0 + id0];
                        countsNT[1] = var.readInfor[id0 + id0 + 1];

                        /*
                         * % control case % ___________ % risk | A | B |Rs1 %
                         * |_____|_____| % wild | C | D |Rs2 % |_____|_____|____
                         * % Cs1 Cs2 N
                         */
                        Arrays.fill(readCountsInt[0], 0);
                        Arrays.fill(readCountsInt[1], 0);
                        readCountsInt[0][0] = (long) (countsT[0]);
                        readCountsInt[0][1] = (long) (countsT[1]);
                        readCountsInt[1][0] = (long) (countsNT[0]);
                        readCountsInt[1][1] = (long) (countsNT[1]);
                        hasLess5 = false;
                        for (int t = 0; t < 2; t++) {
                            for (int j = 0; j < 2; j++) {
                                if (readCountsInt[t][j] < 0) {
                                    readCountsInt[t][j] = 0;
                                }
                                if (readCountsInt[t][j] < 5) {
                                    hasLess5 = true;
                                }
                            }
                        }

                        a00 = readCountsInt[0][0] == 0 ? 0.5 : readCountsInt[0][0];
                        a01 = readCountsInt[0][1] == 0 ? 0.5 : readCountsInt[0][1];
                        a10 = readCountsInt[1][0] == 0 ? 0.5 : readCountsInt[1][0];
                        a11 = readCountsInt[1][1] == 0 ? 0.5 : readCountsInt[1][1];
                        odds = (a01 * a10 / (a00 * a11));
                        if (hasLess5) {
                            p = ContingencyTable.fisherExact22(readCountsInt, 2, 2, 2);
                        } else {
                            p = ContingencyTable.pearsonChiSquared22(readCountsInt);
                            p = Probability.chiSquareComplemented(1, p);
                        }
                        if (p > somatP) {
                            continue;
                        }
                        // orList.add(odds);
                        readEnrichTestInfo.append(setSampleLabelList.get(s)).append(':').append(countsNT[0]).append('/').append(countsNT[1]).append(':').append(countsT[0])
                                .append('/').append(countsT[1]).append(':').append(Util.formatPValue(p)).append(':').append(odds).append('|');
                    } else {
                        readEnrichTestInfo.append(setSampleLabelList.get(s)).append('|');
                    }

                    // the tumor sample is a homozygous: possible
                    // reccessie model
                    if (gtys1[0] == gtys1[1]) {
                        // if
                        // (!controlSharedHomoAllele.contains(gtys1[0]))
                        {
                            // when the tumor sample is a homozygous it
                            // is ver likely the non-tumor is a
                            // heterozygous
                            // this senario is error-prone
                            isNotSomaticVar = false;
                            // System.out.println(gtys[pairID[0]] + " "
                            // + gtys[pairID[1]] + " " + currentLine);
                            somatEvenNum++;
                        }
                    } else // the tumor sample is a heterozygous: possible
                    // dominant model
                    {
                        if (gtys1[0] != gtys0[0] && gtys1[0] != gtys0[1]) // the
                        // charAt(0)
                        // must be
                        // the new
                        // mutant
                        {
                            // if (obsA > 0 && obsU > 0)
                            {
                                // all non-tumor sample have no such new
                                // mutant
                                // dominant model
                                // if
                                // (!controlSharedAllele.contains(gtys1[0]))
                                {
                                    isNotSomaticVar = false;
                                    // System.out.println(gtys[pairID[0]]
                                    // + " " + gtys[pairID[1]] + " " +
                                    // currentLine);
                                    somatEvenNum++;
                                }
                            }
                        } else if (gtys1[1] != gtys0[0] && gtys1[1] != gtys0[1]) {
                            // the charAt(2) must be the new mutant
                            // if (obsA > 0 && obsU > 0)
                            {
                                // all non-tumor sample have no such new
                                // mutant
                                // dominant model
                                // if
                                // (!controlSharedAllele.contains(gtys1[1]))
                                {
                                    isNotSomaticVar = false;
                                    // System.out.println(gtys[pairID[0]]
                                    // + " " + gtys[pairID[1]] + " " +
                                    // currentLine);
                                    somatEvenNum++;
                                }
                            }
                        }
                    }
                }
            }
            if (isNotSomaticVar) {
                continue;
            }
            var.setFeatureValue(geneFeatureNum, String.valueOf(somatEvenNum));
            if (readEnrichTestInfo.length() > 0) {
                var.setFeatureValue(geneFeatureNum + 1, readEnrichTestInfo.substring(0, readEnrichTestInfo.length() - 1));
            } else {
                var.setFeatureValue(geneFeatureNum + 1, null);
            }
            tmpVarList.add(var);
        }
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarList);
        hitNum += tmpVarList.size();
        tmpVarList.clear();

        chromosome.buildVariantIndexMap();
        somaticModelFilter.increaseCount(0, hitNum);
    }

    public void devnoMutationFilterVar(Chromosome chromosome, OpenLongObjectHashMap wahBit, boolean isPhased, List<Individual> subjectList, int[] pedEncodeGytIDMap, int[] controlSetID, List<int[]> setSampleIDList, boolean ingoreHomoInControl,
            FiltrationSummarySet denovoModelFilter) throws Exception {
        int controlNum = controlSetID.length;
        int hardFilteringNum = 0;

        OpenIntIntHashMap mutatedHomoAllele = new OpenIntIntHashMap();
        OpenIntIntHashMap mutatedAllele = new OpenIntIntHashMap();

        OpenIntIntHashMap controlHomoAllele = new OpenIntIntHashMap();
        OpenIntIntHashMap controlAllele = new OpenIntIntHashMap();
        int geneFeatureNum = denovoModelFilter.getAvailableFeatureIndex();
        int[] gtys = null;
        int[] gtys0 = null;
        int[] gtys1 = null;
        int[] gtys2 = null;
        int subID = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();
        int hitNum = 0;

        StringBuilder denoTestInfo = new StringBuilder();
        boolean isNotDevnoMutVar = false;
        String currChr;
        int devoEvenNum = 0;
        int id0, id1, id2 = 0;

        if (chromosome == null) {
            return;
        }

        currChr = chromosome.getName();

        tmpVarList.clear();
        int base = 0;
        int alleleNum = 0;
        boolean[] bits = new boolean[32];
        long startIndex;
        for (Variant var : chromosome.variantList) {
            mutatedHomoAllele.clear();
            mutatedAllele.clear();
            controlHomoAllele.clear();
            controlAllele.clear();

            alleleNum = var.getAltAlleles().length + 1;
            if (isPhased) {
                base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
            } else {
                base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
            }
            if (ingoreHomoInControl) {
                for (int j = controlNum - 1; j >= 0; j--) {
                    subID = controlSetID[j];

                    if (subID < 0) {
                        continue;
                    }
                    subID = pedEncodeGytIDMap[subID];
                    if (var.compressedGty) {
                        startIndex = var.encodedGtyIndex[0] + subID;
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }

                        if (isPhased) {
                            gtys = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                        } else {
                            gtys = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                        }
                    } else if (isPhased) {
                        gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, pedEncodeGytIDMap.length);
                    } else {
                        gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, pedEncodeGytIDMap.length);
                    }

                    if (gtys == null) {
                        continue;
                    }

                    if (gtys[0] == gtys[1]) {
                        controlHomoAllele.put(gtys[0], 0);
                        controlAllele.put(gtys[0], 0);
                    } else {
                        controlAllele.put(gtys[0], 0);
                        controlAllele.put(gtys[1], 0);
                    }
                }
            }
            isNotDevnoMutVar = true;
            devoEvenNum = 0;
            denoTestInfo.delete(0, denoTestInfo.length());

            for (int[] setIDs : setSampleIDList) {
                id0 = setIDs[0];
                id0 = pedEncodeGytIDMap[id0];
                if (id0 < 0) {
                    continue;
                }
                if (var.compressedGty) {
                    startIndex = var.encodedGtyIndex[0] + id0;
                    for (int i = 0; i < base; i++) {
                        bits[i] = wahBit.containsKey(startIndex);
                        startIndex += pedEncodeGytIDMap.length;
                    }

                    if (isPhased) {
                        gtys0 = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, id0);
                    } else {
                        gtys0 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, id0);
                    }
                } else if (isPhased) {
                    gtys0 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, id0, pedEncodeGytIDMap.length);
                } else {
                    gtys0 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, id0, pedEncodeGytIDMap.length);
                }

                if (gtys0 == null) {
                    continue;
                }
                id1 = -1;
                id2 = -1;
                if (setIDs[1] != -9) {
                    id1 = setIDs[1];
                    id1 = pedEncodeGytIDMap[id1];

                    if (id1 < 0) {
                        gtys1 = null;
                    } else if (var.compressedGty) {
                        startIndex = var.encodedGtyIndex[0] + id1;
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        if (isPhased) {
                            gtys1 = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, id1);
                        } else {
                            gtys1 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, id1);
                        }
                    } else if (isPhased) {
                        gtys1 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, id1, pedEncodeGytIDMap.length);
                    } else {
                        gtys1 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, id1, pedEncodeGytIDMap.length);
                    }

                } else {
                    gtys1 = null;
                }

                if (setIDs[2] != -9) {
                    id2 = setIDs[2];
                    id2 = pedEncodeGytIDMap[id2];
                    if (id2 < 0) {
                        gtys2 = null;
                    } else if (var.compressedGty) {
                        startIndex = var.encodedGtyIndex[0] + id2;
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        if (isPhased) {
                            gtys2 = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, id2);
                        } else {
                            gtys2 = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, id2);
                        }
                    } else if (isPhased) {
                        gtys2 = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, id2, pedEncodeGytIDMap.length);
                    } else {
                        gtys2 = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, id2, pedEncodeGytIDMap.length);
                    }
                } else {
                    gtys2 = null;
                }

                if (gtys2 == null && gtys1 == null) {
                    continue;
                }

                if (setIDs[1] != -9 && setIDs[2] != -9 && gtys1 != null && gtys2 != null) {
                    // assumption:
                    if (gtys0[0] == gtys0[1]) {
                        if (ingoreHomoInControl && controlHomoAllele.containsKey(gtys0[0])) {
                            continue;
                        }
                        // homoszygous
                        // only check the father of females
                        if (!currChr.equals("X") || subjectList.get(setIDs[0]).getGender() == 2) {
                            if (gtys0[0] != gtys1[0] && gtys0[0] != gtys1[1]) {
                                //ignore the mutations controls
                                if (controlAllele.containsKey(gtys0[0])) {
                                    continue;
                                }
                                isNotDevnoMutVar = false;
                                // System.out.println(gtys[pairID[0]] + " "
                                // + gtys[pairID[1]] + " " + currentLine);
                                devoEvenNum++;
                                denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                        .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&')
                                        .append(subjectList.get(setIDs[1]).getLabelInChip()).append(":").append(gtys1[0]).append('/').append(gtys1[1]).append(":")
                                        .append((int) var.readInfor[id1 + id1]).append(',').append((int) var.readInfor[id1 + id1 + 1]).append('&')
                                        .append(subjectList.get(setIDs[2]).getLabelInChip()).append(":").append(gtys2[0]).append('/').append(gtys2[1]).append(":")
                                        .append((int) var.readInfor[id2 + id2]).append(',').append((int) var.readInfor[id2 + id2 + 1]).append('|');
                                mutatedHomoAllele.containsKey(gtys0[0]);

                            }
                        }
                        // check mother's genotype
                        if (gtys0[0] != gtys2[0] && gtys0[0] != gtys2[1]) {
                            //ignore the mutations controls
                            if (controlAllele.containsKey(gtys0[0])) {
                                continue;
                            }
                            isNotDevnoMutVar = false;
                            // System.out.println(gtys[pairID[0]] + " " +
                            // gtys[pairID[1]] + " " + currentLine);
                            devoEvenNum++;
                            denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                    .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&')
                                    .append(subjectList.get(setIDs[1]).getLabelInChip()).append(":").append(gtys1[0]).append('/').append(gtys1[1]).append(":")
                                    .append((int) var.readInfor[id1 + id1]).append(',').append((int) var.readInfor[id1 + id1 + 1]).append('&')
                                    .append(subjectList.get(setIDs[2]).getLabelInChip()).append(":").append(gtys2[0]).append('/').append(gtys2[1]).append(":")
                                    .append((int) var.readInfor[id2 + id2]).append(',').append((int) var.readInfor[id2 + id2 + 1]).append('|');
                            mutatedHomoAllele.put(gtys0[0], 0);
                        }
                    } else {
                        // heterozygous
                        if (gtys0[0] != gtys1[0] && gtys0[0] != gtys1[1] && gtys0[0] != gtys2[0] && gtys0[0] != gtys2[1]) {
                            //ignore the mutations controls
                            if (controlAllele.containsKey(gtys0[0])) {
                                continue;
                            }
                            isNotDevnoMutVar = false;
                            devoEvenNum++;
                            denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                    .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&')
                                    .append(subjectList.get(setIDs[1]).getLabelInChip()).append(":").append(gtys1[0]).append('/').append(gtys1[1]).append(":")
                                    .append((int) var.readInfor[id1 + id1]).append(',').append((int) var.readInfor[id1 + id1 + 1]).append('&')
                                    .append(subjectList.get(setIDs[2]).getLabelInChip()).append(":").append(gtys2[0]).append('/').append(gtys2[1]).append(":")
                                    .append((int) var.readInfor[id2 + id2]).append(',').append((int) var.readInfor[id2 + id2 + 1]).append('|');
                            mutatedAllele.put(gtys0[0], 0);
                        }
                        if (gtys0[1] != gtys1[0] && gtys0[1] != gtys1[1] && gtys0[1] != gtys2[0] && gtys0[1] != gtys2[1]) {
                            //ignore the mutations controls
                            if (controlAllele.containsKey(gtys0[1])) {
                                continue;
                            }
                            isNotDevnoMutVar = false;
                            devoEvenNum++;
                            denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                    .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&')
                                    .append(subjectList.get(setIDs[1]).getLabelInChip()).append(":").append(gtys1[0]).append('/').append(gtys1[1]).append(":")
                                    .append((int) var.readInfor[id1 + id1]).append(',').append((int) var.readInfor[id1 + id1 + 1]).append('&')
                                    .append(subjectList.get(setIDs[2]).getLabelInChip()).append(":").append(gtys2[0]).append('/').append(gtys2[1]).append(":")
                                    .append((int) var.readInfor[id2 + id2]).append(',').append((int) var.readInfor[id2 + id2 + 1]).append('|');

                            mutatedAllele.put(gtys0[1], 0);
                        }
                    }
                } else if (setIDs[1] != -9 && gtys1 != null) {
                    // assumption:
                    if (gtys0[0] == gtys0[1]) {
                        if (ingoreHomoInControl && controlHomoAllele.containsKey(gtys0[0])) {
                            continue;
                        }
                        // homoszygous
                        if (!currChr.equals("X") || subjectList.get(setIDs[0]).getGender() == 2) {
                            if (gtys0[0] != gtys1[0] && gtys0[0] != gtys1[1]) {
                                //ignore the mutations controls
                                if (controlAllele.containsKey(gtys0[0])) {
                                    continue;
                                }
                                isNotDevnoMutVar = false;
                                // System.out.println(gtys[pairID[0]] + " "
                                // + gtys[pairID[1]] + " " + currentLine);
                                devoEvenNum++;
                                if (gtys2 != null) {
                                    denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                            .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&')
                                            .append(subjectList.get(setIDs[1]).getLabelInChip()).append(":").append(gtys1[0]).append('/').append(gtys1[1]).append(":")
                                            .append((int) var.readInfor[id1 + id1]).append(',').append((int) var.readInfor[id1 + id1 + 1]).append('&')
                                            .append(subjectList.get(setIDs[2]).getLabelInChip()).append(":").append(gtys2[0]).append('/').append(gtys2[1]).append(":")
                                            .append((int) var.readInfor[id2 + id2]).append(',').append((int) var.readInfor[id2 + id2 + 1]).append('|');

                                } else if (setIDs[2] != -9) {
                                    denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                            .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&')
                                            .append(subjectList.get(setIDs[1]).getLabelInChip()).append(":").append(gtys1[0]).append('/').append(gtys1[1]).append(":")
                                            .append((int) var.readInfor[id1 + id1]).append(',').append((int) var.readInfor[id1 + id1 + 1]).append('&')
                                            .append(subjectList.get(setIDs[2]).getLabelInChip()).append(":").append('.').append('/').append('.').append(":").append(".")
                                            .append('|');
                                } else {
                                    denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                            .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&')
                                            .append(subjectList.get(setIDs[1]).getLabelInChip()).append(":").append(gtys1[0]).append('/').append(gtys1[1]).append(":")
                                            .append((int) var.readInfor[id1 + id1]).append(',').append((int) var.readInfor[id1 + id1 + 1]).append('&').append("Unknown").append(":")
                                            .append('.').append('/').append('.').append(":").append(".").append('|');
                                }
                                mutatedHomoAllele.put(gtys0[0], 0);
                            }
                        }

                    }
                } else if (setIDs[2] != -9 && gtys2 != null) {
                    // assumption:
                    if (gtys0[0] == gtys0[1]) {
                        if (ingoreHomoInControl && controlHomoAllele.containsKey(gtys0[0])) {
                            continue;
                        }

                        // homoszygous
                        if (gtys0[0] != gtys2[0] && gtys0[0] != gtys2[1]) {
                            //ignore the mutations controls
                            if (controlAllele.containsKey(gtys0[0])) {
                                continue;
                            }
                            isNotDevnoMutVar = false;
                            // System.out.println(gtys[pairID[0]] + " " +
                            // gtys[pairID[1]] + " " + currentLine);
                            devoEvenNum++;
                            if (gtys1 != null) {
                                denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                        .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&')
                                        .append(subjectList.get(setIDs[1]).getLabelInChip()).append(":").append(gtys1[0]).append('/').append(gtys1[1]).append(":")
                                        .append((int) var.readInfor[id1 + id1]).append(',').append((int) var.readInfor[id1 + id1 + 1]).append('&')
                                        .append(subjectList.get(setIDs[2]).getLabelInChip()).append(":").append(gtys2[0]).append('/').append(gtys2[1]).append(":")
                                        .append((int) var.readInfor[id2 + id2]).append(',').append((int) var.readInfor[id2 + id2 + 1]).append('|');
                            } else if (setIDs[1] != -9) {
                                denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                        .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&')
                                        .append(subjectList.get(setIDs[1]).getLabelInChip()).append(":").append('.').append('/').append('.').append(":").append('.')
                                        .append('&').append(subjectList.get(setIDs[2]).getLabelInChip()).append(":").append(gtys2[0]).append('/').append(gtys2[1]).append(":")
                                        .append((int) var.readInfor[id2 + id2]).append(',').append((int) var.readInfor[id2 + id2 + 1]).append('|');
                            } else {
                                denoTestInfo.append(subjectList.get(setIDs[0]).getLabelInChip()).append(":").append(gtys0[0]).append('/').append(gtys0[1]).append(":")
                                        .append((int) var.readInfor[id0 + id0]).append(',').append((int) var.readInfor[id0 + id0 + 1]).append('&').append("Unknown").append(":")
                                        .append('.').append('/').append('.').append(":").append('.').append('&').append(subjectList.get(setIDs[2]).getLabelInChip())
                                        .append(":").append(gtys2[0]).append('/').append(gtys2[1]).append(":")
                                        .append((int) var.readInfor[id2 + id2]).append(',').append((int) var.readInfor[id2 + id2 + 1]).append('|');
                            }
                            mutatedHomoAllele.put(gtys0[0], 0);
                        }
                    }
                }

            }

            if (isNotDevnoMutVar) {
                hardFilteringNum++;
                continue;
            }

            if (denoTestInfo.length() > 0) {
                var.setFeatureValue(geneFeatureNum, denoTestInfo.substring(0, denoTestInfo.length() - 1));
            } else {
                var.setFeatureValue(geneFeatureNum, "");
            }
            tmpVarList.add(var);

        }
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarList);
        hitNum += tmpVarList.size();
        tmpVarList.clear();

        chromosome.buildVariantIndexMap();
        // StringBuilder message = new StringBuilder();

        if (hardFilteringNum > 0) {
            // message.append(hardFilteringNum).append(" variants are ignored by de-filtering\n");
        }

        denovoModelFilter.increaseCount(0, hitNum);
        //message.append(hitNum).append(" variant(s) are left after filtration by denovo mutation!");

        //  LOG.info(message.toString());
    }

    public void sumFilterCaseControlVar(Chromosome chromosome, int reqHetA, int reqHomA, int reqHetU, int reqHomU, int minOBSA, int minOBSU,
            FiltrationSummarySet minMissingQCFilter, FloatArrayList mafList) throws Exception {
        // AffectedRefHomGtyNum\tAffectedHetGtyNum\tAffectedAltHomGtyNum\tUnaffectedRefHomGtyNum\tUnaffectedHetGtyNum\tUnaffectedAltHomGtyNum
        int subID = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();
        int hitNum = 0;
        int hetA;
        int homA;
        int hetU;
        int homU;
        int obsA;
        int obsU;
        int ignoredLineNumMinOBSA = 0;
        int ignoredLineNumMinOBSU = 0;
        int ignoredLineNumMinHetHomU = 0;
        int ignoredLineNumMinHetHomA = 0;

        if (chromosome == null) {
            return;
        }
        boolean calculateMAF = false;
        if (mafList != null) {
            calculateMAF = true;
        }
        double maf = 0;
        List<Variant> varList = chromosome.variantList;
        tmpVarList.clear();
        for (Variant var : varList) {
            /*
             * if (var.getAltAlleles().length > 3) { continue; }
             */
            hetA = var.getAffectedHetGtyNum();
            homA = var.getAffectedAltHomGtyNum();
            obsA = var.getAffectedRefHomGtyNum() + hetA + homA;

            hetU = var.getUnaffectedHetGtyNum();
            homU = var.getUnaffectedAltHomGtyNum();
            obsU = var.getUnaffectedRefHomGtyNum() + hetU + homU;

            if (obsA < minOBSA) {
                ignoredLineNumMinOBSA++;
                continue;
            }

            if (obsU < minOBSU) {
                ignoredLineNumMinOBSU++;
                continue;
            }

            if (reqHetU > 0 && reqHomU > 0) {
                if (hetU < reqHetU && homU < reqHomU) {
                    ignoredLineNumMinHetHomU++;
                    continue;
                }
            } else if (reqHetU > 0) {
                if (hetU < reqHetU) {
                    ignoredLineNumMinHetHomU++;
                    continue;
                }
            } else if (reqHomU > 0) {
                if (homU < reqHomU) {
                    ignoredLineNumMinHetHomU++;
                    continue;
                }
            }

            if (reqHetA > 0 && reqHomA > 0) {
                if (hetA < reqHetA && homA < reqHomA) {
                    ignoredLineNumMinHetHomA++;
                    continue;
                }
            } else if (reqHetA > 0) {
                if (hetA < reqHetA) {
                    ignoredLineNumMinHetHomA++;
                    continue;
                }
            } else if (reqHomA > 0) {
                if (homA < reqHomA) {
                    ignoredLineNumMinHetHomA++;
                    continue;
                }
            }
            if (calculateMAF) {
                maf = (homA + 0.5 * hetA + homU + 0.5 * hetU);
                maf = maf / (obsA + obsU);
                if (maf > 0.5) {
                    maf = 1 - maf;
                }
                mafList.add((float) maf);
            }

            tmpVarList.add(var);
        }
        varList.clear();
        varList.addAll(tmpVarList);
        hitNum += tmpVarList.size();
        tmpVarList.clear();

        chromosome.buildVariantIndexMap();

        if (ignoredLineNumMinHetHomA > 0) {
            minMissingQCFilter.increaseCount(0, ignoredLineNumMinHetHomA);
        }

        if (ignoredLineNumMinHetHomU > 0) {
            minMissingQCFilter.increaseCount(1, ignoredLineNumMinHetHomU);
        }

        if (ignoredLineNumMinOBSA > 0) {
            minMissingQCFilter.increaseCount(2, ignoredLineNumMinOBSA);
        }

        if (ignoredLineNumMinOBSU > 0) {
            minMissingQCFilter.increaseCount(3, ignoredLineNumMinOBSU);
        }

        if (reqHetA > 0 || reqHomA > 0 || reqHetU > 0 || reqHomU > 0 || minOBSA > 0 || minOBSU > 0) {
            minMissingQCFilter.increaseCount(4, hitNum);
        }

    }

    public void inheritanceModelFilterVar(Chromosome chromosome, OpenLongObjectHashMap wahBit, boolean isPhased, List<Individual> subjectList, int[] pedEncodeGytIDMap, int[] caseSetID, int[] controlSetID, boolean[] genotypeFilters, String hardFilterModel,
            FiltrationSummarySet inheritanceModelFilter) throws Exception {
        // AffectedRefHomGtyNum\tAffectedHetGtyNum\tAffectedAltHomGtyNum\tUnaffectedRefHomGtyNum\tUnaffectedHetGtyNum\tUnaffectedAltHomGtyNum

        int caseNum = 0;
        if (caseSetID != null) {
            caseNum = caseSetID.length;
        }
        int controlNum = 0;
        if (controlSetID != null) {
            controlNum = controlSetID.length;
        }
        if (caseNum == 0 && controlNum == 0) {
            return;
        }
        int hardFilteringNum = 0;

        OpenIntIntHashMap caseSharedHomoAllele = new OpenIntIntHashMap();
        OpenIntIntHashMap caseSharedHeteAllele = new OpenIntIntHashMap();
        OpenIntIntHashMap controlSharedHomoAllele = new OpenIntIntHashMap();
        OpenIntIntHashMap caseSharedAllele = new OpenIntIntHashMap();
        OpenIntIntHashMap controlSharedAllele = new OpenIntIntHashMap();
        boolean modeMatched = true;

        int[] gtys = null;
        int subID = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();
        int hitNum = 0;

        if (chromosome == null) {
            return;
        }
        tmpVarList.clear();

        int base = 0;
        int alleleNum = 0;
        boolean[] bits = new boolean[32];
        long startIndex;
        int totalSubjectNum = pedEncodeGytIDMap.length;
        for (Variant var : chromosome.variantList) {
            caseSharedHomoAllele.clear();
            caseSharedHeteAllele.clear();
            controlSharedHomoAllele.clear();
            caseSharedAllele.clear();
            controlSharedAllele.clear();
            alleleNum = var.getAltAlleles().length + 1;
            if (isPhased) {
                base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
            } else {
                base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
            }
            modeMatched = true;
            // http://www.price.com.hk/product.php?p=166515


            /*
             * if (var.getAltAlleles().length > 3) { continue; }
             */
 /*
             * 
             * Code\t\tApplicable mode 1\tExclude variants which have
             * heterozygous genotypes in one or more affected family
             * members\tRecessive 2\tExclude variants which have the same
             * homozygous genotypes in both affected and unaffected family
             * members\tRecessive and full penetrance causal mutation(s) in
             * the sample 3\tExclude variants which have reference
             * homozygous genotype in one or more affected family
             * members\tRare dominant and no consanguineous mating or
             * compound-heterozygosity 4\tExclude variants which have the
             * same heterozygous genotypes in both affected and unaffected
             * family members\tDominant; full penetrance causal mutation(s)
             * in the sample 5\tExclude variants which have alternative
             * homozygous genotype in one or more affected family
             * members\tRare dominant and no consanguineous mating or
             * compound-heterozygosity 6\tExclude variants which have NO
             * shared alleles in affected family members
             */
 /*
             if (var.refStartPosition == 56102470) {
             int sss = 0;
             }
             */
            for (int j = 0; j < caseNum; j++) {
                subID = caseSetID[j];
                subID = pedEncodeGytIDMap[subID];
                if (subID < 0) {
                    continue;
                }
                if (var.compressedGty) {
                    startIndex = var.encodedGtyIndex[0] + subID;
                    for (int i = 0; i < base; i++) {
                        bits[i] = wahBit.containsKey(startIndex);
                        startIndex += pedEncodeGytIDMap.length;
                    }
                    if (isPhased) {
                        gtys = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                    } else {
                        gtys = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                    }
                } else if (isPhased) {
                    gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, pedEncodeGytIDMap.length);
                } else {
                    gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, pedEncodeGytIDMap.length);
                }

                if (gtys == null) {
                    continue;
                }

                // for rare monogenic disorders, cases are homozygous
                if (gtys[0] == gtys[1]) {
                    // 3\tExclude variants which have reference homozygous
                    // genotype in one or more affected family members\tRare
                    // dominant and no consanguineous mating or
                    // compound-heterozygosity
                    if (genotypeFilters[2] && gtys[0] == 0) {
                        modeMatched = false;
                        break;
                    }
                    // 5\tExclude variants which have alternative homozygous
                    // genotype in one or more affected family members\tRare
                    // dominant and no consanguineous mating or
                    // compound-heterozygosity
                    if (genotypeFilters[4] && gtys[0] != 0) {
                        modeMatched = false;
                        break;
                    }

                    if (caseSharedHomoAllele.isEmpty()) {
                        caseSharedHomoAllele.put(gtys[0], 0);
                    } else // 6\tExclude variants which have NO shared alleles
                    // in all affected family members
                    {
                        if (genotypeFilters[5]) {
                            if (!caseSharedHomoAllele.containsKey(gtys[0])) {
                                modeMatched = false;
                                break;
                            }
                        } else {
                            // for 9 Exclude variants at which the cases and
                            // controls have the same number of homozygous
                            // mutation types
                            caseSharedHomoAllele.put(gtys[0], 0);
                        }
                    }
                    caseSharedAllele.put(gtys[0], 0);
                } else {
                    // 1\tExclude variants which have heterozygous genotypes
                    // in one or more affected family members\tRecessive
                    if (genotypeFilters[0]) {
                        modeMatched = false;
                        break;
                    }

                    if (caseSharedHeteAllele.isEmpty()) {
                        caseSharedHeteAllele.put(gtys[0], 0);
                        caseSharedHeteAllele.put(gtys[1], 0);
                    } else // 6\tExclude variants which have NO shared alleles
                    // in affected family members
                    {
                        if (genotypeFilters[5]) {
                            if (!caseSharedHeteAllele.containsKey(gtys[0])) {
                                caseSharedHeteAllele.removeKey(gtys[0]);
                            }
                            if (!caseSharedHeteAllele.containsKey(gtys[1])) {
                                caseSharedHeteAllele.removeKey(gtys[1]);
                            }
                            if (caseSharedHeteAllele.isEmpty()) {
                                modeMatched = false;
                                break;
                            }
                        }
                    }
                    caseSharedAllele.put(gtys[0], 0);
                    caseSharedAllele.put(gtys[1], 0);
                }
            }

            // if there is shared alleles in controls
            if (modeMatched) {
                for (int j = controlNum - 1; j >= 0; j--) {
                    subID = controlSetID[j];
                    subID = pedEncodeGytIDMap[subID];
                    if (subID < 0) {
                        continue;
                    }
                    if (var.compressedGty) {
                        startIndex = var.encodedGtyIndex[0] + subID;
                        for (int i = 0; i < base; i++) {
                            bits[i] = wahBit.containsKey(startIndex);
                            startIndex += pedEncodeGytIDMap.length;
                        }
                        if (isPhased) {
                            gtys = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, subID);
                        } else {
                            gtys = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, subID);
                        }
                    } else if (isPhased) {
                        gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, subID, pedEncodeGytIDMap.length);
                    } else {
                        gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, subID, pedEncodeGytIDMap.length);
                    }

                    if (gtys == null) {
                        continue;
                    }

                    if (gtys[0] == gtys[1]) {

                        // 2\tExclude variants which have the same
                        // homozygous genotypes in both affected and
                        // unaffected family members\tRecessive and full
                        // penetrance causal mutation(s) in the sample
                        // here assumes no heterogeneity
                        if (genotypeFilters[1]) {
                            // only consider the non shared allele
                            if (caseSharedHomoAllele.containsKey(gtys[0])) {
                                modeMatched = false;
                                break;
                            }
                        }
                        controlSharedAllele.put(gtys[0], 0);
                    } else {
                        // 4\tExclude variants which have the same
                        // heterozygous genotypes in both affected and
                        // unaffected family members\tDominant; full
                        // penetrance causal mutation(s) in the sample
                        if (genotypeFilters[3]) {
                            /*
                             * if (caseSharedHeteAllele.isEmpty()) {
                             * modeMatched = false; break; } else
                             */
                            if (!caseSharedHeteAllele.isEmpty() && caseSharedHeteAllele.containsKey(gtys[0]) && caseSharedHeteAllele.containsKey(gtys[1])) {
                                modeMatched = false;
                                break;
                            }
                        }
                        controlSharedAllele.put(gtys[0], 0);
                        controlSharedAllele.put(gtys[1], 0);
                    }
                }
            }

            if (modeMatched) {
                tmpVarList.add(var);
            } else {
                hardFilteringNum++;
            }
        }
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarList);
        hitNum += tmpVarList.size();
        tmpVarList.clear();

        // genome.setVarNum(hitNum);
        chromosome.buildVariantIndexMap();
        if (hardFilteringNum > 0) {
            inheritanceModelFilter.increaseCount(0, hardFilteringNum);
        }
        inheritanceModelFilter.increaseCount(1, hitNum);
        /*
         StringBuilder message = new StringBuilder();
         if (hardFilteringNum > 0) {
         message.append(hardFilteringNum).append(" variants are ignored by genotype-based hard-filtering;\n");
         }
         message.append(hitNum).append(" variant(s) are left after filtration according to inheritance mode at genotypes, ").append(hardFilterModel).append("!");
         LOG.info(message.toString());
         */
    }

    public void superDupFilter(Chromosome chromosome, AnnotationSummarySet ass, ReferenceGenome refGenome) throws Exception {

        List<Variant> tmpList = new ArrayList<Variant>();

        if (chromosome == null) {
            return;
        }
        for (Variant var : chromosome.variantList) {
            List<RefDup> cnvList = refGenome.getDupFeature(chromosome.getName(), var, new RNABoundaryIndex(0));
            if (cnvList == null || cnvList.isEmpty()) {
                tmpList.add(var);
            }
        }
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpList);

        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + tmpList.size());
        ass.setAnnotNum(ass.getTotalNum() - ass.getLeftNum());
        tmpList.clear();

        chromosome.buildVariantIndexMap();
        //StringBuilder info = new StringBuilder();
        //chromosome.setVarNum(leftVariantNum);    
    }

    public void filterByGeneVarNum(Chromosome chromosome, AnnotationSummarySet ass, int maxGeneVar, List<String> featureLabels, boolean isPhasedGty) {
        //List<String> featureLabels = genome.getVariantFeatureLabels();
        int driverPredicIndex = -1;

        List<String[]> hitCounts1 = new ArrayList<String[]>();

        // AffectedRefHomGtyNum\tAffectedHetGtyNum\tAffectedAltHomGtyNum\tUnaffectedRefHomGtyNum\tUnaffectedHetGtyNum\tUnaffectedAltHomGtyNum
        for (int i = 0; i < featureLabels.size(); i++) {
            if (driverPredicIndex >= 0) {
                break;
            }
        }

        IntArrayList caseIDs = new IntArrayList();
        IntArrayList controlIDs = new IntArrayList();
        IntArrayList unknownDiseaseIDs = new IntArrayList();

        // cluster counts according to phentoypes
        int caseNum = caseIDs.size();
        int controlNum = controlIDs.size();
        int unknownNum = unknownDiseaseIDs.size();

        Map<String, int[]> geneVarSum = new HashMap<String, int[]>();
        List<String> heads1 = new ArrayList<String>();
        heads1.add("Gene");
        heads1.add("Length");
        List<String> heads2 = new ArrayList<String>();
        heads2.add("Disease");
        heads2.add(".");

        hitCounts1.add(heads1.toArray(new String[0]));
        hitCounts1.add(heads2.toArray(new String[0]));

        //boolean isPhasedGty = genome.isIsPhasedGty();
        String lastGeneSymb = null;
        List<Variant> geneVars = new ArrayList<Variant>();
        List<Variant> tmpChromVars = new ArrayList<Variant>();
        int letfNum = 0;

        if (chromosome == null) {
            return;
        }
        lastGeneSymb = null;
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
                if (lastGeneSymb == null) {
                    geneVars.add(var);
                    lastGeneSymb = var.geneSymb;
                } else if (var.geneSymb.equals(lastGeneSymb)) {
                    geneVars.add(var);
                } else {
                    if (geneVars.size() < maxGeneVar) {
                        tmpChromVars.addAll(geneVars);
                    }
                    geneVars.clear();
                    geneVars.add(var);
                    lastGeneSymb = var.geneSymb;
                }
            } else {
                tmpChromVars.add(var);
            }
        }

        // the last gene
        if (!geneVars.isEmpty()) {
            if (geneVars.size() < maxGeneVar) {
                tmpChromVars.addAll(geneVars);
            }
            geneVars.clear();
        }
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + tmpChromVars.size());

        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpChromVars);
        letfNum += tmpChromVars.size();
        tmpChromVars.clear();

        chromosome.buildVariantIndexMap();
        //genome.buildVariantIndexMapOnChromosomes();

    }

    public void keepGenesInSet(Chromosome chromosome, AnnotationSummarySet ass, Set<String> geneSet) {
        //int varNum = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();

        if (chromosome == null) {
            return;
        }
        String geneFeatureAnnot, geneSymb;
        StringBuilder sb = new StringBuilder();
        int index;
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb != null) {
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
                            if (geneSet.contains(geneSymb)) {
                                tmpVarList.add(var);
                                break;
                            }
                        }
                    }
                    sb.delete(0, sb.length());
                }
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVarList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVarList.size());
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarList);
        tmpVarList.clear();
        //varNum += chromosome.variantList.size();

        chromosome.buildVariantIndexMap();
        //genome.buildVariantIndexMapOnChromosomes();
    }

    public void keepGenesOutSet(Chromosome chromosome, AnnotationSummarySet ass, Set<String> geneSet) {
//        int varNum = 0;
        List<Variant> tmpVarList = new ArrayList<Variant>();

        if (chromosome == null) {
            return;
        }
        String geneFeatureAnnot, geneSymb;
        StringBuilder sb = new StringBuilder();
        int index;
        boolean isNotContain = false;
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb == null) {
                tmpVarList.add(var);
                continue;
            }
            if (!geneSet.contains(var.geneSymb)) {
                tmpVarList.add(var);

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
                isNotContain = true;
                if (sb.length() > 0) {
                    String[] items = Util.tokenize(sb.substring(0, sb.length() - 1), ';');
                    for (String item : items) {
                        index = item.indexOf(':');
                        if (index > 0) {
                            geneSymb = item.substring(0, index);
                            if (geneSet.contains(geneSymb)) {
                                isNotContain = false;
                                break;
                            }
                        }
                    }
                    sb.delete(0, sb.length());
                }
                if (isNotContain) {
                    tmpVarList.add(var);
                }
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVarList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVarList.size());

        chromosome.variantList.clear();

        chromosome.variantList.addAll(tmpVarList);

        tmpVarList.clear();
//        varNum += chromosome.variantList.size();

        chromosome.buildVariantIndexMap();
//        genome.buildVariantIndexMapOnChromosomes();

    }

    public void filterBy4ModelPValue(Chromosome chromosome, AnnotationSummarySet ass, double[] pValueThresholds, Genome genome) {
        int[] pValueIndexes = new int[4];
        List<String> featureLabels = genome.getVariantFeatureLabels();
        Arrays.fill(pValueIndexes, -1);

        for (int i = 0; i < featureLabels.size(); i++) {
            if (featureLabels.get(i).equals("PValueAllelic")) {
                pValueIndexes[0] = i;
            } else if (featureLabels.get(i).equals("PValueDom")) {
                pValueIndexes[1] = i;
            } else if (featureLabels.get(i).equals("PValueRec")) {
                pValueIndexes[2] = i;
            } else if (featureLabels.get(i).equals("PValueGeno")) {
                pValueIndexes[3] = i;
            }
        }

        List<Variant> tmpVariantList = new ArrayList<Variant>();
        boolean atLeastOnePass = false;
        int leftVarNum = genome.getVarNum();

        if (chromosome == null) {
            return;
        }

        int len1 = chromosome.variantList.size();
        for (int j = 0; j < len1; j++) {
            Variant var = chromosome.variantList.get(j);
            atLeastOnePass = false;
            for (int t = 0; t < 4; t++) {
                double p = Double.parseDouble(var.getFeatureValues()[pValueIndexes[t]]);
                if (p <= pValueThresholds[t]) {
                    atLeastOnePass = true;
                    break;
                }
            }
            if (atLeastOnePass) {
                tmpVariantList.add(var);
            } else {
                leftVarNum--;
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVariantList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVariantList.size());

        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVariantList);
        tmpVariantList.clear();

        chromosome.buildVariantIndexMap();
        //genome.buildVariantIndexMapOnChromosomes();
//        genome.setVarNum(leftVarNum);
//        String info = leftVarNum + " variant(s) with p valus <= the above p-value cutoffs are left!";
//        LOG.info(info);
    }

    public void hardFilterInByRegions(Chromosome chromosome, AnnotationSummarySet ass, Map<String, int[]> regions) {

        boolean checkRegion = false;
        boolean isInRegion = false;
        boolean isWholeChrom = false;

        int varNum = 0;
        StringBuilder strRegion = new StringBuilder();

        if (chromosome == null) {
            return;
        }
        int[] physicalRegions = regions.get(chromosome.getName());
        if (physicalRegions == null) {
            chromosome.variantList.clear();
            return;
        }

        int regionNum = physicalRegions.length / 2;

        int j = 0;
        isWholeChrom = false;
        for (j = 0; j < regionNum; j++) {
            if (physicalRegions[j * 2] != -9 || physicalRegions[j * 2 + 1] != -9) {
                checkRegion = true;
                if (physicalRegions[j * 2 + 1] == -9) {
                    physicalRegions[1] = Integer.MAX_VALUE;
                }
                strRegion.append(" chr");
                strRegion.append(chromosome.getName());
                strRegion.append("[");
                strRegion.append(physicalRegions[j * 2]);
                strRegion.append(",");
                strRegion.append(physicalRegions[j * 2 + 1]);
                strRegion.append("]bp ");
            } else if (physicalRegions[j * 2] == -9 && physicalRegions[j * 2 + 1] == -9) {
                isWholeChrom = true;
                physicalRegions[1] = Integer.MAX_VALUE;
                strRegion.append("chr").append(chromosome.getName()).append(" ");
            }
        }

        if (isWholeChrom) {
            varNum += chromosome.variantList.size();
            return;
        }

        List<Variant> tmpVarList = new ArrayList();

        int vi = 0;
        for (Variant var1 : chromosome.variantList) {
            // System.out.println(currChr); 
            //filter physical region
            if (checkRegion) {
                isInRegion = false;
                for (int k = 0; k < regionNum; k++) {
                    if (var1.refStartPosition >= physicalRegions[k * 2] && var1.refStartPosition <= physicalRegions[k * 2 + 1]) {
                        isInRegion = true;
                        break;
                    }
                }

                if (isInRegion) {
                    tmpVarList.add(var1);
                    varNum++;
                }

            }
            vi++;
        }
        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVarList);

        ass.setAnnotNum(ass.getAnnotNum() + tmpVarList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVarList.size());

//        genome.setVarNum(varNum);
//        genome.buildVariantIndexMapOnChromosomes();
        chromosome.buildVariantIndexMap();
        // String info = varNum + " variant(s) are left in the regions " + strRegion.toString();
        // LOG.info(info);

    }

    public void hardFilterOutByRegions(Chromosome chromosome, AnnotationSummarySet ass, Map<String, int[]> regions) {

        boolean checkRegion = false;
        boolean isInRegion = false;
        boolean isWholeChrom = false;

        int varNum = 0;
        StringBuilder strRegion = new StringBuilder();

        if (chromosome == null) {
            return;
        }
        int[] physicalRegions = regions.get(chromosome.getName());
        if (physicalRegions == null) {
            varNum += chromosome.variantList.size();
            return;
        }

        int regionNum = physicalRegions.length / 2;

        int j = 0;
        isWholeChrom = false;
        for (j = 0; j < regionNum; j++) {
            if (physicalRegions[j * 2] != -9 || physicalRegions[j * 2 + 1] != -9) {
                checkRegion = true;
                if (physicalRegions[j * 2 + 1] == -9) {
                    physicalRegions[1] = Integer.MAX_VALUE;
                }
                strRegion.append(" [");
                strRegion.append(physicalRegions[j * 2]);
                strRegion.append(",");
                strRegion.append(physicalRegions[j * 2 + 1]);
                strRegion.append("]bp ");
            } else if (physicalRegions[j * 2] == -9 && physicalRegions[j * 2 + 1] == -9) {
                isWholeChrom = true;
                physicalRegions[1] = Integer.MAX_VALUE;
                strRegion.append("chr").append(chromosome.getName()).append(" ");
            }
        }

        if (isWholeChrom) {
            chromosome.variantList.clear();
            return;
        }

        List<Variant> tmpVarList = chromosome.variantList;
        chromosome.variantList.clear();
        int vi = 0;
        for (Variant var1 : tmpVarList) {
            // System.out.println(currChr); 
            //filter physical region
            if (checkRegion) {
                isInRegion = false;
                for (int k = 0; k < regionNum; k++) {
                    if (var1.refStartPosition >= physicalRegions[k * 2] && var1.refStartPosition <= physicalRegions[k * 2 + 1]) {
                        isInRegion = true;
                        break;
                    }
                }

                if (!isInRegion) {
                    tmpVarList.add(var1);
                    varNum++;
                }

            }
            vi++;
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVarList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVarList.size());

//        genome.setVarNum(varNum);
//        genome.buildVariantIndexMapOnChromosomes();
        chromosome.buildVariantIndexMap();
        String info = varNum + " variant(s) are left beyond the regions " + strRegion.toString();
        LOG.info(info);
    }

    public void filterByAlleleFreqExcModel(Chromosome chromosome, AnnotationSummarySet ass, float minAlleleFreq, FloatArrayList mafList) {

//        int leftVarNum = genome.getVarNum();
        int hitNum = 0;
        List<Variant> tmpVariantList = new ArrayList<Variant>();

        //else maxAlleleFreq >0
        if (chromosome == null) {
            return;
        }
        boolean needMAF = false;
        if (mafList != null) {
            needMAF = true;
        }
        for (Variant var : chromosome.variantList) {
            if (var.altAF == -1) {
                //not exist in any database
                //keep it anyhow
                tmpVariantList.add(var);
            } else if (Float.isNaN(var.altAF)) {
                //exist in a database but have no frequence informaion
                //not exist in any database
                tmpVariantList.add(var);
                /*
                 if (minAlleleFreq == 0f) {
                 //                    leftVarNum--;
                 hitNum++;
                 } else {
                 tmpVariantList.add(var);
                 }*/
            } else //exist in a database and have freq infor
             if ((var.altAF >= minAlleleFreq)) {
//                    leftVarNum--;
                    hitNum++;
                } else {
                    tmpVariantList.add(var);
                }
            if (needMAF) {
                if (var.altAF == -1 || Float.isNaN(var.altAF)) {
                    mafList.add(0);
                } else if (var.altAF <= 0.5) {
                    mafList.add(var.altAF);
                } else {
                    mafList.add(1 - var.altAF);
                }
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVariantList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVariantList.size());

        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVariantList);
        tmpVariantList.clear();

//        genome.setVarNum(leftVarNum);
//        genome.buildVariantIndexMapOnChromosomes();
        chromosome.buildVariantIndexMap();
//        String info = leftVarNum + " variant(s) with minor allele frequency [" + 0 + ", " + minAlleleFreq + ") in the reference datasets above are left!";
//        LOG.info(info);
    }

    public void filterByAlleleFreqIncModel(Chromosome chromosome, AnnotationSummarySet ass, float minAlleleFreq, float maxAlleleFreq, FloatArrayList mafList) {
//        Chromosome[] chroms = genome.getChromosomes();
//        int leftVarNum = genome.getVarNum();
        int hitNum = 0;
        List<Variant> tmpVariantList = new ArrayList<Variant>();

        //else maxAlleleFreq >0
        if (chromosome == null) {
            return;
        }
        boolean needMAF = false;
        if (mafList != null) {
            needMAF = true;
        }
        int len1 = chromosome.variantList.size();
        for (int j = 0; j < len1; j++) {
            Variant var = chromosome.variantList.get(j);
            if (var.altAF == -1) {
                //remove is anyhow
//                leftVarNum--;
                hitNum++;
            } else if (Float.isNaN(var.altAF)) {
                // exist in a database but have no frequence informaion
                //remove is anyhow
//                leftVarNum--;
                hitNum++;
            } else //exist in a database and have freq infor
             if ((var.altAF < minAlleleFreq || var.altAF > maxAlleleFreq)) {
//                    leftVarNum--;
                    hitNum++;
                } else {
                    tmpVariantList.add(var);
                }
            if (needMAF) {
                if (var.altAF == -1 || Float.isNaN(var.altAF)) {
                    mafList.add(0);
                } else if (var.altAF <= 0.5) {
                    mafList.add(var.altAF);
                } else {
                    mafList.add(1 - var.altAF);
                }
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + tmpVariantList.size());
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - tmpVariantList.size());

        chromosome.variantList.clear();
        chromosome.variantList.addAll(tmpVariantList);
        tmpVariantList.clear();

//        genome.buildVariantIndexMapOnChromosomes();
//        genome.setVarNum(leftVarNum);
        chromosome.buildVariantIndexMap();
//        String info = leftVarNum + " variant(s) with minor allele frequency [" + minAlleleFreq + ", " + maxAlleleFreq + "] in the reference datasets above are left!";
//        LOG.info(info);
    }
}
