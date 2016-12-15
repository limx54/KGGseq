/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

<<<<<<< HEAD
import cern.colt.list.DoubleArrayList;
=======
>>>>>>> origin/master
import cern.colt.map.OpenLongObjectHashMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
<<<<<<< HEAD
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
=======
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
>>>>>>> origin/master
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.GZIPOutputStream;
import org.apache.log4j.Logger;
import org.cobi.kggseq.Constants;
import static org.cobi.kggseq.Constants.STAND_CHROM_NAMES;
import org.cobi.kggseq.GlobalManager;
import static org.cobi.kggseq.GlobalManager.PLUGIN_PATH;
import org.cobi.kggseq.entity.Individual;
import org.cobi.kggseq.entity.Variant;
<<<<<<< HEAD
import org.cobi.util.plot.PValuePainter;
import org.cobi.util.stat.MultipleTestingMethod;
=======
>>>>>>> origin/master
import org.cobi.util.text.LocalExcelFile;
import org.cobi.util.thread.Task;

/**
 *
 * @author mxli
 */
public class RVTest {

    private static final Logger LOG = Logger.getLogger(RVTest.class);
    String rvtestFolder;
    String inputFileFolder;
    String outputFileFolder;
    String rvtestOutputNameTmp;
    HashMap<String, String[]> hmpRGroup;
    HashMap<String, String[]> hmpRSingle;
    int N = 4;
    int single = -1;
    int group = -1;
    String strCov;
    String strPhe;
    HashMap<String, String[]> hmpM2ST;//Single variant tests.
    HashMap<String, String[]> hmpM2BT;//Burden tests.
    HashMap<String, String[]> hmpM2TT;//Variable threshold models.
    HashMap<String, String[]> hmpM2KT;//Kernel models.
    HashMap<String, String[]> hmpM2MT;//Meta-analysis models.
    String strCommand;
//
//    enum SetTest {
//        cmc, kbac, price, skat
//    }

    public RVTest(String inputFileFolder, String outputFileFolder) throws Exception {
        this.rvtestFolder = PLUGIN_PATH + "/rvtests-master/executable/rvtest";

        File outFile = new File(outputFileFolder);
        File f = new File(outFile.getCanonicalPath());
        if (!f.getParentFile().exists()) {
            f.getParentFile().mkdir();
        }
        if (inputFileFolder == null) {
            this.inputFileFolder = outFile.getCanonicalPath() + ".flt.simple.vcf.gz";
        } else {
            File inFile = new File(inputFileFolder);
            this.inputFileFolder = inFile.getCanonicalPath();
        }

        this.rvtestOutputNameTmp = outFile.getCanonicalPath() + ".rvtestTMP" + (int) (Math.random() * 10000);
        this.outputFileFolder = outFile.getCanonicalPath();

        //Single variant tests
        hmpM2ST = new HashMap<String, String[]>();
        String strTemp[] = new String[]{".SingleScore.assoc", "PVALUE", "12"};
        hmpM2ST.put("score", strTemp);
        strTemp = new String[]{".SingleWald.assoc", "Pvalue", "8"};
        hmpM2ST.put("wald", strTemp);
        strTemp = new String[]{".FisherExact.assoc", "PvalueTwoSide", "11"};
        hmpM2ST.put("exact", strTemp);
        strTemp = new String[]{".DominantFisherExact.assoc", "PvalueTwoSide", "11"};
        hmpM2ST.put("dominantExact", strTemp);
        strTemp = new String[]{".FamLRT.assoc", "Pvalue", "8"};
        hmpM2ST.put("famLRT", strTemp);
        strTemp = new String[]{".FamScore.assoc", "Pvalue", "8"};
        hmpM2ST.put("famScore", strTemp);
        strTemp = new String[]{".FamGrammarGamma.assoc", "Pvalue", "8"};
        hmpM2ST.put("famGrammarGamma", strTemp);
        strTemp = new String[]{".SingleFirth.assoc", "Pvalue", "8"};
        hmpM2ST.put("firth", strTemp);

        //Burden tests
        hmpM2BT = new HashMap<String, String[]>();
        strTemp = new String[]{".CMC.assoc", "Pvalue", "6"};
        hmpM2BT.put("cmc", strTemp);
        strTemp = new String[]{".Zeggini.assoc", "Pvalue", "5"};
        hmpM2BT.put("zeggini", strTemp);
        strTemp = new String[]{".MadsonBrowning.assoc", "Pvalue", "10"};
        hmpM2BT.put("mb", strTemp);
        strTemp = new String[]{".Fp.assoc", "Pvalue", "5"};
        hmpM2BT.put("fp", strTemp);
        strTemp = new String[]{".CMCFisherExact.assoc", "PvalueTwoSide", "9"};
        hmpM2BT.put("exactCMC", strTemp);
        strTemp = new String[]{".CMCWald.assoc", "Pvalue", "8"};
        hmpM2BT.put("cmcWald", strTemp);
        strTemp = new String[]{".RareCover.assoc", "PermPvalue", "11"};
        hmpM2BT.put("rarecover", strTemp);
        strTemp = new String[]{".CMAT.assoc", "PermPvalue", "10"};//to be checked!
        hmpM2BT.put("cmat", strTemp);
        strTemp = new String[]{".FamCMC.assoc", "Pvalue", "10"};
        hmpM2BT.put("famcmc", strTemp);
        strTemp = new String[]{".FamZeggini.assoc", "Pvalue", "10"};
        hmpM2BT.put("famzeggini", strTemp);

        //variable threshold models
        hmpM2TT = new HashMap<String, String[]>();
        strTemp = new String[]{".VariableThresholdPrice.assoc", "PermPvalue", "13"};
        hmpM2TT.put("price", strTemp);
        strTemp = new String[]{".AnalyticVT.assoc", "Pvalue", "12"};
        hmpM2TT.put("analytic", strTemp);
        strTemp = new String[]{".FamAnalyticVT.assoc", "Pvalue", "12"};
        hmpM2TT.put("famAnalytic", strTemp);

        //kernel models
        hmpM2KT = new HashMap<String, String[]>();
        strTemp = new String[]{".Skat.assoc", "Pvalue", "6"};
        hmpM2KT.put("skat", strTemp);
        strTemp = new String[]{".SkatO.assoc", "Pvalue", "7"};
        hmpM2KT.put("skato", strTemp);
        strTemp = new String[]{".Kbac.assoc", "Pvalue", "5"};
        hmpM2KT.put("kbac", strTemp);
        strTemp = new String[]{".FamSkat.assoc", "Pvalue", "6"};
        hmpM2KT.put("famSkat", strTemp);

        //Meta-analysis models//Not be used currently. 
        hmpM2MT = new HashMap<String, String[]>();
        strTemp = new String[]{".Skat.assoc", "Pvalue", "6"};
        hmpM2MT.put("skat", strTemp);
        strTemp = new String[]{".SkatO.assoc", "Pvalue", "7"};
        hmpM2MT.put("skato", strTemp);
        strTemp = new String[]{".Kbac.assoc", "Pvalue", "5"};
        hmpM2MT.put("kbac", strTemp);
        strTemp = new String[]{".FamSkat.assoc", "Pvalue", "6"};
        hmpM2MT.put("famSkat", strTemp);

    }

    public void setCommand(String strCommand) {
        this.strCommand = strCommand;
    }

    public String getCmd(String strMarker) {
        String strCmdSingle = getCmdItem(hmpM2ST, "--single");
        String strCmdBurden = getCmdItem(hmpM2BT, "--burden");
        String strCmdVT = getCmdItem(hmpM2TT, "--vt");
        String strCmdKernel = getCmdItem(hmpM2KT, "--kernel");
        String strCmdGroup = strCmdBurden + strCmdVT + strCmdKernel;
        if (!strCmdSingle.isEmpty() && !strCmdGroup.isEmpty()) {
            System.out.println("rvTest cannot support both single variant and region based tests! The region based tests are carried out!");
            LOG.info("rvTest cannot support both single variant and region based tests! The region based tests are carried out!");
            return strCmdGroup;
        }
        if (strCmdSingle.isEmpty() && strCmdGroup.isEmpty()) {
            if (strMarker.contains("variant")) {
                this.setCommand("score,wald");
                return "--single score,wald";
            } else {
                this.setCommand("cmc,analytic,skato[nPerm=1000:alpha=0.001:beta1=1:beta2=20]");
                return "--burden cmc --vt analytic --kernel skato[nPerm=1000:alpha=0.001:beta1=1:beta2=20]";
            }
        }
        return strCmdSingle + strCmdBurden + strCmdVT + strCmdKernel;
    }

    public void setPheno(String phenoName) {
        strPhe = "--pheno-name " + phenoName;
    }

    public String getCmdItem(HashMap<String, String[]> hmpTemp, String strPara) {
        String strCmd = "";
        String[] strItem = strCommand.split(",");
        for (String str : hmpTemp.keySet()) {
            for (int i = 0; i < strItem.length; i++) {
                if (strItem[i].replaceAll("\\[.*\\]", "").equals(str)) {
                    strCmd += strItem[i] + ",";
                    break;
                }
            }
        }
        if (strCmd != "") {
            strCmd = strPara + " " + strCmd.substring(0, strCmd.length() - 1) + " ";
        }
        return strCmd;
    }

    public void setCov(String[] covItem) {
        if (covItem != null) {
            // strCov = " --covar " + pedFile + " --covar-name ";
            strCov = " --covar-name ";
            for (int i = 0; i < covItem.length; i++) {
                strCov += covItem[i] + ",";
            }
            strCov = strCov.substring(0, strCov.length() - 1);
        }
    }

    public void runTabix() {
        String bgizpFolder = PLUGIN_PATH + "/rvtests-master/third/tabix/";
        String cmd = bgizpFolder + "tabix -f " + inputFileFolder;
        try {

            //LocalFileFunc.gunzipFile(inputFileFolder, inputFileFolder.substring(0, inputFileFolder.length() - 3));
            String line;

            Process pr = Runtime.getRuntime().exec(cmd);
            try {
                BufferedReader inputError = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
                while (((line = inputError.readLine()) != null)) {
                    //  System.out.println(line);
                }
                int exitVal = pr.waitFor();
                pr.destroy();
                inputError.close();
                if (exitVal != 0) {
                    LOG.info("Failed to run the command:" + cmd);
                }
            } catch (Exception ex) {
                LOG.error(ex + "\nFailed to run the command:" + cmd);
            }
        } catch (IOException ex) {
            LOG.error(ex + "\nFailed to run the command:" + cmd);
        } catch (Exception ex) {
            LOG.error(ex + "\nFailed to run the command:" + cmd);
        }

    }

    public void collect(HashMap<String, String[]> hmpRGroup, String key, HashMap<String, String[]> hmpTemp, int order, boolean boolKeep, int sg, String strMarker) throws IOException {
        for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
            String strIn = rvtestOutputNameTmp + (strMarker.equals("geneset") ? ".set" : ".chr" + STAND_CHROM_NAMES[i]) + hmpTemp.get(key)[0];
            File file = new File(strIn);
<<<<<<< HEAD
            //  System.out.println(file.getCanonicalPath());
=======
>>>>>>> origin/master
            if (!file.exists()) {
                continue;
            }
            buildMap(hmpRGroup, boolKeep, file, order, hmpTemp.get(key)[1], sg, strMarker);
        }
    }

<<<<<<< HEAD
    public File drwaQQPlot(HashMap<String, String[]> hmpTemp, ArrayList<String> altTitle, String multipleTestingMethodLabel, double pValueThreshold, String outpath) throws Exception {
        PValuePainter pvPainter = new PValuePainter(450, 450);
        File plotFile2 = new File(outpath);

        String multiCorrMethodName = null;
        if (multipleTestingMethodLabel == null || multipleTestingMethodLabel.equals("no")) {
            multiCorrMethodName = "";
        } else if (multipleTestingMethodLabel.equals("benfdr")) {
            multiCorrMethodName = " by Benjamini-Hochberg FDR";
        } else if (multipleTestingMethodLabel.equals("bonhol")) {
            multiCorrMethodName = " by Bonferroni-Holm procedure";
        } else if (multipleTestingMethodLabel.equals("bonf")) {
            multiCorrMethodName = " by Standard Bonferroni";
        }

        StringBuilder message = new StringBuilder();
        if (multiCorrMethodName.length() > 0) {
            message.append("Significance level of p value cutoffs for RVTESTS p-values for the overal error rate ").append(pValueThreshold).append(multiCorrMethodName + ":\n");
        } else {
            message.append("Significance level of a fixed p value cutoff for RVTESTS p-values ").append(pValueThreshold).append(" : ");
        }
        DoubleArrayList[] genePVList = new DoubleArrayList[altTitle.size()];
        String tmpStr;
        for (int i = 0; i < altTitle.size(); i++) {
            double threshold = -1;
            genePVList[i] = new DoubleArrayList();

            for (Map.Entry<String, String[]> pMap : hmpTemp.entrySet()) {
                tmpStr = pMap.getValue()[i];
                if (tmpStr.equals("NA")) {
                    continue;
                }
                genePVList[i].add(Double.parseDouble(tmpStr));
                //System.out.println(Double.parseDouble(pMap.getValue()[i]));
            }
            genePVList[i].quickSort();
            message.append(" ");
            if (multipleTestingMethodLabel == null || multipleTestingMethodLabel.equals("no")) {
                threshold = pValueThreshold;
            } else if (multipleTestingMethodLabel.equals("benfdr")) {
                threshold = MultipleTestingMethod.benjaminiHochbergFDR(pValueThreshold, genePVList[i]);
            } else if (multipleTestingMethodLabel.equals("bonhol")) {
                threshold = MultipleTestingMethod.bonferroniHolmFWE(pValueThreshold, genePVList[i]);
            } else if (multipleTestingMethodLabel.equals("bonf")) {
                threshold = pValueThreshold / genePVList[i].size();
            }
            message.append(altTitle.get(i));
            message.append(": ");
            message.append(threshold).append("; ");
        }

        LOG.info(message.substring(0, message.length() - 1));
        pvPainter.drawMultipleQQPlot(Arrays.asList(genePVList), altTitle, null, plotFile2.getCanonicalPath(), 1E-10);
        String info = "The QQ plot saved in " + plotFile2.getCanonicalPath();
        LOG.info(info);
        return plotFile2;
    }

=======
>>>>>>> origin/master
    public void outputResult(HashMap<String, String[]> hmpTemp, ArrayList<String> altTitle, String[] affix, boolean excel, int sg) throws IOException, Exception {
        File fleOutput = new File(this.outputFileFolder + affix[0] + (excel ? "xlsx" : "txt"));
        if (!fleOutput.getParentFile().exists()) {
            fleOutput.getParentFile().mkdir();
        }
        if (!excel) {
            BufferedWriter bw = new BufferedWriter(new FileWriter(fleOutput));
            String strTemp = "";
            for (int i = 0; i < altTitle.size(); i++) {
                strTemp += altTitle.get(i) + "\t";
            }
            strTemp = strTemp.substring(0, strTemp.length() - 1);
            strTemp += "\n";
//            bw.write("Gene\tCMC\tKbac\tVariableThresholdPrice\tSkat\n");
            bw.write(strTemp);
            for (Map.Entry<String, String[]> entry : hmpTemp.entrySet()) {
                strTemp = entry.getKey();
                for (int i = 0; i < sg; i++) {
                    strTemp += "\t" + entry.getValue()[i];
                }
                strTemp += "\n";
                bw.write(strTemp);
            }
            bw.close();
        } else {
            List<String[]> arrays = new ArrayList<String[]>();
//            String[] titles = new String[]{"Gene", "CMC", "Kbac", "VariableThresholdPrice", "Skat"};
            String[] titles = new String[altTitle.size()];
            titles = altTitle.toArray(titles);

            for (Map.Entry<String, String[]> entry : hmpTemp.entrySet()) {
                String strTemp = entry.getKey();
                String[] cells = new String[sg + 1];
                cells[0] = strTemp;
                for (int i = 0; i < sg; i++) {
                    cells[i + 1] = entry.getValue()[i];
                }
                arrays.add(cells);
            }
            LocalExcelFile.writeArray2XLSXFile(fleOutput.getCanonicalPath(), titles, arrays);
        }
        String infor = null;
        if (affix[1].equals("Gene")) {
            infor = "The association analysis results of gene-base association are saved in " + fleOutput.getCanonicalPath();
        } else if (affix[1].equals("Variant")) {
            infor = "The association analysis results of variant-base association are saved in " + fleOutput.getCanonicalPath();
        } else if (affix[1].equals("GeneSet")) {
            infor = "The association analysis results of geneset-base association are saved in " + fleOutput.getCanonicalPath();
        }

        LOG.info(infor);
    }

<<<<<<< HEAD
    public ArrayList<String> collectResultGene(boolean keepGrp, boolean keepTmp, boolean excel) throws Exception {
=======
    public void collectResultGene(boolean boolKeep, boolean excel) throws Exception {
>>>>>>> origin/master
        String[] strItem = strCommand.split(",");
        ArrayList<String> altSingleTitle = new ArrayList<String>();
        altSingleTitle.add("Variant");
        ArrayList<String> altGroupTitle = new ArrayList<String>();
        altGroupTitle.add("Gene");

        for (int i = 0; i < strItem.length; i++) {
            if (hmpM2ST.keySet().contains(strItem[i].replaceAll("\\[.*\\]", ""))) {
                altSingleTitle.add(strItem[i].replaceAll("\\[.*\\]", ""));
            } else if (hmpM2BT.keySet().contains(strItem[i].replaceAll("\\[.*\\]", ""))) {
                altGroupTitle.add(strItem[i].replaceAll("\\[.*\\]", ""));
            } else if (hmpM2TT.keySet().contains(strItem[i].replaceAll("\\[.*\\]", ""))) {
                altGroupTitle.add(strItem[i].replaceAll("\\[.*\\]", ""));
            } else if (hmpM2KT.keySet().contains(strItem[i].replaceAll("\\[.*\\]", ""))) {
                altGroupTitle.add(strItem[i].replaceAll("\\[.*\\]", ""));
            }
        }
        single = altSingleTitle.size() - 1;
        group = altGroupTitle.size() - 1;

        if (group > 0) {
            int intCount = 0;
            hmpRGroup = new HashMap<String, String[]>();
            for (int i = 1; i < altGroupTitle.size(); i++) {
                if (hmpM2BT.keySet().contains(altGroupTitle.get(i))) {
<<<<<<< HEAD
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2BT, intCount, keepTmp, group, "gene");
                    intCount++;
                } else if (hmpM2TT.keySet().contains(altGroupTitle.get(i))) {
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2TT, intCount, keepTmp, group, "gene");
                    intCount++;
                } else if (hmpM2KT.keySet().contains(altGroupTitle.get(i))) {
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2KT, intCount, keepTmp, group, "gene");
=======
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2BT, intCount, boolKeep, group, "gene");
                    intCount++;
                } else if (hmpM2TT.keySet().contains(altGroupTitle.get(i))) {
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2TT, intCount, boolKeep, group, "gene");
                    intCount++;
                } else if (hmpM2KT.keySet().contains(altGroupTitle.get(i))) {
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2KT, intCount, boolKeep, group, "gene");
>>>>>>> origin/master
                    intCount++;
                }
            }
            outputResult(hmpRGroup, altGroupTitle, new String[]{".rvtest.gene.", "Gene"}, excel, group);
<<<<<<< HEAD

//            hmpRGroup.clear();
=======
            hmpRGroup.clear();
>>>>>>> origin/master
        }

        if (group <= 0 && single > 0) {
            int intCount = 0;
            hmpRSingle = new HashMap<String, String[]>();
            for (int i = 1; i < altSingleTitle.size(); i++) {
                if (hmpM2ST.keySet().contains(altSingleTitle.get(i))) {
<<<<<<< HEAD
                    collect(hmpRSingle, altSingleTitle.get(i), hmpM2ST, intCount, keepTmp, single, "var");
=======
                    collect(hmpRSingle, altSingleTitle.get(i), hmpM2ST, intCount, boolKeep, single, "var");
>>>>>>> origin/master
                    intCount++;
                }
            }
            outputResult(hmpRSingle, altSingleTitle, new String[]{".rvtest.variant.", "Variant"}, excel, single);
        }

        //remove the temporary products. 
        for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
            File fleRM = new File(rvtestOutputNameTmp + ".chr" + STAND_CHROM_NAMES[i] + ".log");
            if (fleRM.exists()) {
                fleRM.delete();
            }
        }
<<<<<<< HEAD
        if (!keepGrp) {
=======
        if (!boolKeep) {
>>>>>>> origin/master
            for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
                File fleRM = new File(outputFileFolder + ".chr" + STAND_CHROM_NAMES[i] + ".gene.rvtest.grp.gz");
                if (fleRM.exists()) {
                    fleRM.delete();
                }
            }
        }
<<<<<<< HEAD
        if (group <= 0 && single > 0) {
            return altSingleTitle;
        } else {
            return altGroupTitle;
        }
    }

    public ArrayList<String> collectResultGeneset(boolean keepGrp, boolean keepTmp, boolean excel) throws Exception {
=======

    }

    public void collectResultGeneset(boolean boolKeep, boolean excel) throws Exception {
>>>>>>> origin/master

        String[] strItem = strCommand.split(",");
        ArrayList<String> altGroupTitle = new ArrayList<String>();
        altGroupTitle.add("GeneSet");

        for (int i = 0; i < strItem.length; i++) {
            if (hmpM2BT.keySet().contains(strItem[i].replaceAll("\\[.*\\]", ""))) {
                altGroupTitle.add(strItem[i].replaceAll("\\[.*\\]", ""));
            } else if (hmpM2TT.keySet().contains(strItem[i].replaceAll("\\[.*\\]", ""))) {
                altGroupTitle.add(strItem[i].replaceAll("\\[.*\\]", ""));
            } else if (hmpM2KT.keySet().contains(strItem[i].replaceAll("\\[.*\\]", ""))) {
                altGroupTitle.add(strItem[i].replaceAll("\\[.*\\]", ""));
            }
        }

        group = altGroupTitle.size() - 1;
        if (group > 0) {
            int intCount = 0;
            hmpRGroup = new HashMap<String, String[]>();
            for (int i = 1; i < altGroupTitle.size(); i++) {
                if (hmpM2BT.keySet().contains(altGroupTitle.get(i))) {
<<<<<<< HEAD
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2BT, intCount, keepTmp, group, "geneset");
                    intCount++;
                } else if (hmpM2TT.keySet().contains(altGroupTitle.get(i))) {
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2TT, intCount, keepTmp, group, "geneset");
                    intCount++;
                } else if (hmpM2KT.keySet().contains(altGroupTitle.get(i))) {
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2KT, intCount, keepTmp, group, "geneset");
=======
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2BT, intCount, boolKeep, group, "geneset");
                    intCount++;
                } else if (hmpM2TT.keySet().contains(altGroupTitle.get(i))) {
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2TT, intCount, boolKeep, group, "geneset");
                    intCount++;
                } else if (hmpM2KT.keySet().contains(altGroupTitle.get(i))) {
                    collect(hmpRGroup, altGroupTitle.get(i), hmpM2KT, intCount, boolKeep, group, "geneset");
>>>>>>> origin/master
                    intCount++;
                }
            }
            outputResult(hmpRGroup, altGroupTitle, new String[]{".rvtest.geneset.", "GeneSet"}, excel, group);
<<<<<<< HEAD

        }

=======
            hmpRGroup.clear();
        }

//        hmpRGroup = new HashMap<String, double[]>();
//        String strIn = rvtestOutputNameTmp + ".set.CMC.assoc";
//        File file = new File(strIn);
//        if (file.exists()) {
//            buildMap(hmpRGroup, boolKeep, file, SetTest.cmc.ordinal(), 6,group,"group");            
//        } 
//
//        strIn = rvtestOutputNameTmp + ".set.Kbac.assoc";
//        file = new File(strIn);
//        if (file.exists()) {
//            buildMap(hmpRGroup, boolKeep, file, SetTest.kbac.ordinal(), 5,group,"group");
//        }
//        
//        strIn = rvtestOutputNameTmp + ".set.VariableThresholdPrice.assoc";
//        file = new File(strIn);
//        if (file.exists()) {
//            buildMap(hmpRGroup, boolKeep, file, SetTest.price.ordinal(), 13,group,"group");
//        }
//        strIn = rvtestOutputNameTmp + ".set.Skat.assoc";
//        file = new File(strIn);
//        if (file.exists()) {
//            buildMap(hmpRGroup, boolKeep, file, SetTest.skat.ordinal(), 6,group,"group");
//        }
        //output the result. 
//        File fleOutput = new File(this.outputFileFolder + ".rvtest.geneset." + (excel ? "xlsx" : "txt"));
//        if (!fleOutput.getParentFile().exists()) {
//            fleOutput.getParentFile().mkdir();
//        }
//        if (!excel) {
//            BufferedWriter bw = new BufferedWriter(new FileWriter(fleOutput));
//            bw.write("Set\tCMC\tKbac\tVariableThresholdPrice\tSkat\n");
//            for (Map.Entry<String, double[]> entry : hmpRGroup.entrySet()) {
//                String strTemp = entry.getKey();
//                for (int i = 0; i < N; i++) {
//                    strTemp += "\t" + (entry.getValue()[i] == -1 ? "." : entry.getValue()[i]);
//                }
//                strTemp += "\n";
//                bw.write(strTemp);
//            }
//            bw.close();
//        } else {
//            List<String[]> arrays = new ArrayList<String[]>();
//            String[] titles = new String[]{"Set", "CMC", "Kbac", "VariableThresholdPrice", "Skat"};
//
//            for (Map.Entry<String, double[]> entry : hmpRGroup.entrySet()) {
//                String strTemp = entry.getKey();
//                String[] cells = new String[5];
//                cells[0] = strTemp;
//                for (int i = 0; i < N; i++) {
//                    cells[i + 1] = (entry.getValue()[i] == -1 ? "." : String.valueOf(entry.getValue()[i]));
//                }
//                arrays.add(cells);
//            }
//            LocalExcelFile.writeArray2XLSXFile(fleOutput.getCanonicalPath(), titles, arrays);
//        }
//        hmpRGroup.clear();
//        String infor = "The association analysis results of geneset-base association are saved in " + fleOutput.getCanonicalPath();
//        LOG.info(infor);
>>>>>>> origin/master
        //remove the temporary products. 
        File fleRM = new File(rvtestOutputNameTmp + ".set.log");
        if (fleRM.exists()) {
            fleRM.delete();
        }
<<<<<<< HEAD
        if (!keepGrp) {
=======
        if (!boolKeep) {
>>>>>>> origin/master
            fleRM = new File(outputFileFolder + ".geneset.rvtest.grp.gz");
            if (fleRM.exists()) {
                fleRM.delete();
            }
        }
<<<<<<< HEAD
        return altGroupTitle;
=======
>>>>>>> origin/master
    }

    public void buildMap(HashMap<String, String[]> hmpRV, boolean boolKeep, File fleInput, int intP1, String strP2, int sg, String strMarker) throws FileNotFoundException, IOException {
        if (!fleInput.exists()) {
            return;
        }
        BufferedReader br = new BufferedReader(new FileReader(fleInput));
        //System.out.println(fleInput.getCanonicalFile());
        String strLine = br.readLine();
        String[] strHead = strLine.split("\t");
        int intP2;
        for (intP2 = 0; intP2 < strHead.length; intP2++) {
            if (strHead[intP2].contains(strP2)) {
                break;
            }
        }
        int intP3;
        for (intP3 = 0; intP3 < strHead.length; intP3++) {
            if (strHead[intP3].contains("CHROM")) {
                break;
            }
        }
        String key;
        while ((strLine = br.readLine()) != null) {
            String strItem[] = strLine.split("\t");
            if (!strMarker.equals("var")) {
                key = strItem[0];
            } else {
                key = strItem[intP3] + ":" + strItem[intP3 + 1];
            }
            if (hmpRV.containsKey(key)) {
                try {
                    if (!hmpRV.get(key)[intP1].isEmpty()) {
                        continue;
                    }
//                    hmpRV.get(key)[intP1]+=strItem[intP2]+";";
                    hmpRV.get(key)[intP1] = strItem[intP2];
                } catch (Exception e) {
                    hmpRV.get(key)[intP1] = ".";
                }
            } else {
                String[] intItem = new String[sg];
                Arrays.fill(intItem, "");
                try {
//                    intItem[intP1] = strItem[intP2]+";";
                    intItem[intP1] = strItem[intP2];
                } catch (Exception e) {
                    intItem[intP1] = ".";
                }
                hmpRV.put(key, intItem);
            }
        }
        br.close();
        if (!boolKeep) {
            fleInput.delete();
        }
<<<<<<< HEAD

    }

    public void releaseHmpRGroup() {
        this.hmpRGroup.clear();
    }

    public void releaseHmpRSingle() {
        this.hmpRSingle.clear();
=======
        return;
>>>>>>> origin/master
    }

    public class CallRVTestTask extends Task implements Callable<String>, Constants {

        String phenoFilePath;
        String chrName;
        String strMarker;

        public CallRVTestTask(String phenoFilePath, String chrName, String strMarker) {
            this.phenoFilePath = phenoFilePath;
            this.chrName = chrName;
            this.strMarker = strMarker;
        }

        @Override
        public String call() throws Exception {
<<<<<<< HEAD

            File f1 = new File(phenoFilePath);
            phenoFilePath = f1.getCanonicalPath();
//                String[] params = new String[15];
            String[] params;
            if (strMarker.contains("variant")) {
                params = new String[8];
                params[0] = rvtestFolder;
                params[1] = " --inVcf";
                params[2] = inputFileFolder;
                params[3] = "--pheno";
                params[4] = phenoFilePath;
                params[5] = "--out";
                params[6] = rvtestOutputNameTmp + ".chr" + chrName;
                params[7] = getCmd(strMarker);
            } else {
                params = new String[10];
                params[0] = rvtestFolder;
                params[1] = " --inVcf";
                params[2] = inputFileFolder;
                params[3] = "--pheno";
                params[4] = phenoFilePath;
                params[5] = "--setFile";
                params[6] = outputFileFolder + ".chr" + chrName + ".gene.rvtest.grp.gz";
                params[7] = "--out";
                params[8] = rvtestOutputNameTmp + ".chr" + chrName;
                //                params[9] = "--burden";
                //                params[10] = "cmc";
                //                params[11] = "--vt";
                //                params[12] = "price";
                //                params[13] = "--kernel";
                //                params[14] = "skat[nPerm=1000:alpha=0.001:beta1=1:beta2=20],kbac";

                params[9] = getCmd(strMarker);
            }

            StringBuilder comInfor = new StringBuilder();
            for (String param : params) {
                comInfor.append(param);
                comInfor.append(" ");
            }
            if (strPhe != null) {
                comInfor.append(strPhe);
                comInfor.append(" ");
            }
            if (strCov != null) {
                comInfor.append(strCov);
            }
            //example command
            //plugin/rvtests/executable/rvtest --inVcf test1.flt.vcf.gz --pheno assoc.ped --out output --setFile test1.gene.rvtest.grp.gz --burden cmc --vt price --kernel skat[nPerm=100:alpha=0.001:beta1=1:beta2=20],kbac
            LOG.info("Running RVTests:\n" + comInfor.toString());
            Process pr = Runtime.getRuntime().exec(comInfor.toString());
            /*
             try
            {            
                FileOutputStream fos = new FileOutputStream(rvtestOutputNameTmp+".out");

                // any error message?
                StreamGobbler errorGobbler = new 
                    StreamGobbler(pr.getErrorStream(), "ERROR");            

                // any output?
                StreamGobbler outputGobbler = new 
                    StreamGobbler(pr.getInputStream(), "OUTPUT", fos);

                // kick them off
                errorGobbler.start();
                outputGobbler.start();

                // any error???
                int exitVal = pr.waitFor();
                System.out.println("ExitValue: " + exitVal);
                fos.flush();
                fos.close();        
            } catch (Throwable t)
              {
                t.printStackTrace();
              }                    
             */

            String line;
            StringBuilder errorMsg = new StringBuilder();
            try {
                try (BufferedReader inputOut = new BufferedReader(new InputStreamReader(pr.getInputStream()))) {
                    while (((line = inputOut.readLine()) != null)) {
                        // System.out.println(line);
                    }
                }

                int exitVal;
                try (BufferedReader inputError = new BufferedReader(new InputStreamReader(pr.getErrorStream()))) {
                    while (((line = inputError.readLine()) != null)) {
                        errorMsg.append(line);
                        errorMsg.append("\n");
                    }
                }
                exitVal = pr.waitFor();
                pr.destroy();
                if (exitVal != 0) {
                    LOG.info("Rvtest failed to run by the command: " + comInfor.toString() + "\n" + errorMsg);
                }
            } catch (Exception ex) {
                LOG.error(ex + "\nRvtest failed to run!");
            }

            return "";
        }
=======
            try {
                File f1 = new File(phenoFilePath);
                phenoFilePath = f1.getCanonicalPath();
//                String[] params = new String[15];
                String[] params;
                if (strMarker.contains("variant")) {
                    params = new String[8];
                    params[0] = rvtestFolder;
                    params[1] = " --inVcf";
                    params[2] = inputFileFolder;
                    params[3] = "--pheno";
                    params[4] = phenoFilePath;
                    params[5] = "--out";
                    params[6] = rvtestOutputNameTmp + ".chr" + chrName;
                    params[7] = getCmd(strMarker);
                } else {
                    params = new String[10];
                    params[0] = rvtestFolder;
                    params[1] = " --inVcf";
                    params[2] = inputFileFolder;
                    params[3] = "--pheno";
                    params[4] = phenoFilePath;
                    params[5] = "--setFile";
                    params[6] = outputFileFolder + ".chr" + chrName + ".gene.rvtest.grp.gz";
                    params[7] = "--out";
                    params[8] = rvtestOutputNameTmp + ".chr" + chrName;
                    //                params[9] = "--burden";
                    //                params[10] = "cmc";
                    //                params[11] = "--vt";
                    //                params[12] = "price";
                    //                params[13] = "--kernel";
                    //                params[14] = "skat[nPerm=1000:alpha=0.001:beta1=1:beta2=20],kbac";

                    params[9] = getCmd(strMarker);
                }

                StringBuilder comInfor = new StringBuilder();
                for (String param : params) {
                    comInfor.append(param);
                    comInfor.append(" ");
                }
                if (strPhe != null) {
                    comInfor.append(strPhe);
                    comInfor.append(" ");
                }
                if (strCov != null) {
                    comInfor.append(strCov);
                }
                //example command
                //plugin/rvtests/executable/rvtest --inVcf test1.flt.vcf.gz --pheno assoc.ped --out output --setFile test1.gene.rvtest.grp.gz --burden cmc --vt price --kernel skat[nPerm=100:alpha=0.001:beta1=1:beta2=20],kbac
                LOG.info("Running RVTests:\n" + comInfor.toString());
                Process pr = Runtime.getRuntime().exec(comInfor.toString());

                String line;
                StringBuilder errorMsg = new StringBuilder();
                try {
                    BufferedReader inputOut = new BufferedReader(new InputStreamReader(pr.getInputStream()));
                    while (((line = inputOut.readLine()) != null)) {
                        //  System.out.println(line);
                    }

                    int exitVal;
                    try (BufferedReader inputError = new BufferedReader(new InputStreamReader(pr.getErrorStream()))) {
                        while (((line = inputError.readLine()) != null)) {
                            errorMsg.append(line);
                            errorMsg.append("\n");
                        }
                        exitVal = pr.waitFor();
                        pr.destroy();
                    }

                    if (exitVal != 0) {
                        LOG.info("Rvtest failed to run by the command: " + comInfor.toString() + "\n" + errorMsg);
                    }
                } catch (Exception ex) {
                    LOG.error(ex + "\nRvtest failed to run!");
                }

            } catch (IOException ex) {
                LOG.error(ex);
            } catch (Exception ex) {
                LOG.error(ex);
            }
            return "";
        }

>>>>>>> origin/master
    }

    public void runGeneAssoc(String phenoFilePath, int maxThreadNum, String strMarker) throws Exception {
        ExecutorService exec = Executors.newFixedThreadPool(maxThreadNum);
        final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
        int runningThread = 0;
        if (strMarker.contains("gene")) {
<<<<<<< HEAD
            for (String STAND_CHROM_NAMES1 : STAND_CHROM_NAMES) {
                File f = new File(outputFileFolder + ".chr" + STAND_CHROM_NAMES1 + ".gene.rvtest.grp.gz");
                if (!f.exists()) {
                    continue;
                }
                CallRVTestTask task = new CallRVTestTask(phenoFilePath, STAND_CHROM_NAMES1, strMarker);
=======
            for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
                File f = new File(outputFileFolder + ".chr" + STAND_CHROM_NAMES[i] + ".gene.rvtest.grp.gz");

                if (!f.exists()) {
                    continue;
                }

                CallRVTestTask task = new CallRVTestTask(phenoFilePath, STAND_CHROM_NAMES[i], strMarker);
>>>>>>> origin/master
                serv.submit(task);
                runningThread++;
            }
        } else {
<<<<<<< HEAD
            for (String STAND_CHROM_NAMES1 : STAND_CHROM_NAMES) {
                File f = new File(outputFileFolder + ".chr" + STAND_CHROM_NAMES1 + ".var.rvtest.grp.gz");
                if (!f.exists()) {
                    continue;
                }
                CallRVTestTask task = new CallRVTestTask(phenoFilePath, STAND_CHROM_NAMES1, strMarker);
=======
            for (int i = 0; i < STAND_CHROM_NAMES.length; i++) {
                File f = new File(outputFileFolder + ".chr" + STAND_CHROM_NAMES[i] + ".var.rvtest.grp.gz");

                if (!f.exists()) {
                    continue;
                }

                CallRVTestTask task = new CallRVTestTask(phenoFilePath, STAND_CHROM_NAMES[i], strMarker);
>>>>>>> origin/master
                serv.submit(task);
                runningThread++;
            }
        }

        for (int s = 0; s < runningThread; s++) {
            Future<String> task = serv.take();
            String infor = task.get();
            //  System.out.println(infor);
        }
        exec.shutdown();
        //collect data
    }

    public void generateGenesetAssocGroup(Map<String, List<Variant>> geneVars, BufferedWriter bw, boolean hasChrLabel) {
        String chrName;
        String lb;
        try {
            for (Map.Entry<String, List<Variant>> gVars : geneVars.entrySet()) {
                List<Variant> snps = gVars.getValue();
                int size = snps.size();
                if (size == 0) {
                    continue;
                }

                /*
                bw.write(gVars.getKey());
                bw.write(" ");
                //bw.write(snps.get(0).label);                
                chrName = STAND_CHROM_NAMES[snps.get(0).chrID];
                if (hasChrLabel) {
                    chrName = "chr" + chrName;
                }
              
                lb = chrName + ":" + snps.get(0).refStartPosition + "-" + (snps.get(0).refStartPosition);
                bw.write(lb);
                for (int i = 1; i < size; i++) {
                    bw.write(",");
                    chrName = STAND_CHROM_NAMES[snps.get(i).chrID];
                    lb = chrName + ":" + snps.get(i).refStartPosition + "-" + (snps.get(i).refStartPosition);
                    bw.write(lb);
                }
                bw.write("\n");
                 */
                for (int i = 0; i < size; i++) {
                    bw.write(gVars.getKey());
                    bw.write(" ");
                    chrName = STAND_CHROM_NAMES[snps.get(i).chrID];
                    if (hasChrLabel) {
                        chrName = "chr" + chrName;
                    }
                    lb = chrName + ":" + snps.get(i).refStartPosition + "-" + (snps.get(i).refStartPosition);
                    bw.write(lb);
                    bw.write("\n");
                }
            }
        } catch (Exception ex) {
            LOG.error(ex);
        }

    }

    public void runGenesetAssoc(String phenoFilePath, File groupFile) throws Exception {

        File f1 = new File(phenoFilePath);
        phenoFilePath = f1.getCanonicalPath();
        String[] params = new String[10];
        params[0] = rvtestFolder;
        params[1] = " --inVcf";
        params[2] = inputFileFolder;
        params[3] = "--pheno";
        params[4] = phenoFilePath;
        params[5] = "--setFile";
        params[6] = groupFile.getCanonicalPath();
        params[7] = "--out";
        params[8] = rvtestOutputNameTmp + ".set";
//        params[9] = "--burden";
//        params[10] = "cmc";
//        params[11] = "--vt";
//        params[12] = "price";
//        params[13] = "--kernel";
//        params[14] = "skat[nPerm=1000:alpha=0.001:beta1=1:beta2=20],kbac";
        params[9] = getCmd("geneset");
        StringBuilder comInfor = new StringBuilder();
        for (String param : params) {
            comInfor.append(param);
            comInfor.append(" ");
        }
        if (strPhe != null) {
            comInfor.append(strPhe);
            comInfor.append(" ");
        }
        if (strCov != null) {
            comInfor.append(strCov);
        }
        //example command
        //plugin/rvtests/executable/rvtest --inVcf test1.flt.vcf.gz --pheno assoc.ped --out output --setFile test1.gene.rvtest.grp.gz --burden cmc --vt price --kernel skat[nPerm=100:alpha=0.001:beta1=1:beta2=20],kbac
        //System.out.println(comInfor.toString());
        LOG.info("Running RVTests:\n" + comInfor.toString());
        Process pr = Runtime.getRuntime().exec(comInfor.toString());
<<<<<<< HEAD
        /*
        try {
         
            FileOutputStream fos = new FileOutputStream(rvtestOutputNameTmp + ".out");

            // any error message?
            StreamGobbler errorGobbler = new StreamGobbler(pr.getErrorStream(), "ERROR");

            // any output?
            StreamGobbler outputGobbler = new StreamGobbler(pr.getInputStream(), "OUTPUT", fos);

            // kick them off
            errorGobbler.start();
            outputGobbler.start();

            // any error???
            int exitVal = pr.waitFor();
            System.out.println("ExitValue: " + exitVal);
            fos.flush();
            fos.close();
        } catch (Throwable t) {
            t.printStackTrace();
        }
         */
=======

>>>>>>> origin/master
        String line;
        StringBuilder errorMsg = new StringBuilder();
        try {
            BufferedReader inputOut = new BufferedReader(new InputStreamReader(pr.getInputStream()));
            while (((line = inputOut.readLine()) != null)) {
<<<<<<< HEAD
                // System.out.println(line);
=======
                //  System.out.println(line);
>>>>>>> origin/master
            }

            int exitVal;
            try (BufferedReader inputError = new BufferedReader(new InputStreamReader(pr.getErrorStream()))) {
                while (((line = inputError.readLine()) != null)) {
                    errorMsg.append(line);
                    errorMsg.append("\n");
                }
<<<<<<< HEAD
            }
            exitVal = pr.waitFor();
            pr.destroy();
=======
                exitVal = pr.waitFor();
                pr.destroy();
            }

>>>>>>> origin/master
            if (exitVal != 0) {
                LOG.info("Rvtest failed to run by the command: " + comInfor.toString() + "\n" + errorMsg);
            }
        } catch (Exception ex) {
            LOG.error(ex + "\nRvtest failed to run!");
        }
        //collect data
    }

<<<<<<< HEAD
    public void summarizeVarCountsBySubject(Map<String, List<Variant>> geneVars,  List<Individual> subjectList, int[] pedEncodeGytIDMap, boolean isPhasedGty,
=======
    public void summarizeVarCountsBySubject(Map<String, List<Variant>> geneVars, OpenLongObjectHashMap wahBit, List<Individual> subjectList, int[] pedEncodeGytIDMap, boolean isPhasedGty,
>>>>>>> origin/master
            Map<String, Integer> phenotypeColID, String exportPath, int scoreIndex, boolean outGZ) throws Exception {
        BufferedWriter bwPed = null;
        if (geneVars == null) {
            return;
        }
        if (outGZ) {
            bwPed = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(exportPath + ".gz"))));
        } else {
            bwPed = new BufferedWriter(new FileWriter(exportPath));
        }

        int[] gtys = null;
        int alleleNum, base = 2;
        int subID = -1;
        int gtyID = 0;
        double count;
        bwPed.write("SubjectID");
        String[] traitNames = new String[phenotypeColID.size()];
        for (Map.Entry<String, Integer> items : phenotypeColID.entrySet()) {
            traitNames[items.getValue()] = items.getKey();
        }
        for (int i = 0; i < traitNames.length; i++) {
            bwPed.write("\t");
            bwPed.write(traitNames[i]);
        }
        for (Map.Entry<String, List<Variant>> items : geneVars.entrySet()) {
            bwPed.write("\t");
            bwPed.write(items.getKey());
        }
        bwPed.write("\n");
        boolean[] bits = new boolean[32];
<<<<<<< HEAD
        int startIndex;

        /*
        //restore the compressed genotypes
        int bitNum, byteNum, code, i1, i2;
        byte[] tempB;
        for (Map.Entry<String, List<Variant>> items : geneVars.entrySet()) {
            List<Variant> vars = items.getValue();
            if (vars == null || vars.isEmpty()) {
                continue;
            }
            for (Variant var : vars) {
                if (var.compressedGtyLabel >= 0) {
                    alleleNum = var.getAltAlleles().length + 1;
                    if (isPhasedGty) {
                        base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
                    } else {
                        base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
                    }
                    bitNum = alleleNum * pedEncodeGytIDMap.length;
                    byteNum = (bitNum) / 8;
                    if (bitNum % 8 != 0) {
                        byteNum++;
                    }
                    tempB = new byte[byteNum];
                    for (int j = 0; j < var.compressedGtyLabel; j += 3) {
                        code = (((var.encodedGty[j] & 0xFF) << 16) | ((var.encodedGty[j + 1] & 0xFF) << 8) | (var.encodedGty[j] & 0xFF));
                        i1 = code / 8;
                        i2 = code / 8;
                        tempB[i1] = (byte) (tempB[i1] | GlobalManager.byteOpers[i2]);
                    }
                }
            }
        }
         */
=======
        long startIndex;
>>>>>>> origin/master
        for (Individual indiv : subjectList) {
            if (indiv == null) {
                continue;
            }

            subID++;
            gtyID = pedEncodeGytIDMap[subID];
            if (gtyID < 0) {
                continue;
            }

            bwPed.write(indiv.getLabelInChip());
            double[] traits = indiv.getTraits();
            for (int i = 0; i < traits.length; i++) {
                bwPed.write("\t");
                bwPed.write(String.valueOf(traits[i]));
            }

            for (Map.Entry<String, List<Variant>> items : geneVars.entrySet()) {
                bwPed.write("\t");
                List<Variant> vars = items.getValue();
                if (vars == null || vars.isEmpty()) {
                    bwPed.write("0");
                    continue;
                }
                count = 0;
                for (Variant var : vars) {
                    alleleNum = var.getAltAlleles().length + 1;
                    if (isPhasedGty) {
                        base = GlobalManager.phasedAlleleBitMap.get(alleleNum);
<<<<<<< HEAD
                        if (var.compressedGtyLabel > 0) {
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
                                startIndex += pedEncodeGytIDMap.length;
                            }
                            gtys = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, gtyID);
                        } else if (var.compressedGtyLabel == 0) {
                            Arrays.fill(bits, 0, base, false);
                            gtys = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, gtyID);
                        } else {
                            gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, gtyID, pedEncodeGytIDMap.length);
                        }
                    } else {
                        base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
                        if (var.compressedGtyLabel > 0) {
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
                                startIndex += pedEncodeGytIDMap.length;
                            }
                            gtys = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, gtyID);
                        } else if (var.compressedGtyLabel == 0) {
                            Arrays.fill(bits, 0, base, false);
                            gtys = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, gtyID);
=======
                        if (var.compressedGty) {
                            startIndex = var.encodedGtyIndex[0] + gtyID;
                            for (int i = 0; i < base; i++) {
                                bits[i] = wahBit.containsKey(startIndex);
                                startIndex += pedEncodeGytIDMap.length;
                            }
                            gtys = BinaryGtyProcessor.getPhasedGtyBool(bits, alleleNum, base, gtyID);
                        } else {
                            gtys = BinaryGtyProcessor.getPhasedGtyAt(var.encodedGty, alleleNum, base, gtyID, pedEncodeGytIDMap.length);
                        }

                    } else {
                        base = GlobalManager.unphasedAlleleBitMap.get(alleleNum);
                        if (var.compressedGty) {
                            startIndex = var.encodedGtyIndex[0] + gtyID;
                            for (int i = 0; i < base; i++) {
                                bits[i] = wahBit.containsKey(startIndex);
                                startIndex += pedEncodeGytIDMap.length;
                            }
                            gtys = BinaryGtyProcessor.getUnphasedGtyBool(bits, alleleNum, base, gtyID);
>>>>>>> origin/master
                        } else {
                            gtys = BinaryGtyProcessor.getUnphasedGtyAt(var.encodedGty, alleleNum, base, gtyID, pedEncodeGytIDMap.length);
                        }
                    }

                    if (gtys != null) {
                        if (gtys[0] != 0) {
                            //   count+=(Double.parseDouble(var.featureValues[scoreIndex]));
                            count++;
                        }
                        if (gtys[1] != 0) {
                            //   count+=(Double.parseDouble(var.featureValues[scoreIndex]));
                            count++;
                        }
                    }
                }
                bwPed.write(String.valueOf(count));
            }
            bwPed.write("\n");
        }
        bwPed.close();
    }
<<<<<<< HEAD

    class StreamGobbler extends Thread {

        InputStream is;
        String type;
        OutputStream os;

        StreamGobbler(InputStream is, String type) {
            this(is, type, null);
        }

        StreamGobbler(InputStream is, String type, OutputStream redirect) {
            this.is = is;
            this.type = type;
            this.os = redirect;
        }

        public void run() {
            try {
                PrintWriter pw = null;
                if (os != null) {
                    pw = new PrintWriter(os);
                }

                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);
                String line = null;
                while ((line = br.readLine()) != null) {
                    if (pw != null) {
                        pw.println(line);
                    }
                    System.out.println(type + ">" + line);
                }
                if (pw != null) {
                    pw.flush();
                }
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
        }
    }

    public HashMap<String, String[]> getHmpRGroup() {
        return hmpRGroup;
    }

    public HashMap<String, String[]> getHmpRSingle() {
        return hmpRSingle;
    }

=======
>>>>>>> origin/master
}
