/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import org.apache.log4j.Logger;
import java.util.HashMap;
import java.util.List;
import static org.cobi.kggseq.GlobalManager.PLUGIN_PATH;
import static org.cobi.kggseq.GlobalManager.osCode;
import org.cobi.kggseq.entity.AnnotationSummarySet;
import org.cobi.kggseq.entity.Chromosome;
import org.cobi.kggseq.entity.Variant;
import org.cobi.util.file.LocalFileFunc;

/**
 *
 * @author JiangLi
 */
public class Phenolyzer {

    private static final Logger LOG = Logger.getLogger(Phenolyzer.class);
    String strPerlPath;
    String strCMD;
    String fleOutput = "./phenolyzer";
    String strPhenolyzerPath;
    HashMap<String, String> hmpPhenolyzer = null;
    BufferedReader br = null;
    BufferedWriter bw = null;
    List<String> altSearchTerms = null;
//    String strURL = KGGSeq_URL + "download/lib/phenolyzer-master.zip";

    public Phenolyzer(List<String> searchTerms, String outPath) {
        if (osCode == 1) {
            this.strPerlPath = "C:/Strawberry/perl/bin/perl.exe";
        } else if (osCode == 2) {
            this.strPerlPath = "perl";
        } else if (osCode == 3) {
            this.strPerlPath = "perl";
        } else if (osCode == 4) {
            this.strPerlPath = "perl";
        }
        fleOutput = outPath;
        this.strPhenolyzerPath = PLUGIN_PATH + "phenolyzer-master/disease_annotation.pl";
        altSearchTerms = searchTerms;
    }

    public void runPhenolyzer() {
        String strItems = "";
        StringBuilder comInfor = new StringBuilder();
        for (int i = 0; i < altSearchTerms.size(); i++) {
            strItems += altSearchTerms.get(i) + ";";
        }
        strItems = strItems.substring(0, strItems.length() - 1);
        File f = new File(fleOutput);
        if (!f.getParentFile().exists()) {
            f.getParentFile().mkdir();
        }
        String[] params = new String[12];
        params[0] = strPerlPath;
        params[1] = strPhenolyzerPath;
        params[2] = "\"" + strItems + "\"";
        params[3] = "-p";
        params[4] = "-ph";
        params[5] = "-logistic";
        params[6] = "-out";
        params[7] = fleOutput;
        params[8] = "-addon";
        params[9] = "DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE,DB_GENECARDS_GENE_DISEASE_SCORE";
        params[10] = "-addon_weight";
        params[11] = "0.25";

        // strCMD = strPerlPath + " " + strPhenolyzerPath + " \"" + strItems + "\" -p -ph -logistic -out " + fleOutput + " -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE,DB_GENECARDS_GENE_DISEASE_SCORE -addon_weight 0.25";
        // System.out.println(strCMD);
//            String strTest="C:\\Strawberry\\perl\\bin\\perl.exe D:\\01WORK\\KGGseq\\software\\KGGseq\\plugin\\phenolyzer-master\\disease_annotation.pl 'alzheimer;crohn' -p -ph -logistic -out testPhenolyzer_New\\testPhenolyzerphenolyzer\\out";
//            System.out.println(strTest);
        // Process pr = Runtime.getRuntime().exec(strCMD);
        //Process pr = new ProcessBuilder(strCMD).start();
        String line;

        try {
            Process pr = Runtime.getRuntime().exec(params);

            BufferedReader inputError = new BufferedReader(new InputStreamReader(pr.getErrorStream()));
            while (((line = inputError.readLine()) != null)) {
                //  System.out.println(line);
            }

            int exitVal = pr.waitFor();
            pr.destroy();
            inputError.close();

            if (exitVal != 0) {
                for (String param : params) {
                    comInfor.append(param);
                    comInfor.append(" ");
                }
                LOG.info("Phenolyzer failed to run by the command: " + comInfor.toString());
            }
        } catch (Exception ex) {
            LOG.error("Phenolyzer failed to run by the command: " + comInfor.toString() + "\n" + ex.toString());
        }

    }

    public void parseResult() {
        try {
//            if(fleInput==null || !fleInput.exists()) return;
            hmpPhenolyzer = new HashMap<String, String>();
            br = new BufferedReader(new FileReader(fleOutput + ".final_gene_list"));

            String strLine;
//            while((strLine=br.readLine())!=null){
//                if(strLine.contains("Normalized score")){
//                    String[] strItems=strLine.split("\t");
//                    double[] dblScore=new double[2];
//                    dblScore[0]=Double.parseDouble(strItems[1].substring(18));
//                    dblScore[1]=Double.parseDouble(strItems[2].substring(11));
//                    hmpPhenolyzer.put(strItems[0], dblScore);
//                }
//            }

            while ((strLine = br.readLine()) != null) {
                String[] strItems = strLine.split("\t");
                hmpPhenolyzer.put(strItems[1], strItems[3]);
            }
            br.close();
            File f = new File(fleOutput);

            if (f.getParentFile().exists()) {
                LocalFileFunc.delAll(f.getParentFile());
            }
        } catch (FileNotFoundException ex) {
            LOG.info(ex);

        } catch (IOException ex) {
            LOG.info(ex);
        }
    }

    public HashMap<String, String> getHashMap() {
        return this.hmpPhenolyzer;
    }

    public void addScore(Chromosome chromosome, AnnotationSummarySet ass, HashMap<String, String> hmpPhenolyzer) {
        int intNum = 0;
        if (chromosome == null) {
            return;
        }

        int varFeatureNum = ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            if (var.geneSymb == null) {
                var.setFeatureValue(varFeatureNum, null);
//                var.setFeatureValue(ass.getAvailableFeatureIndex() + 1, null);
                continue;
            }
            String strGene = var.geneSymb.toUpperCase();
            if (hmpPhenolyzer.containsKey(strGene)) {
                var.setFeatureValue(varFeatureNum, hmpPhenolyzer.get(strGene));
                intNum++;
            } else {
                var.setFeatureValue(varFeatureNum, ".");
            }
        }

        ass.setAnnotNum(ass.getAnnotNum() + intNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.mRNAList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.mRNAList.size() - intNum);
    }

}
