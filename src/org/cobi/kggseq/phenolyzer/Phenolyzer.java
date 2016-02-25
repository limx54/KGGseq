/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.phenolyzer;

import java.io.IOException;
import org.apache.log4j.Logger;
import static org.cobi.kggseq.GlobalManager.*;
import org.cobi.kggseq.Options;


/**
 *
 * @author JiangLi
 */
public class Phenolyzer {
     private static final Logger LOG = Logger.getLogger(Phenolyzer.class);
    String strPerlPath;
    String strPhenolyzerPath=PLUGIN_PATH+"phenolyzer-master/disease_annotation.pl";
    String strCMD;

    public Phenolyzer(String strPerlPath, String strPhenolyzerPath) {
        this.strPerlPath = strPerlPath;
        this.strPhenolyzerPath = strPhenolyzerPath;
    }

    public Phenolyzer(String strPerlPath) {
        this.strPerlPath = strPerlPath;
    }
    
    public void runPhenolyzer(){
        try {
            Runtime.getRuntime().exec(strCMD);
        } catch (IOException ex) {
            LOG.error(ex);
           
        }
    }
    
    
}
