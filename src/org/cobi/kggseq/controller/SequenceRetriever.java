/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.cobi.kggseq.Constants;
import org.cobi.kggseq.GlobalManager;
import org.cobi.kggseq.Options;
import org.cobi.kggseq.entity.AnnotationSummarySet;
import org.cobi.kggseq.entity.Chromosome;
import org.cobi.kggseq.entity.Genome;
import org.cobi.kggseq.entity.Variant;
import org.cobi.util.net.NetUtils;
import org.cobi.util.text.LocalFile;

/**
 *
 * @author mxli
 */
public class SequenceRetriever implements Constants {

    public void addFlankingSequences(Chromosome chromosome, AnnotationSummarySet ass, int extendLen, String refGenomeVersion) throws Exception {
        int hitNum=0;
        if (chromosome == null) {
            return;
        }
        
        File fastAFile = new File(GlobalManager.RESOURCE_PATH + "/" + refGenomeVersion + "/chr" + chromosome.getName() + ".fa.gz");
        StringBuilder sequences = LocalFile.retrieveData(fastAFile.getCanonicalPath(), 1);
        StringBuilder flankSeq = new StringBuilder();
        
        int varFeatureNum=ass.getAvailableFeatureIndex();
        for (Variant var : chromosome.variantList) {
            flankSeq.delete(0, flankSeq.length());
            int left = var.refStartPosition - extendLen - 1;
            if (left < 0) {
                left = 0;
            }
            flankSeq.append(sequences.substring(left, var.refStartPosition - 1));
            flankSeq.append('[');
            flankSeq.append(var.getRefAllele());

            for (String seq : var.getAltAlleles()) {
                flankSeq.append('/');
                flankSeq.append(seq);
            }
            flankSeq.append(']');
            int right = var.refStartPosition + extendLen;

            if (right > sequences.length()) {
                right = sequences.length();
            }
            flankSeq.append(sequences.substring(var.refStartPosition, right));
//            var.addFeatureValue(flankSeq.toString());
            var.setFeatureValue(varFeatureNum, flankSeq.toString());
            hitNum++;
        }
        ass.setAnnotNum(ass.getAnnotNum() + hitNum);
        ass.setTotalNum(ass.getTotalNum() + chromosome.variantList.size());
        ass.setLeftNum(ass.getLeftNum() + chromosome.variantList.size() - hitNum);
    }
}
