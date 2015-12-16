/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.kggseq.controller;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
//import net.sf.picard.liftover.LiftOver;
//import net.sf.picard.util.Interval;
//import net.sf.samtools.util.AsciiLineReader;
//import net.sf.samtools.util.CompressedFileReader;
//import net.sf.samtools.util.LineReader;
import org.apache.log4j.Logger;
import org.cobi.kggseq.entity.Genome;
import org.cobi.kggseq.entity.Variant;
import org.cobi.util.file.LocalFileFunc;
import org.cobi.util.text.Util;

/**
 *
 * @author mxli
 */
public class PileupFormatParser {

    //  private static final Log LOG = Log.getInstance(PileupFormatParser.class);
    private static final Logger LOG = Logger.getLogger(PileupFormatParser.class);
    /*
     * Function to read the output by samtool pileup command.
     An example of the comand is like:
     * ./samtools pileup -vcf chr1.fa A_1.chr1.sorted.bam > A_1.chr1.gty.raw
     ./samtools.pl varFilter -d 8  A_1.chr1.gty.raw | awk '$5>=30' > A_1.chr1.gty
    
     * The output
     chr1	59133	A	G	125	141	36	39	GGGGGGGGGGGGGG.GGGGGGGGGGGGGGGGGGG^FG^HG^6G^LG^HG	E#-D?A3??=EBEDFGB@EGEBCDGBFF5GFFFDBGDGG
     chr1	716832	g	S	63	63	57	14	....,,cc,c,ccc	:D=E@BBE@CFC?-
     chr1	716890	C	A	25	54	55	10	,aaaaaaaaa	GGEGDGGGG:
     chr1	798785	G	A	175	175	52	49	AAaAaAAAAAAAAAAAAAAAaaaaaaAaaAaaaAAaaaaaaAAaaaaaA	-#G?G>@9@##@?CC#AC@#?E>ACEB==EC0:B?BB?=BCGE@@?=AG
     chr1	798791	C	Y	173	173	52	52	Tt.t..TTTTTTT......t,,,,,Ttt.tt,.T,,,,tt.Tttt,,.,,..	#G?G?>A?##C??D#AC?#CDC=GFC=?GC?:FDB@?BDCGE??@=BG##GC
    
     */

    /**
     * @return @throws Exception
     * @pdOid f0621cff-9d97-421e-a77a-765bd0938dfb
     */
//    public Genome readVariantUnique(String vAFile, double avgSeqQualityThrehsold, double gtyQualityThrehsold, int coverageT, boolean needProgressionIndicator) throws Exception {
//        int indexCHROM = 0;
//        int indexPOS = 1;
//        int indexREF = 2;
//        int indexALT = 3;
//        int indexSeqQuality = 4;
//        int indexGtyQuality = 5;
//        int indexMapQuality = 6;
//        int indexDP = 7;
//        int mapQualityThrehsold = 10;
//        
//        int maxColNum = indexCHROM;
//        maxColNum = Math.max(maxColNum, indexPOS);
//        maxColNum = Math.max(maxColNum, indexREF);
//        maxColNum = Math.max(maxColNum, indexALT);
//        maxColNum = Math.max(maxColNum, indexSeqQuality);
//        maxColNum = Math.max(maxColNum, indexGtyQuality);
//        maxColNum = Math.max(maxColNum, indexMapQuality);
//        maxColNum = Math.max(maxColNum, indexDP);
//        
//        String currentLine = null;
//        String currChr = null;
//        StringBuilder tmpBuffer = new StringBuilder();
//        long lineCounter = 0;
//        
//        File dataFile = new File(vAFile);
//        if (!dataFile.exists()) {
//            throw new Exception("No such a pileup file: " + dataFile.getCanonicalPath());
//        }
//        
//        BufferedReader br = null;
//        if (dataFile.exists() && dataFile.getName().endsWith(".zip")) {// should we support zip format?
//            br = LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
//        } else {
//            if (dataFile.exists() && dataFile.getName().endsWith(".gz")) {
//                br = LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
//            } else {
//                if (dataFile.exists()) {
//                    br = LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
//                } else {
//                    throw new Exception("No input file: " + dataFile.getCanonicalPath());
//                }
//            }
//        }
//        int charIndex1 = 0;
//        int charIndex2 = 0;
//        int makerPostion = 0;
//        
//        String ref = ".";
//        String alt = null;
//        boolean incomplete = true;
//        int acceptVarNum = 0;
//        
//        int ignoredLowQualGtyNum = 0;
//        double avgSeqQuality = 0;
//        double gtyQuality = 0;
//        double mapQuality = 0;
//        int coverage = 0;
//        int iCol = 0;
//        final String REF_HOM = "0/0";
//        final CharSequence REF_ALLE = "0";
//        //temp variables
//        Genome genome = new Genome(dataFile.getName());
//        int indelNum = 0;
//        int index1 = 0;
//        int index2 = 0;
//        int colonNum = 0;
//        int len = tmpBuffer.length();
//        
//        try {
//            iCol = 0;
//
//            //  bw.write(currentLine);                bw.write("\n");
//            while ((currentLine = br.readLineS()) != null) {
//                lineCounter++;
//                if (needProgressionIndicator && lineCounter % 10000 == 0) {
//                    String prog = String.valueOf(lineCounter);
//                    System.out.print(prog);
//                    char[] backSpaces = new char[prog.length()];
//                    Arrays.fill(backSpaces, '\b');
//                    System.out.print(backSpaces);
//                }
//                StringTokenizer st = new StringTokenizer(currentLine.trim());
//                //initialize varaibles
//                incomplete = true;
//                coverage = 0;
//                avgSeqQuality = 0;
//                currChr = null;
//                makerPostion = -1;
//                
//                ref = ".";
//                alt = null;
//
////  chr1	59133	A	G	125	141	36	39	GGGGGGGGGGGGGG.GGGGGGGGGGGGGGGGGGG^FG^HG^6G^LG^HG	E#-D?A3??=EBEDFGB@EGEBCDGBFF5GFFFDBGDGG
//                //  chr1	716832	g	S	63	63	57	14	....,,cc,c,ccc	:D=E@BBE@CFC?-
//                //  chr1	716890	C	A	25	54	55	10	,aaaaaaaaa	GGEGDGGGG:
//                for (iCol = 0; iCol <= maxColNum; iCol++) {
//                    if (st.hasMoreTokens()) {
//                        tmpBuffer.delete(0, tmpBuffer.length());
//                        tmpBuffer.append(st.nextToken().trim());
//                        if (iCol == indexCHROM) {
//                            currChr = tmpBuffer.toString();
//                        } else if (iCol == indexPOS) {
//                            makerPostion = Util.parseInt(tmpBuffer.toString());
//                        } else if (iCol == indexREF) {
//                            ref = tmpBuffer.toString();
//                        } else if (iCol == indexALT) {
//                            alt = tmpBuffer.toString();
//                        } else if (iCol == indexSeqQuality) {
//                            avgSeqQuality = Double.parseDouble(tmpBuffer.toString());
//                        } else if (iCol == indexGtyQuality) {
//                            gtyQuality = Double.parseDouble(tmpBuffer.toString());
//                        } else if (iCol == indexDP) {
//                            coverage = Util.parseInt(tmpBuffer.toString());
//                        }
//                    } else {
//                        break;
//                    }
//                    if (iCol >= maxColNum) {
//                        incomplete = false;
//                    }
//                }
//                
//                if (incomplete) {
//                    continue;
//                }
//                
//                if (ref.equals("*")) {
//                    indelNum++;
//                    continue;
//                }
//                if (avgSeqQuality < avgSeqQualityThrehsold) {
//                    ignoredLowQualGtyNum++;
//                    continue;
//                }
//                if (coverage < coverageT) {
//                    ignoredLowQualGtyNum++;
//                    continue;
//                }
//                if (gtyQuality < gtyQualityThrehsold) {
//                    ignoredLowQualGtyNum++;
//                    continue;
//                }
//                
//                switch (Character.toUpperCase(alt.charAt(0))) {
//                    case 'A': //AA 0010
//                        alt = "A";
//                        break;
//                    case 'T': //TT 0011
//                        alt = "T";
//                        break;
//                    case 'G': //GG 0100
//                        alt = "G";
//                        break;
//                    case 'C': //CC 0101
//                        alt = "C";
//                        break;
//                    case 'R': // R AG 0110
//                        if (ref.charAt(0) == 'A') {
//                            alt = "G";
//                        } else {
//                            alt = "A";
//                        }
//                        break;
//                    case 'Y': //YCT 0111
//                        if (ref.charAt(0) == 'C') {
//                            alt = "T";
//                        } else {
//                            alt = "C";
//                        }
//                        break;
//                    case 'K': //KGT 1000
//                        if (ref.charAt(0) == 'G') {
//                            alt = "T";
//                        } else {
//                            alt = "G";
//                        }
//                        break;
//                    case 'M': //MAC 1001
//                        if (ref.charAt(0) == 'A') {
//                            alt = "C";
//                        } else {
//                            alt = "A";
//                        }
//                        break;
//                    case 'S': //SGC 1010
//                        if (ref.charAt(0) == 'G') {
//                            alt = "C";
//                        } else {
//                            alt = "G";
//                        }
//                        break;
//                    case 'W': //WAT  1011
//                        if (ref.charAt(0) == 'A') {
//                            alt = "T";
//                        } else {
//                            alt = "A";
//                        }
//                        break;
//                    default:
//                        String infor = "Unknown gentype " + alt + " when parsing at line " + lineCounter + ": " + currentLine;
//                        System.out.println(infor);
//                        continue;
//                }
//
//                //Note: without pedigree information, I do not know who are affected and unaffected
//                // so this  can be changed for real case unique or control uniuqe variants
//                Variant var = new Variant(makerPostion, ref, new String[]{alt});
//                genome.addVariantFullChromName(currChr, var);
//                //   bw.write(currentLine);
//                //    bw.write("\n");
//                acceptVarNum++;
//            }
//            
//        } catch (NumberFormatException nex) {
//            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
//            // LOG.error(nex, info);
//            throw new Exception(info);
//        }
//        br.close();
//        //  bw.close();
//        String prog = String.valueOf(lineCounter);
//        System.out.println(prog);
//        genome.buildVariantIndexMapOnChromosomes();
//        
//        StringBuilder message = new StringBuilder();
//        
//        if (ignoredLowQualGtyNum > 0) {
//            message.append(ignoredLowQualGtyNum).append(" genotypes in valid variants are ignroed "
//                    + "due to sequencing quality < ").append(avgSeqQualityThrehsold).
//                    append(" or genotype quality < ").append(gtyQualityThrehsold).
//                    append(" or depth < ").append(coverageT).append('\n');
//        }
//        if (indelNum > 0) {
//            message.append(indelNum).append(" indels are not considered.").append('\n');
//        }
//        message.append('\n').append(lineCounter).append(" lines are scanned;  and ").append(acceptVarNum).append(" variants ").append(" are valid in ").append(vAFile);
//        LOG.info(message);
//        if (acceptVarNum == 0) {
//            System.exit(1);
//        }
//        return genome;
//    }
//    
//    public Genome readVariantSimpleFormat(String vAFile, boolean genotypeAlt, Set<String> caseSet, boolean needProgressionIndicator) throws Exception {
//        int indexCHROM = 0;
//        int indexPOS = 1;
//        int indexREF = 2;
//        
//        int maxColNum = indexCHROM;
//        maxColNum = Math.max(maxColNum, indexPOS);
//        maxColNum = Math.max(maxColNum, indexREF);
//        
//        String currentLine = null;
//        String currChr = null;
//        StringBuilder tmpBuffer = new StringBuilder();
//        long lineCounter = 0;
//        
//        File dataFile = new File(vAFile);
//        if (!dataFile.exists()) {
//            throw new Exception("No such a file: " + dataFile.getCanonicalPath());
//        }
//        LineReader br = null;
//        if (dataFile.exists() && dataFile.getName().endsWith(".zip")) {
//            br = new CompressedFileReader(dataFile);
//        } else {
//            if (dataFile.exists() && dataFile.getName().endsWith(".gz")) {
//                br = new CompressedFileReader(dataFile);
//            } else {
//                if (dataFile.exists()) {
//                    br = new AsciiLineReader(dataFile);
//                } else {
//                    throw new Exception("No input file: " + dataFile.getCanonicalPath());
//                }
//            }
//        }
//        
//        int makerPostion = 0;
//        
//        String ref = null;
//        String alt = null;
//        boolean incomplete = true;
//        int acceptVarNum = 0;
//        
//        Genome genome = new Genome(dataFile.getName());
//        
//        genome.addVariantFeatureLabel("AffectedRefHomGtyNum");
//        genome.addVariantFeatureLabel("AffectedHetGtyNum");
//        genome.addVariantFeatureLabel("AffectedAltHomGtyNum");
//        
//        int iCol = 0;
//        //temp variables
//        int indelNum = 0;
//        boolean isIndel = false;
//        int homoNum = 0;
//        
//        try {
//            iCol = 0;
////skip head line
//            currentLine = br.readLineS();
//            String[] cells = currentLine.split("\t");
//            int indivIndex = -1;
//            for (int i = 0; i < cells.length; i++) {
//                if (caseSet.contains(cells[i])) {
//                    indivIndex = i;
//                    break;
//                }
//            }
//            if (indivIndex == -1) {
//                throw (new Exception("No affected sample ID found " + caseSet.toString()));
//            }
//            maxColNum = Math.max(maxColNum, indivIndex);
//            //  bw.write(currentLine);                bw.write("\n");
//            while ((currentLine = br.readLineS()) != null) {
//                lineCounter++;
//                if (needProgressionIndicator && lineCounter % 10000 == 0) {
//                    String prog = String.valueOf(lineCounter);
//                    System.out.print(prog);
//                    char[] backSpaces = new char[prog.length()];
//                    Arrays.fill(backSpaces, '\b');
//                    System.out.print(backSpaces);
//                }
//                StringTokenizer st = new StringTokenizer(currentLine.trim());
//                //initialize varaibles
//                incomplete = true;
//                currChr = null;
//                makerPostion = -1;
//                
//                ref = null;
//                alt = null;
//                isIndel = false;
//                homoNum = 0;
//                //  chr1	59133	A	G	125	141	36	39	GGGGGGGGGGGGGG.GGGGGGGGGGGGGGGGGGG^FG^HG^6G^LG^HG	E#-D?A3??=EBEDFGB@EGEBCDGBFF5GFFFDBGDGG
//                //  chr1	716832	g	S	63	63	57	14	....,,cc,c,ccc	:D=E@BBE@CFC?-
//                //  chr1	716890	C	A	25	54	55	10	,aaaaaaaaa	GGEGDGGGG:
//
//                for (iCol = 0; iCol <= maxColNum; iCol++) {
//                    if (st.hasMoreTokens()) {
//                        tmpBuffer.delete(0, tmpBuffer.length());
//                        tmpBuffer.append(st.nextToken().trim());
//                        if (iCol == indexCHROM) {
//                            currChr = tmpBuffer.toString();
//                        } else if (iCol == indexPOS) {
//                            makerPostion = Util.parseInt(tmpBuffer.toString());
//                        } else if (iCol == indexREF) {
//                            ref = tmpBuffer.toString();
//                        } else if (iCol == indivIndex) {
//                            alt = tmpBuffer.toString();
//                        }
//                        
//                    } else {
//                        break;
//                    }
//                    if (iCol >= maxColNum) {
//                        incomplete = false;
//                    }
//                }
//                
//                if (incomplete) {
//                    continue;
//                }
//                if (alt.equals("N")) {
//                    continue;
//                }
//                if (ref.indexOf('-') >= 0 || alt.indexOf('-') >= 0) {
//                    indelNum++;
//                    isIndel = true;
//                }
//                if (genotypeAlt && !isIndel) {
//                    switch (Character.toUpperCase(alt.charAt(0))) {
//                        case 'A': //AA 0010
//                            alt = "A";
//                            homoNum = 1;
//                            break;
//                        case 'T': //TT 0011
//                            alt = "T";
//                            homoNum = 1;
//                            break;
//                        case 'G': //GG 0100
//                            alt = "G";
//                            homoNum = 1;
//                            break;
//                        case 'C': //CC 0101
//                            alt = "C";
//                            homoNum = 1;
//                            break;
//                        case 'R': // R AG 0110
//                            if (ref.charAt(0) == 'A') {
//                                alt = "G";
//                            } else {
//                                alt = "A";
//                            }
//                            break;
//                        case 'Y': //YCT 0111
//                            if (ref.charAt(0) == 'C') {
//                                alt = "T";
//                            } else {
//                                alt = "C";
//                            }
//                            break;
//                        case 'K': //KGT 1000
//                            if (ref.charAt(0) == 'G') {
//                                alt = "T";
//                            } else {
//                                alt = "G";
//                            }
//                            break;
//                        case 'M': //MAC 1001
//                            if (ref.charAt(0) == 'A') {
//                                alt = "C";
//                            } else {
//                                alt = "A";
//                            }
//                            break;
//                        case 'S': //SGC 1010
//                            if (ref.charAt(0) == 'G') {
//                                alt = "C";
//                            } else {
//                                alt = "G";
//                            }
//                            break;
//                        case 'W': //WAT  1011
//                            if (ref.charAt(0) == 'A') {
//                                alt = "T";
//                            } else {
//                                alt = "A";
//                            }
//                            break;
//                        default:
//                            String infor = "Unknown gentype " + alt + " when parsing at line " + lineCounter + ": " + currentLine;
//                            System.out.println(infor);
//                            continue;
//                    }
//                } else if (isIndel) {
//                    //the alt one is an allele not a genotype 
//                    makerPostion--;
//                    //insertion
//                    if (ref.charAt(0) == '-') {
//                        alt = "+" + alt;
//                    } else {
//                        //deletion
//                        alt = "-" + ref;
//                        ref = "-";
//                    }
//                }
//
//                //Note: without pedigree information, I do not know who are affected and unaffected
//                // so this  can be changed for real case unique or control uniuqe variants
//                Variant var = new Variant(makerPostion, ref, new String[]{alt});
//                var.setIsIndel(isIndel);
//                var.addFeatureValue("0");
//                var.addFeatureValue(String.valueOf(1 - homoNum));
//                var.addFeatureValue(String.valueOf(homoNum));
//                genome.addVariant(currChr, var);
//
//                //   bw.write(currentLine);
//                //    bw.write("\n");
//                acceptVarNum++;
//            }
//            
//        } catch (NumberFormatException nex) {
//            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
//            // LOG.error(nex, info);
//            throw new Exception(info);
//        }
//        br.close();
//        //  bw.close();
//        String prog = String.valueOf(lineCounter);
//        System.out.println(prog);
//        genome.buildVariantIndexMapOnChromosomes();
//        StringBuilder message = new StringBuilder();
//        genome.setVarNum(acceptVarNum);
//        
//        if (indelNum > 0) {
//            message.append(indelNum).append(" indels are not considered.").append('\n');
//        }
//        message.append('\n').append(lineCounter).append(" lines are scanned;  and ").append(acceptVarNum).append(" variants ").append(" are valid in ").append(vAFile);
//        LOG.info(message);
//        if (acceptVarNum == 0) {
//            System.exit(1);
//        }
//        return genome;
//    }
//    
//    public Genome readVariantMyFilterFormat(String vAFile, boolean needProgressionIndicator) throws Exception {
//        File dataFile = new File(vAFile);
//        Genome genome = new Genome(dataFile.getName());
//        
//        int indexChrom = 0;
//        int indexPosition = 1;
//        int indexREF = 2;
//        int indexALT = 3;
//        int indexMAF = 4;
//        
//        LineReader br = null;
//        if (dataFile.exists() && dataFile.getName().endsWith(".zip")) {
//            br = new CompressedFileReader(dataFile);
//        } else {
//            if (dataFile.exists() && dataFile.getName().endsWith(".gz")) {
//                br = new CompressedFileReader(dataFile);
//            } else {
//                if (dataFile.exists()) {
//                    br = new AsciiLineReader(dataFile);
//                } else {
//                    throw new Exception("No input file: " + dataFile.getCanonicalPath());
//                }
//            }
//        }
//        
//        int maxColNum = indexChrom;
//        maxColNum = Math.max(maxColNum, indexPosition);
//        maxColNum = Math.max(maxColNum, indexREF);
//        maxColNum = Math.max(maxColNum, indexALT);
//        maxColNum = Math.max(maxColNum, indexMAF);
//        
//        String currentLine;
//        int lineCounter = 0;
//        boolean incomplete = true;
//        String currChr = null;
//        int position = -1;
//        StringBuilder tmpBuffer = new StringBuilder();
//        String ref;
//        String alt;
//        float maf = 0;
//        boolean is1KGIndel = false;
//        
//        char[] backSpaces = null;
//        int varNum = 0;
//        while ((currentLine = br.readLineS()) != null) {
//            StringTokenizer st = new StringTokenizer(currentLine.trim());
//            //initialize varaibles
//            incomplete = true;
//            position = -1;
//            ref = null;
//            alt = null;
//            currChr = null;
//            for (int iCol = 0; iCol <= maxColNum; iCol++) {
//                if (st.hasMoreTokens()) {
//                    tmpBuffer.delete(0, tmpBuffer.length());
//                    tmpBuffer.append(st.nextToken().trim());
//                    if (iCol == indexChrom) {
//                        //the format is like this 1       469     C       G       0.150 
//                        currChr = tmpBuffer.toString();
//                    } else if (iCol == indexPosition) {
//                        position = Util.parseInt(tmpBuffer.toString());
//                    } else if (iCol == indexREF) {
//                        ref = tmpBuffer.toString();
//                    } else if (iCol == indexALT) {
//                        alt = tmpBuffer.toString();
//                    } else if (iCol == indexMAF) {
//                        maf = Util.parseFloat(tmpBuffer.toString());
//                    }
//                } else {
//                    break;
//                }
//                if (iCol == maxColNum) {
//                    incomplete = false;
//                }
//            }
//            if (incomplete) {
//                continue;
//            }
//
//            // if (!currChr.equals("6")) continue;
//            lineCounter++;
//            if (needProgressionIndicator && lineCounter % 20000 == 0) {
//                String prog = String.valueOf(lineCounter);
//                System.out.print(prog);
//                backSpaces = new char[prog.length()];
//                Arrays.fill(backSpaces, '\b');
//                System.out.print(backSpaces);
//            }
//            is1KGIndel = false;
//            //format:1	45113	-	0TATGG	0.715732
/////1	53599	CTA	3	0.890916
////1	223450	CT	2	0.207385
//            if (alt.charAt(0) - '0' <= 9 && alt.charAt(0) - '0' >= 0) {
//                is1KGIndel = true;
//                position--;
//            }
//
//            /*
//             Interval int1 = new Interval("chr" + currChr, position, position);
//             Interval int2 = liftOver.liftOver(int1);
//             if (int2 != null) {
//             position = int2.getStart();
//             } else {
//             failtedMapVarNum++;
//             continue;
//             }
//             * 
//             */
////Note: without pedigree information, I do not know who are affected and unaffected
//            // so this  can be changed for real case unique or control uniuqe variants
//            Variant var = new Variant(position, ref, new String[]{alt});
//            var.setIsIndel(is1KGIndel);
//            var.altAF = maf;
//            genome.addVariant(currChr, var);
//            varNum++;
//        }
//        br.close();
//        genome.buildVariantIndexMapOnChromosomes();
//        String info = varNum + " variant(s) are left after filtered by " + lineCounter + " variant(s) in " + dataFile.getName();
//        if (needProgressionIndicator) {
//            backSpaces = new char[7];
//            Arrays.fill(backSpaces, '\b');
//            System.out.print(backSpaces);
//        }
//        
//        LOG.info(info);
//        return genome;
//    }
//    
    public Genome readVariantAnnovarFormat(String vAFile, boolean needProgressionIndicator) {
        int indexCHROM = 0;
        int indexPosStart = 1;
        int indexPosEnd = 2;
        int indexREF = 3;
        int indexALT = 4;
        int indexComments = 5;

        int maxColNum = indexCHROM;
        maxColNum = Math.max(maxColNum, indexPosStart);
        maxColNum = Math.max(maxColNum, indexPosEnd);
        maxColNum = Math.max(maxColNum, indexREF);
        maxColNum = Math.max(maxColNum, indexALT);

        String currentLine = null;
        String currChr = null;
        StringBuilder tmpBuffer = new StringBuilder();
        long lineCounter = 0;
        int makerPostionStart = 0;
        int makerPostionEnd = 0;

        String ref = null;
        String alt = null;
        boolean incomplete = true;
        int acceptVarNum = 0;

        Genome genome = null;
        boolean hasCommentsCol = false;
        boolean hasHead = false;
        int iCol = 0;
        //temp variables
        int indelNum = 0;
        boolean isIndel = false;
        int homoNum = 0;

        StringBuilder sb = new StringBuilder();
        String comments = null;

        try {
            File dataFile = new File(vAFile);
            if (!dataFile.exists()) {
                throw new Exception("No ANNOVAR a file: " + dataFile.getCanonicalPath());
            }
            BufferedReader br = null;
            br = LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
//            if (dataFile.exists() && dataFile.getName().endsWith(".zip")) { // should we support .zip format
//                br = LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
//            } else {
//                if (dataFile.exists() && dataFile.getName().endsWith(".gz")) {
//                   br = LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
//                } else {
//                    if (dataFile.exists()) {
//                        br = LocalFileFunc.getBufferedReader(dataFile.getCanonicalPath());
//                    } else {
//                        throw new Exception("No input file: " + dataFile.getCanonicalPath());
//                    }
//                }
//            }

            genome = new Genome("AnnovarGenome", "annovar");
            genome.removeTempFileFromDisk();
            //currentLine = br.readLineS();
            currentLine = br.readLine();
            if (currentLine == null) {
                return genome;
            }
            char delimiter = '\t';
            int geneFeatureStartIndex = genome.getVariantFeatureNum();
            String[] cells = Util.tokenizeIngoreConsec(currentLine, delimiter);
            iCol = cells.length;
            if (iCol < 4) {
                delimiter = ' ';
                cells = Util.tokenizeIngoreConsec(currentLine, delimiter);
                iCol = cells.length;
            }
            if (iCol < 4) {
                String info = "Invalid format at line " + currentLine + " in the " + vAFile;
                throw new Exception(info);
            }

            if (iCol > 5) {
                maxColNum++;
                hasCommentsCol = true;
            }
            if (Util.isNumeric(cells[indexPosStart])) {
                if (hasCommentsCol) {
                    genome.addVariantFeatureLabel("Comments");
                }
            } else {
                //has head
                hasHead = true;
                for (int t = 5; t < iCol; t++) {
                    genome.addVariantFeatureLabel(cells[t]);
                }
                lineCounter++;
                //currentLine = br.readLineS();
                if (currentLine == null) {
                    return genome;
                }
            }
            //skip the head row
            if (hasHead) {
                currentLine = br.readLine();
            }
            //  bw.write(currentLine);                bw.write("\n");
            do {
                lineCounter++;
                if (needProgressionIndicator && lineCounter % 10000 == 0) {
                    String prog = String.valueOf(lineCounter);
                    System.out.print(prog);
                    char[] backSpaces = new char[prog.length()];
                    Arrays.fill(backSpaces, '\b');
                    System.out.print(backSpaces);
                }

                cells = Util.tokenizeIngoreConsec(currentLine, delimiter);
                iCol = cells.length;
                //initialize varaibles
                incomplete = true;
                currChr = null;
                makerPostionStart = -1;
                makerPostionEnd = -1;
                ref = null;
                alt = null;
                isIndel = false;
                comments = null;
                homoNum = 0;
                /*
                 Chr 	Start 	End 	Ref 	Obs 	Comments
                 1 	161003087 	161003087 	C 	T 	comments: rs1000050, a SNP in Illumina SNP arrays
                 1 	84647761 	84647761 	C 	T 	comments: rs6576700 or SNP_A-1780419, a SNP in Affymetrix SNP arrays
                 1 	13133880 	13133881 	TC 	- 	comments: rs59770105, a 2-bp deletion
                 1 	11326183 	11326183 	- 	AT 	comments: rs35561142, a 2-bp insertion
                 1 	105293754 	105293754 	A 	ATAAA 	comments: rs10552169, a block substitution
                 16 	49303427 	49303427 	C 	T 	comments: rs2066844 (R702W), a non-synonymous SNP in NOD2
                 16 	49314041 	49314041 	G 	C 	comments: rs2066845 (G908R), a non-synonymous SNP in NOD2
                 16 	49321279 	49321279 	- 	C 	comments: rs2066847 (c.3016_3017insC), a frameshift SNP in NOD2
                 13 	19661686 	19661686 	G 	- 	comments: rs1801002 (del35G), a frameshift mutation in GJB2, associated with hearing loss
                 13 	19695176 	20003944 	0 	- 	comments: a 342kb deletion encompassing GJB6, associated with hearing loss
                
                 */

                currChr = cells[indexCHROM];
                makerPostionStart = Util.parseInt(cells[indexPosStart]);
                makerPostionEnd = Util.parseInt(cells[indexPosEnd]);
                ref = cells[indexREF];
                alt = cells[indexALT];
                if (hasCommentsCol && !hasHead) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    for (int t = 5; t < iCol; t++) {
                        tmpBuffer.append(' ');
                        tmpBuffer.append(cells[t]);
                    }
                    //ignore the first space
                    comments = tmpBuffer.substring(1);
                }

                //new format on KGGSeq
                //for Indel TTCC TT--
                //for Insertion T +TTTT
                if (ref.indexOf('-') >= 0 || alt.indexOf('-') >= 0) {
                    indelNum++;
                    isIndel = true;
                }
                if (isIndel) {
                    //insertion
                    if (ref.charAt(0) == '-') {
                        alt = "+" + alt;
                        ref = "?";
                    } else {
                        sb.delete(0, sb.length());
                        makerPostionStart--;
                        //deletion
                        sb.append('?');
                        for (int t = ref.length(); t > 0; t--) {
                            sb.append('-');
                        }
                        alt = sb.toString();
                        ref = "?" + ref;
                    }
                }
                //format in KGG
                if (currChr.toLowerCase().startsWith("chr")) {
                    currChr = currChr.substring(3);
                }
                if (currChr.equals("23")) {
                    currChr = "X";
                } else if (currChr.equals("24")) {
                    currChr = "Y";
                } else if (currChr.equals("25")) {
                    currChr = "M";
                }
                //Note: without pedigree information, I do not know who are affected and unaffected
                // so this  can be changed for real case unique or control uniuqe variants
                Variant var = new Variant(makerPostionStart, ref, new String[]{alt});
                var.setIsIndel(isIndel);

                if (hasCommentsCol && hasHead) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    for (int t = 5; t < iCol; t++) {
                        var.setFeatureValue(geneFeatureStartIndex + t - 5, cells[t]);
                    }
                } else if (hasCommentsCol) {
                    var.setFeatureValue(geneFeatureStartIndex, comments);
                }

                genome.addVariant(currChr, var);
                //   bw.write(currentLine);
                //    bw.write("\n");
                acceptVarNum++;
            } while ((currentLine = br.readLine()) != null);

            br.close();
            //  bw.close();
            if (needProgressionIndicator) {
                String prog = String.valueOf(lineCounter);
                System.out.println(prog);
            }

            genome.buildVariantIndexMapOnChromosomes();
            StringBuilder message = new StringBuilder();
            genome.setVarNum(acceptVarNum);
            message.append('\n').append(lineCounter).append(" lines are scanned;  and ").append(acceptVarNum).append(" variants (").append(indelNum).append(" indels) are valid in ").append(vAFile);
            LOG.info(message);
            if (acceptVarNum == 0) {
                System.exit(1);
            }
        } catch (Exception nex) {
            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
            LOG.fatal(info, nex);

        }

        return genome;
    }

//    public Genome readCancerGenomeVariantFormat(String vAFile, boolean needProgressionIndicator) throws Exception {
//        
//        String currentLine = null;
//        String currChr = null;
//        StringBuilder tmpBuffer = new StringBuilder();
//        long lineCounter = 0;
//        
//        File dataFile = new File(vAFile);
//        if (!dataFile.exists()) {
//            throw new Exception("No ANNOVAR a file: " + dataFile.getCanonicalPath());
//        }
//        LineReader br = null;
//        if (dataFile.exists() && dataFile.getName().endsWith(".zip")) {
//            br = new CompressedFileReader(dataFile);
//        } else {
//            if (dataFile.exists() && dataFile.getName().endsWith(".gz")) {
//                br = new CompressedFileReader(dataFile);
//            } else {
//                if (dataFile.exists()) {
//                    br = new AsciiLineReader(dataFile);
//                } else {
//                    throw new Exception("No input file: " + dataFile.getCanonicalPath());
//                }
//            }
//        }
//        
//        int makerPostionStart = 0;
//        int makerPostionEnd = 0;
//        
//        String ref = null;
//        String alt = null;
//        boolean incomplete = true;
//        int acceptVarNum = 0;
//        
//        int iCol = 0;
//        //temp variables
//        int indelNum = 0;
//        boolean isIndel = false;
//        
//        StringBuilder sb = new StringBuilder();
//        String comments = null;
//        
//        Map<String, Integer> varCountMap = new HashMap<String, Integer>();
//        
//        try {
//            iCol = 0;
////skip head line
//            currentLine = br.readLineS();
//            
//            int indexCHROM = -1;
//            int indexPosStart = -1;
//            int indexPosEnd = -1;
//            int indexREF = -1;
//            int indexALT = -1;
//            int indexType = -1;
//            
//            StringTokenizer st = new StringTokenizer(currentLine.trim());
//            while (st.hasMoreTokens()) {
//                String ts = st.nextToken().trim();
//                if (ts.equals("chr")) {
//                    indexCHROM = iCol;
//                } else if (ts.equals("pos")) {
//                    indexPosStart = iCol;
//                    indexPosEnd = iCol;
//                } else if (ts.equals("ref_allele")) {
//                    indexREF = iCol;
//                } else if (ts.equals("newbase")) {
//                    indexALT = iCol;
//                } else if (ts.equals("classification")) {
//                    indexType = iCol;
//                }
//                iCol++;
//            }
//            
//            int maxColNum = indexCHROM;
//            maxColNum = Math.max(maxColNum, indexPosStart);
//            maxColNum = Math.max(maxColNum, indexPosEnd);
//            maxColNum = Math.max(maxColNum, indexREF);
//            maxColNum = Math.max(maxColNum, indexALT);
//            maxColNum = Math.max(maxColNum, indexType);
//            
//            currentLine = br.readLineS();
//            if (currentLine == null) {
//                return null;
//            }            
//            do {
//                lineCounter++;
//                if (needProgressionIndicator && lineCounter % 10000 == 0) {
//                    String prog = String.valueOf(lineCounter);
//                    System.out.print(prog);
//                    char[] backSpaces = new char[prog.length()];
//                    Arrays.fill(backSpaces, '\b');
//                    System.out.print(backSpaces);
//                }
//                st = new StringTokenizer(currentLine.trim());
//               // System.out.println(currentLine);
//                
//                //initialize varaibles
//                incomplete = true;
//                currChr = null;
//                makerPostionStart = -1;
//                makerPostionEnd = -1;
//                ref = null;
//                alt = null;
//                isIndel = false;
//                comments = null;
//
//                /*
//                 0	1	2	3	4	5	6	7	8	9	10
//                 ttype	patient	gene	classification	type	chr	pos	ref_allele	newbase	context65	cons46
//                 BRCA	BR-0001	Unknown	SNP	IGR	1	102094413	T	C	62	50
//                 BRCA	BR-0001	Unknown	SNP	IGR	1	105351609	T	C	49	50
//                 BRCA	BR-0001	Unknown	SNP	IGR	1	105803301	T	G	49	51
//                 BRCA	BR-0001	Unknown	SNP	IGR	1	106065970	G	A	34	51
//                 BRCA	BR-0001	Unknown	SNP	IGR	1	106147918	T	C	49	52
//                 BRCA	BR-0001	Unknown	SNP	IGR	1	111047687	A	G	1	55
//                 BRCA	BR-0001	Unknown	DEL	IGR	1	111460251	TTTTTGTTTTTTTG	-	64	27
//                 BRCA	BR-0001	RSBN1	SNP	Intron	1	114350413	G	A	48	39
//                 BRCA	BR-0001	DENND2C	DEL	Intron	1	115091329	A	-	16	50
//                
//                 */
//                for (iCol = 0; iCol <= maxColNum; iCol++) {
//                    if (st.hasMoreTokens()) {
//                        tmpBuffer.delete(0, tmpBuffer.length());
//                        tmpBuffer.append(st.nextToken().trim());
//                        if (iCol == indexCHROM) {
//                            currChr = tmpBuffer.toString();
//                        } else if (iCol == indexPosStart) {
//                            makerPostionStart = Util.parseInt(tmpBuffer.toString());
//                            makerPostionEnd = makerPostionStart;
//                        } else if (iCol == indexPosEnd) {
//                        } else if (iCol == indexREF) {
//                            ref = tmpBuffer.toString();
//                        } else if (iCol == indexALT) {
//                            alt = tmpBuffer.toString();
//                        } else if (iCol == indexType) {
//                            comments = tmpBuffer.toString();
//                            if (!comments.equals("SNP")) {
//                                break;
//                            }
//                        }
//                        
//                    } else {
//                        break;
//                    }
//                    if (iCol >= maxColNum) {
//                        incomplete = false;
//                    }
//                }
//                
//                if (incomplete) {
//                    continue;
//                }
//
////Ingore Indel
//                /*
//                 //new format on KGGSeq
//                 //for Indel TTCC TT--
//                 //for Insertion T +TTTT
//                 if (ref.indexOf('-') >= 0 || alt.indexOf('-') >= 0) {
//                 indelNum++;
//                 isIndel = true;
//                 }
//                 if (isIndel) {
//                 //insertion
//                 if (ref.charAt(0) == '-') {
//                 alt = "+" + alt;
//                 ref = "?";
//                 } else {
//                 sb.delete(0, sb.length());
//                 makerPostionStart--;
//                 //deletion
//                 sb.append('?');
//                 for (int t = ref.length(); t > 0; t--) {
//                 sb.append('-');
//                 }
//                 alt = sb.toString();
//                 ref = "?" + ref;
//                 }
//                 }
//                 */
//                //format in KGG
//                if (currChr.toLowerCase().startsWith("chr")) {
//                    currChr = currChr.substring(3);
//                }
//                if (currChr.equals("23")) {
//                    currChr = "X";
//                } else if (currChr.equals("24")) {
//                    currChr = "Y";
//                } else if (currChr.equals("25")) {
//                    currChr = "M";
//                }
//                if (ref.equals("-") || alt.equals("-")) {
//                    continue;
//                }
//               
//                String label = currChr + ":" + makerPostionStart + ":" + ref + ":" + alt;
//                Integer count = varCountMap.get(label);
//                if (count == null) {
//                    varCountMap.put(label, 1);
//                } else {
//                    varCountMap.put(label, count + 1);
//                }
//                acceptVarNum++;
//            } while ((currentLine = br.readLineS()) != null);
//            
//        } catch (NumberFormatException nex) {
//            String info = nex.toString() + " when parsing at line " + lineCounter + ": " + currentLine;
//            // LOG.error(nex, info);
//            throw new Exception(info);
//        }
//        br.close();
//        
//        StringBuilder message = new StringBuilder();
//        message.append('\n').append(lineCounter).append(" lines are scanned;  and ").append(acceptVarNum).append(" variants (").append(indelNum).append(" indels) are valid in ").append(vAFile);
//        LOG.info(message);
//        if (acceptVarNum == 0) {
//            System.exit(1);
//        }
//        
//        Genome genome = new Genome(dataFile.getName());
//        genome.addVariantFeatureLabel("Comments");
//        
//        for (Map.Entry<String, Integer> item : varCountMap.entrySet()) {
//            String[] varInfo = item.getKey().split(":");
//            makerPostionStart = Util.parseInt(varInfo[1]);
//            Variant var = new Variant(makerPostionStart, varInfo[2], new String[]{varInfo[3]});
//            var.addFeatureValue(item.getValue().toString());
//            genome.addVariant(varInfo[0], var);
//        }
//        
//        genome.buildVariantIndexMapOnChromosomes();
//        genome.setVarNum(acceptVarNum);
//        varCountMap.clear();
//        return genome;
//    }
//    
//    public static void main(String[] args) {
//        try {
//            
//            BufferedReader br = new BufferedReader(new FileReader("./testdata/crh.flt.annovar"));
//            BufferedWriter bw = new BufferedWriter(new FileWriter("./testdata/crh.flt.pileup.hg19"));
//
//            // BufferedReader br = new BufferedReader(new FileReader("./testdata/sca.flt.annovar"));
//            // BufferedWriter bw = new BufferedWriter(new FileWriter("./testdata/sca.flt.pileup.hg19"));
//            File CHAIN_FILE = new File("./resources/hg18ToHg19.over.chain");
//            LiftOver liftOver = new LiftOver(CHAIN_FILE);
//            
//            String line = null;
//            int failedVar = 0;
//            while ((line = br.readLine()) != null) {
//                if (line.trim().length() == 0) {
//                    continue;
//                }
//                // System.out.println(line);
//                String[] cells = line.split(" ");
//                Interval int1 = new Interval("chr" + cells[0], Util.parseInt(cells[1]), Util.parseInt(cells[1]));
//                Interval int2 = liftOver.liftOver(int1);
//                if (int2 != null) {
//                    bw.write("chr" + cells[0] + " " + int2.getStart() + " " + cells[3] + " " + cells[4] + "\n");
//                } else {
//                    failedVar++;
//                }
//            }
//            br.close();
//            bw.close();
//            System.out.println(failedVar + " convertions are failed!");
//        } catch (Exception ex) {
//            ex.printStackTrace();
//        }
//        
//    }
}
