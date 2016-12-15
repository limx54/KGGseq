/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import static org.cobi.util.text.BGZFInputStream.BUF_SIZE;
import org.cobi.util.thread.Task;

/**
 *
 * @author mxli
 */
public class BZPartBuilder extends Task implements Callable<String> {

    long longCurr;
    long longEnd;
    long longRemainSize;
    int spiderID;
    int intChrom;
    int intStart;
    RandomAccessFile raf = null;
    List<String> altIndex = null;
    BZPartReader bzpr = null;
    int[] intSep;
    byte[] buf = new byte[BUF_SIZE];
    String prefixName;
    int intBufPos = 0;
    int n_buf = -1;
    int intContent = 0;
    boolean loadBuf = true;
    boolean boolFlag = false;

    public BZPartBuilder(File input, int spiderID, long longStart, long longEnd, int intChrom, int intStart, boolean hasHead, String prefixName) throws FileNotFoundException, IOException {
        this.spiderID = spiderID;
        this.longCurr = longStart;
        this.longEnd = longEnd;
        this.longRemainSize = (longEnd - longStart);
        this.intChrom = intChrom;
        this.intStart = intStart;
        raf = new RandomAccessFile(input.getCanonicalPath(), "r");
        raf.seek(longStart);
        intSep = new int[Math.max(intStart + 1, intChrom + 1) + 1];
        bzpr = new BZPartReader(input, spiderID, longStart, longEnd, intSep, hasHead);
        altIndex = new ArrayList<String>();
        this.prefixName = prefixName;
    }

    @Override
    public String call() {
        try {
            long longNext = 0;
            byte[] currentLine;
            String[] poss = new String[2];
            StringBuilder sb = new StringBuilder();
            int intA, intLen;
            BufferedWriter bw = new BufferedWriter(new FileWriter(prefixName + "." + spiderID));
//                System.out.println("CP2."+spiderID+" Start1 ----->"+new Date());
            do {
                longCurr += longNext;
                currentLine = bzpr.getStartLine();
                //   String strTemp = new String(currentLine);
                if (currentLine != null) {
                    // String[] strItem = new String(currentLine).split("\t");
                    // String strIndex=strItem[intChrom]+"\t"+strItem[intStart]+"\t"+String.valueOf(longCurr)+"\tTD"+this.spiderID+"\n";
                    // String strIndex = strItem[intChrom] + "\t" + strItem[intStart] + "\t" + String.valueOf(longCurr) + "\n";
                    // System.out.print(strIndex);
                    if (intChrom >= 0) {
                        intA = intChrom == 0 ? 0 : (intSep[intChrom] + 1);
                        intLen = intSep[intChrom + 1] - intA;
                        poss[0] = new String(currentLine, intA, intLen);
                        sb.append(poss[0]);
                        sb.append('\t');
                    }

                    intA = intStart == 0 ? 0 : (intSep[intStart] + 1);
                    intLen = intSep[intStart + 1] - intA;
                    poss[1] = new String(currentLine, intA, intLen);
                    sb.append(poss[1]);
                    sb.append('\t');
                    sb.append(longCurr);
                    //altIndex.add(sb.toString());
//                        System.out.println("CP2."+spiderID+" Start2 ----->"+new Date());
                    bw.write(sb.toString());
//                        System.out.println(sb.toString());
                    bw.write("\n");
                    sb.delete(0, sb.length());
//                        System.out.println("CP2."+spiderID+" Start3 ----->"+new Date());
                } else {
                    break;
                }
//                    System.out.println("CP2."+spiderID+" Start4 ----->"+new Date());
                longNext = getNextStart(raf, longCurr);
//                    System.out.println("CP2."+spiderID+" Start5 ----->"+new Date());
                if (longNext == -1 || longRemainSize <= 0) {
                    break;
                }
//         System.out.println(longRemainSize);
            } while (bzpr.jumpTo(longNext, intSep));
            raf.close();
            bzpr.closeInputStream();
            bw.close();
//                System.out.println("TD" + this.spiderID + " is finished!");
//                System.out.println("CP2."+spiderID+"End ----->"+new Date());
            return "TD" + this.spiderID + " is finished!";
        } catch (Exception e) {
            e.printStackTrace();
            return "TD" + this.spiderID + " is wrong!";
        }
    }

    public void directWrite(BufferedWriter bw) {
        try {
            long longNext = 0;
            byte[] currentLine;
            String[] poss = new String[2];
            StringBuilder sb = new StringBuilder();
            int intA, intLen;
            do {
                longCurr += longNext;
                currentLine = bzpr.getStartLine();
                //                String strTemp = new String(currentLine);
                if (currentLine != null) {
                    // String[] strItem = new String(currentLine).split("\t");
                    // String strIndex=strItem[intChrom]+"\t"+strItem[intStart]+"\t"+String.valueOf(longCurr)+"\tTD"+this.spiderID+"\n";
                    // String strIndex = strItem[intChrom] + "\t" + strItem[intStart] + "\t" + String.valueOf(longCurr) + "\n";
                    // System.out.print(strIndex);
                    if (intChrom >= 0) {
                        intA = intChrom == 0 ? 0 : (intSep[intChrom] + 1);
                        intLen = intSep[intChrom + 1] - intA;
                        poss[0] = new String(currentLine, intA, intLen);
                        sb.append(poss[0]);
                        sb.append('\t');
                    }

                    intA = intStart == 0 ? 0 : (intSep[intStart] + 1);
                    intLen = intSep[intStart + 1] - intA;
                    poss[1] = new String(currentLine, intA, intLen);
                    sb.append(poss[1]);
                    sb.append('\t');
                    sb.append(longCurr);
                    sb.append('\n');
                    bw.write(sb.toString());
                    sb.delete(0, sb.length());
                } else {
                    break;
                }
                longNext = getNextStart(raf, longCurr);
                if (longNext == -1 || longRemainSize <= 0) {
                    break;
                }
//                    System.out.println(longRemainSize);
            } while (bzpr.jumpTo(longNext, intSep));
            raf.close();
            bzpr.closeInputStream();
            System.out.println("TD" + this.spiderID + " is finished!");

        } catch (Exception e) {
            e.printStackTrace();
            //return "TD" + this.spiderID + " is wrong!";
        }
    }

    private long getNextStart(RandomAccessFile raf, long longCurr) throws IOException {
        long longNext = 0;
//            raf.seek(longCurr + 1);
        boolean boolNoFound = true;
        do {
            if (loadBuf) {
                raf.seek(longCurr + longNext);
                intContent = raf.read(buf);//when a head was recognized, it will be read in once. So it's time-consuming, especiall when file is large. 
                intBufPos = 0;
                if (intContent == -1) {
                    longRemainSize = 0;
                    return -1;
                }
//                    if(boolFlag){
//                        n_buf=-1;
//                        boolFlag=false;
//                    }
//                    n_buf++;
                loadBuf = false;
            }

            /*
                for (int id = 0; id < intContent - 3; id++) {
                    if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4&& buf[id + 4] == 0) { //This is not a accurate method to identify the head. 
                        longNext += (id + n_buf * buf.length + 1);
                        boolNoFound = false;
                        break;
                    }
                }
             */
            for (int id = intBufPos + 1; id < intContent - 4; id++) {
                if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4 && buf[id + 4] == 0) { //This is not a accurate method to identify the head. 
                    longNext += (id - intBufPos);
                    intBufPos = id;
                    boolNoFound = false;
                    boolFlag = true;
                    break;
                }
            }
            if (boolNoFound) {
                longNext += intContent - 4 - intBufPos;
                intBufPos = intContent - 4;
                loadBuf = true;
            }
        } while (boolNoFound);
        longRemainSize -= longNext;
        return longNext;
    }
}
