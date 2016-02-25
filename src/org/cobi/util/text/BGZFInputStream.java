/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.text;

import cern.colt.list.IntArrayList;
import cern.colt.list.LongArrayList;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.GZIPInputStream;
import org.cobi.util.thread.Task;

public class BGZFInputStream {

    static final int BUF_SIZE = 1 * 1024 * 1024;

    File input;
    File index;
    int threadNum = 1;
    long pos[];
    RandomAccessFile raf;
    byte[] buf = new byte[BUF_SIZE];
    boolean notBGZF = false;
    int intFileFormat = -1;
    public BZPartReader[] spider;
    public BZPartBuilder[] builder;

    IntArrayList ialtIndex;
    LongArrayList laltIndex;

    public int getThreadNum() {
        return threadNum;
    }

    public BGZFInputStream(String strFile, int n_thread) {
        this.input = new File(strFile);
        this.threadNum = n_thread;
        pos = new long[n_thread + 1];
        if (strFile.endsWith(".gz")) {
            intFileFormat = 0;
        } else {
            intFileFormat = 1;
        }

    }

    public BGZFInputStream(String strFile, int n_thread, boolean isGz) {
        this.input = new File(strFile);
        this.threadNum = n_thread;
        pos = new long[n_thread + 1];
        if (isGz) {
            intFileFormat = 0;
        } else {
            intFileFormat = 1;
        }
    }

    public void adjustPos() throws IOException {
        for (int i = 0; i < threadNum; i++) {
            pos[i] = (input.length() / threadNum * i);
        }
        pos[pos.length - 1] = input.length();
        if (intFileFormat == 1) {
            //no need split the file when its size is small
            if (input.length() < 1024 * 1024) {
                threadNum = 1;
                pos[1] = pos[pos.length - 1];
            }
            return;
        }

        raf = new RandomAccessFile(input.getCanonicalFile(), "r");
        for (int i = 1; i < threadNum; i++) {
            raf.seek(pos[i]);
            boolean boolNoFound = true;
            int n_buf = -1;//if n_buf>0,we can consider this file is not bgzf format, but the gzip format. 
            do {
                raf.read(buf);
                n_buf++;
                if (n_buf > 0) {
                    notBGZF = true;
                    threadNum = 1;
                    pos[1] = pos[pos.length - 1];
                    // System.out.println("The file is gzip-format, not bgzip-format!");
                    break;
                }
                for (int id = 0; id < buf.length - 3; id++) {
                    if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4) { //This should be used unsigned number or others. 
                        pos[i] += (id + n_buf * buf.length);
                        boolNoFound = false;
                        break;
                    }
                }
            } while (boolNoFound);
            if (notBGZF) {
                break;
            }
        }
        raf.close();
        //For file with small size and many threads. 
        for (int i = 0; i < (threadNum - 1); i++) {
            if (pos[i] == pos[i + 1]) {
                threadNum = 1;
                pos[1] = pos[pos.length - 1];
            }
        }
    }

    public void adjustPos(long start, long end) throws IOException {
        long dist = end - start;
        for (int i = 0; i < threadNum; i++) {
            pos[i] = start + (dist / threadNum * i);
        }
        pos[pos.length - 1] = end;

        raf = new RandomAccessFile(input.getCanonicalFile(), "r");
        for (int i = 1; i < threadNum; i++) {
            raf.seek(pos[i]);
            boolean boolNoFound = true;
            int n_buf = -1;//if n_buf>0,we can consider this file is not bgzf format, but the gzip format. 
            do {
                raf.read(buf);
                n_buf++;
                if (n_buf > 0) {
                    notBGZF = true;
                    threadNum = 1;
                    pos[1] = pos[pos.length - 1];
                    // System.out.println("The file is gzip-format, not bgzip-format!");
                    break;
                }
                for (int id = 0; id < buf.length - 3; id++) {
                    if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4) { //This should be used unsigned number or others. 
                        pos[i] += (id + n_buf * buf.length);
                        boolNoFound = false;
                        break;
                    }
                }
            } while (boolNoFound);
            if (notBGZF) {
                break;
            }
        }
        raf.close();
        //For file with small size and many threads. 
        for (int i = 0; i < (threadNum - 1); i++) {
            if (pos[i] == pos[i + 1]) {
                threadNum = 1;
                pos[1] = pos[pos.length - 1];
            }
        }
        if (threadNum == 1) {
            // System.out.println("Notice: Only one thread is created because of the short of file(region).");
        }
    }

    public void creatSpider() throws IOException {
        spider = new BZPartReader[this.threadNum];
        for (int i = 0; i < this.threadNum; i++) {
            spider[i] = new BZPartReader(i, intFileFormat, this.pos[i], this.pos[i + 1], '\t');
            // System.out.println("BZPartReader" + i + ": created!");
        }
    }

    public void creatSpider(boolean boolRFL) throws IOException {
        spider = new BZPartReader[this.threadNum];
        for (int i = 0; i < this.threadNum; i++) {
            spider[i] = new BZPartReader(i, intFileFormat, this.pos[i], this.pos[i + 1], '\t', boolRFL);
            // System.out.println("BZPartReader" + i + ": created!");
        }
    }

    public void creatSpider(long longStart, long longEnd) throws IOException {
        spider = new BZPartReader[1];
        spider[0] = new BZPartReader(0, 0, longStart, longEnd, '\t');
    }

    public void creatSpider(long longPos[]) throws IOException {
        spider = new BZPartReader[longPos.length - 1];
        for (int i = 0; i < spider.length; i++) {
            spider[i] = new BZPartReader(i, 0, this.pos[i], this.pos[i + 1], '\t');
            // System.out.println("BZPartReader" + i + ": created!");
        }
    }

    private void creatBuilders(int intChrom, int intStart, boolean hasHead) throws IOException {
        builder = new BZPartBuilder[pos.length - 1];
        for (int i = 0; i < (pos.length - 1); i++) {
            builder[i] = new BZPartBuilder(i, this.pos[i], this.pos[i + 1], intChrom, intStart, hasHead);
        }
    }

    public class BZPartReader {

        int spiderID;
        int intFormat;
        long longRemainSize = -1;
        long longEnd = -1;
        InputStream inputStream;
        LimitInputStream cis;

        int intRead = -1;
        byte[] bytBuffer = new byte[BUF_SIZE];
        int intLineStart = 0;
        int intCurrPos = 0;
        int intTempBufferLengthACC = 0;
        byte[] bytTempBuffer = new byte[BUF_SIZE];
        byte[] bytStartLine = null;
        byte[] bytEndLine = null;
        byte bytDelimiter = 9;
        int[] intSep;

        public BZPartReader(int spiderID, int intFormat, long longStart, long longEnd, char chrDelimiter) throws IOException {
            longRemainSize = (longEnd - longStart);
            this.spiderID = spiderID;
            this.intFormat = intFormat;
            //Plan one           
//            RandomAccessFile raf=new RandomAccessFile(fleInput.getCanonicalFile(), "r");
//            raf.seek(longStart);
//            input=new GZIPInputStream(Channels.newInputStream(raf.getChannel()),BUF_SIZE);
//            this.bytDelimiter=(byte)chrDelimiter;
//            if(this.spiderID!=0)    bytStartLine=this.readLine(new int[3000]);

            //Plan two. 
//            RandomAccessFile raf=new RandomAccessFile(fleInput.getCanonicalFile(), "r");
//            MappedByteBuffer mbb = raf.getChannel().map(FileChannel.MapMode.READ_ONLY, longStart,longRemainSize); 
//            for(int i=0;i<300;i++)   System.out.print(mbb.get());
//            System.out.println();
            //Plan Three. 
//            InputStream is = new BufferedInputStream(new FileInputStream(input));
//            LimitInputStream cis = new LimitInputStream(is, longEnd);
//            cis.skip(longStart);
//            ReadableByteChannel inChannel = Channels.newChannel(cis);
//            input= new GZIPInputStream(Channels.newInputStream(inChannel));
//            input=new GZIPInputStream(cis);
            if (intFormat == 0) {
                InputStream is = new BufferedInputStream(new FileInputStream(input));
                cis = new LimitInputStream(is, longEnd);
                cis.skip(longStart);
                inputStream = new GZIPInputStream(cis);
            } else {
                InputStream is = new BufferedInputStream(new FileInputStream(input));
                cis = new LimitInputStream(is, longEnd);
                cis.skip(longStart);
                inputStream = cis;
            }

            if (this.spiderID != 0) {
                bytStartLine = this.readLine();
                spider[this.spiderID - 1].bytEndLine = bytStartLine;
                // System.out.println(new String(bytStartLine));
            }
            bytDelimiter = (byte) chrDelimiter;
        }

        public BZPartReader(int spiderID, int intFormat, long longStart, long longEnd, char chrDelimiter, boolean boolRFL) throws IOException {
            longRemainSize = (longEnd - longStart);
            this.spiderID = spiderID;
            this.intFormat = intFormat;

            if (intFormat == 0) {
                InputStream is = new BufferedInputStream(new FileInputStream(input));
                cis = new LimitInputStream(is, longEnd);
                cis.skip(longStart);
                inputStream = new GZIPInputStream(cis);
            } else {
                InputStream is = new BufferedInputStream(new FileInputStream(input));
                cis = new LimitInputStream(is, longEnd);
                cis.skip(longStart);
                inputStream = cis;
            }

            if (this.spiderID == 0 && boolRFL) {//boolRFL<=>longStart!=0
                bytStartLine = this.readLine();
            }

            if (this.spiderID != 0) {
                bytStartLine = this.readLine();
                spider[this.spiderID - 1].bytEndLine = bytStartLine;
            }
            bytDelimiter = (byte) chrDelimiter;
        }

        public BZPartReader(int spiderID, long longStart, long longEnd, int[] intSep, boolean hasHead) throws IOException {
            this.longEnd = longEnd;
            longRemainSize = (longEnd - longStart);
            this.spiderID = spiderID;
            this.intSep = intSep;

            InputStream is = new BufferedInputStream(new FileInputStream(input));
            cis = new LimitInputStream(is, longEnd);
            cis.skip(longStart);
            cis.mark((int) longRemainSize);
            if (intFormat == 0) {
                inputStream = new GZIPInputStream(cis);
            } else {
                inputStream = cis;
            }

            if (this.spiderID != 0) {
                this.readLine(intSep);
                bytStartLine = readLine(intSep);

            } else {
                if (hasHead) {
                    this.readLine(intSep);
                }
                bytStartLine = readLine(intSep);
            }
            cis.reset();
        }

        public void closeInputStream() throws IOException {
            inputStream.close();
        }

        public synchronized byte[] readLine(int[] intDelPos) throws IOException {
            int intDelPosMarker = 0;
            intDelPos[intDelPosMarker] = 0;
            int len = intDelPos.length;
            do {
                if (intRead == -1) {
                    intRead = (inputStream.read(bytBuffer));
                    if (intRead == -1) {
                        //The end of the block is not a complete line. 
                        if (intTempBufferLengthACC != 0) {
                            bytBuffer = new byte[intTempBufferLengthACC + (bytEndLine == null ? 0 : bytEndLine.length) + 1];
                            System.arraycopy(bytTempBuffer, 0, bytBuffer, 0, intTempBufferLengthACC);
                            if (bytEndLine != null) {
                                System.arraycopy(this.bytEndLine, 0, bytBuffer, intTempBufferLengthACC, this.bytEndLine.length);
                            }
                            bytBuffer[this.bytBuffer.length - 1] = (byte) '\n';
                            intTempBufferLengthACC = 0;
                            intDelPosMarker = 0;
                        } else {
                            return null;
                        }
                    }
                }

                intLineStart = intCurrPos;
                while (intCurrPos != intRead) {
                    if (bytBuffer[intCurrPos] == bytDelimiter) {
                        ++intDelPosMarker;
                        if (intDelPosMarker < len) {
                            intDelPos[intDelPosMarker] = intCurrPos - intLineStart + intTempBufferLengthACC;
                        }
                    } else if (bytBuffer[intCurrPos] == 10) {
                        //parse the line. 
                        int intLineLength = intCurrPos - intLineStart;//don't contaion \n
                        byte[] bytLine = null;
                        if (intTempBufferLengthACC != 0) {
                            // bytLine = new byte[intTempBufferLengthACC + intLineLength];
                            //System.arraycopy(bytTempBuffer, 0, bytLine, 0, intTempBufferLengthACC);                        
                            bytLine = Arrays.copyOfRange(bytTempBuffer, 0, intTempBufferLengthACC + intLineLength);
                            System.arraycopy(bytBuffer, intLineStart, bytLine, intTempBufferLengthACC, intLineLength);
                            intTempBufferLengthACC = 0;
                        } else {
                            // bytLine = new byte[intLineLength];
                            // System.arraycopy(bytBuffer, intLineStart, bytLine, 0, intLineLength);
                            bytLine = Arrays.copyOfRange(bytBuffer, intLineStart, intCurrPos);
                        }
                        intCurrPos++;
                        if (intCurrPos == intRead) {
                            intRead = -1;
                            intCurrPos = 0;
                        }

                        ++intDelPosMarker;
                        if (intDelPosMarker < len) {
                            //return 13 new line 10
                            //the new line or return line symbol is not included
                            if (intCurrPos > 1 && bytBuffer[intCurrPos - 2] == 13) {
                                intDelPos[intDelPosMarker] = bytLine.length - 1;
                            } else {
                                intDelPos[intDelPosMarker] = bytLine.length;
                            }
                            intDelPos[0] = intDelPosMarker;
                        } else {
                            intDelPos[0] = len;
                        }

                        return bytLine;
                    }
                    intCurrPos++;
                }

                //The buffer ends with imcomplete line. 
                int intTempBufferLength = intCurrPos - intLineStart;
                if ((bytTempBuffer.length - intTempBufferLengthACC) < intTempBufferLength) {
                    bytTempBuffer = Arrays.copyOf(bytTempBuffer, bytTempBuffer.length + intTempBufferLength * 2);
                }
                System.arraycopy(bytBuffer, intLineStart, bytTempBuffer, intTempBufferLengthACC, intTempBufferLength);
                intTempBufferLengthACC += intTempBufferLength;
                intRead = -1;
                intCurrPos = 0;
            } while (intTempBufferLengthACC != 0);

            return null;
        }

        public synchronized byte[] readLine() throws IOException {
            do {
                if (intRead == -1) {
                    intRead = (inputStream.read(bytBuffer));
                    if (intRead == -1) {
                        //The end of the block is not a complete line. 
                        if (intTempBufferLengthACC != 0) {
                            bytBuffer = new byte[intTempBufferLengthACC + (bytEndLine == null ? 0 : bytEndLine.length) + 1];
                            System.arraycopy(bytTempBuffer, 0, bytBuffer, 0, intTempBufferLengthACC);
                            if (bytEndLine != null) {
                                System.arraycopy(this.bytEndLine, 0, bytBuffer, intTempBufferLengthACC, this.bytEndLine.length);
                            }
                            bytBuffer[this.bytBuffer.length - 1] = (byte) '\n';
                            intTempBufferLengthACC = 0;
                        } else {
                            return null;
                        }
                    }
                }

                intLineStart = intCurrPos;
                while (intCurrPos != intRead) {
                    if (bytBuffer[intCurrPos] == 10) {
                        //parse the line. 
                        int intLineLength = intCurrPos - intLineStart;//don't contaion \n
                        byte[] bytLine = null;
                        if (intTempBufferLengthACC != 0) {
                            // bytLine = new byte[intTempBufferLengthACC + intLineLength];
                            //System.arraycopy(bytTempBuffer, 0, bytLine, 0, intTempBufferLengthACC);                        
                            bytLine = Arrays.copyOfRange(bytTempBuffer, 0, intTempBufferLengthACC + intLineLength);
                            System.arraycopy(bytBuffer, intLineStart, bytLine, intTempBufferLengthACC, intLineLength);
                            intTempBufferLengthACC = 0;
                        } else {
                            // bytLine = new byte[intLineLength];
                            // System.arraycopy(bytBuffer, intLineStart, bytLine, 0, intLineLength);
                            bytLine = Arrays.copyOfRange(bytBuffer, intLineStart, intCurrPos);
                        }
                        intCurrPos++;
                        if (intCurrPos == intRead) {
                            intRead = -1;
                            intCurrPos = 0;
                        }
                        return bytLine;
                    }
                    intCurrPos++;
                }

                //The buffer ends with imcomplete line. 
                int intTempBufferLength = intCurrPos - intLineStart;
                if ((bytTempBuffer.length - intTempBufferLengthACC) < intTempBufferLength) {
                    bytTempBuffer = Arrays.copyOf(bytTempBuffer, bytTempBuffer.length + intTempBufferLength * 2);
                }
                System.arraycopy(bytBuffer, intLineStart, bytTempBuffer, intTempBufferLengthACC, intTempBufferLength);
                intTempBufferLengthACC += intTempBufferLength;
                intRead = -1;
                intCurrPos = 0;
            } while (intTempBufferLengthACC != 0);

            return null;
        }

        public byte[] getStartLine() {
            return bytStartLine;
        }

        public long getNextBlock() throws IOException {
            long longNext = 0;
            //cis.skip(longCurr + 1);
            boolean boolNoFound = true;
            int n_buf = -1;
            int intContent = 0;
            do {
                intContent = cis.read(buf);
                if (intContent == -1) {
                    return -1;
                }
                n_buf++;
                for (int id = 0; id < intContent - 3; id++) {
                    if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4) { //This should be used unsigned number or others. 
                        longNext += (id + n_buf * buf.length + 1);
                        boolNoFound = false;
                        break;
                    }
                }
            } while (boolNoFound);
            return longNext;
        }

        public boolean jumpTo(long longDict, int[] intSepTmp) throws IOException {
            longRemainSize -= longDict;
            // System.out.println(longRemainSize);
            if (longRemainSize <= 0) {
                bytStartLine = null;
                return false;
            }
            while (longDict > 0) {
                longDict -= cis.skip(longDict);
            }
            cis.mark((int) longRemainSize);
//            cis.skip(longDict);

            if (intFormat == 0) {
                inputStream = new GZIPInputStream(cis);
            } else {
                inputStream = cis;
            }

            intRead = -1;
            intTempBufferLengthACC = 0;
            intCurrPos = 0;
            readLine(intSepTmp);
            bytStartLine = readLine(intSepTmp);
            // System.out.println(new String(bytStartLine));
            cis.reset();
//            inputStream.reset();
            return true;
        }

        public boolean prepareRead(long longDict) throws IOException {
            longRemainSize -= longDict;
            // System.out.println(longRemainSize);
            if (longRemainSize <= 28) {
                bytStartLine = null;
                return false;
            }
            while (longDict > 0) {
                longDict -= cis.skip(longDict);
//                longDict-=inputStream.skip(longDict);
            }
            cis.mark((int) longRemainSize);
            cis.skip(longDict);
//            inputStream.mark((int) longRemainSize);

            if (intFormat == 0) {
                // inputStream = new GZIPInputStream(cis);
            } else {
                //inputStream = cis;
            }

            intRead = -1;
            intTempBufferLengthACC = 0;
            intCurrPos = 0;

            return true;
        }

    }

    public class BZPartBuilder extends Task implements Callable<String> {

        long longCurr;
        long longEnd;
        long longRemainSize;
        int spiderID;
        int intChrom;
        int intStart;
        RandomAccessFile raf = null;
        ArrayList<String> altIndex = null;
        BZPartReader bzpr = null;
        int[] intSep;
        byte[] buf = new byte[BUF_SIZE];

        public BZPartBuilder(int spiderID, long longStart, long longEnd, int intChrom, int intStart, boolean hasHead) throws FileNotFoundException, IOException {
            this.spiderID = spiderID;
            this.longCurr = longStart;
            this.longEnd = longEnd;
            this.longRemainSize = (longEnd - longStart);
            this.intChrom = intChrom;
            this.intStart = intStart;
            raf = new RandomAccessFile(input.getCanonicalPath(), "r");
            raf.seek(longStart);
            intSep = new int[Math.max(intStart + 1, intChrom + 1) + 1];
            bzpr = new BZPartReader(spiderID, longStart, longEnd, intSep, hasHead);
            altIndex = new ArrayList<String>();
        }

        @Override
        public String call() {
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
                        altIndex.add(sb.toString());
                        sb.delete(0, sb.length());
                    } else {
                        break;
                    }
                    longNext = getNextStart(raf, longCurr);
                    if (longNext == -1 || longRemainSize <= 0) {
                        break;
                    }
                } while (bzpr.jumpTo(longNext, intSep));
                raf.close();
                bzpr.closeInputStream();
                return "TD" + this.spiderID + " is finished!";
            } catch (Exception e) {
                return "TD" + this.spiderID + " is wrong!";
            }
        }

        public long getNextStart(RandomAccessFile raf, long longCurr) throws IOException {
            long longNext = 0;
            raf.seek(longCurr + 1);
            boolean boolNoFound = true;
            int n_buf = -1;
            int intContent = 0;
            do {
                intContent = raf.read(buf);
                if (intContent == -1) {
                    longRemainSize = 0;
                    return -1;
                }
                n_buf++;
                for (int id = 0; id < intContent - 3; id++) {
                    if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4) { //This is not a accurate method to identify the head. 
                        longNext += (id + n_buf * buf.length + 1);
                        boolNoFound = false;
                        break;
                    }
                }
            } while (boolNoFound);
            longRemainSize -= longNext;
            return longNext;
        }

    }

    public class ParseTask extends Task implements Callable<String> {

        BZPartReader bzprSpider;
        int intChrom;
        int intStart;
        long longCurr;

        public ParseTask(BZPartReader bzprSpider, int intChrom, int intStart, long longCurr) {
            this.bzprSpider = bzprSpider;
            this.intChrom = intChrom;
            this.intStart = intStart;
        }

        List<String[]> indexes = new ArrayList<String[]>();

        @Override
        public String call() throws Exception {
            byte[] currentLine = null;
            String tmpStr;
            int[] intSep = new int[Math.max(intStart + 1, intChrom + 1) + 1];
            long longNext = 0;
            bzprSpider.prepareRead(longNext);
            // while (bzprSpider.jumpTo(longNext, intSep)) {
            while ((currentLine = bzprSpider.readLine(intSep)) != null) {
                // String strTemp = new String(currentLine);
                if (currentLine != null) {
                    String[] poss = new String[3];
                    int intA = intChrom == 0 ? 0 : (intSep[intChrom] + 1);
                    int intLen = intSep[intChrom + 1] - intA;
                    poss[0] = new String(currentLine, intA, intLen);

                    intA = intStart == 0 ? 0 : (intSep[intStart] + 1);
                    intLen = intSep[intStart + 1] - intA;
                    poss[1] = new String(currentLine, intA, intLen);
                    poss[2] = String.valueOf(longCurr);
                    indexes.add(poss);
                    System.out.println(longCurr);
                } else {
                    break;
                }
                longNext = bzprSpider.getNextBlock();
                longCurr += longNext;
            }

            bzprSpider.closeInputStream();
            return "";
        }
    }

    public int buildIndex(String strFile, int intChrom, int intStart, boolean hasHead) throws IOException, InterruptedException, ExecutionException {
        File fleInput = new File(strFile);
        if (!fleInput.exists()) {
            return 1;
        }
        File fleOutput = new File(strFile + ".idx");
//        BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(fleOutput));      
        BufferedWriter bw = new BufferedWriter(new FileWriter(fleOutput));
        this.creatBuilders(intChrom, intStart, hasHead);

        ExecutorService exec = Executors.newFixedThreadPool(builder.length);
        final CompletionService<String> serv = new ExecutorCompletionService<String>(exec);
        int runningThread = 0;
        for (int s = 0; s < builder.length; s++) {
            serv.submit(builder[s]);
            runningThread++;
        }

        for (int s = 0; s < runningThread; s++) {
            Future<String> task = serv.take();
            String infor = task.get();
            // System.out.println(infor);
        }
        exec.shutdown();

        for (BZPartBuilder builder1 : builder) {
            for (int j = 0; j < builder1.altIndex.size(); j++) {
                bw.write(builder1.altIndex.get(j));
            }
//            altIndex.addAll(builder[i].altIndex);
        }
        bw.close();
        index = fleOutput;
        return 0;
    }

    public long getNextStart(RandomAccessFile raf, long longCurr) throws IOException {
        long longNext = 0;
        raf.seek(longCurr + 1);
        boolean boolNoFound = true;
        int n_buf = -1;
        int intContent = 0;
        do {
            intContent = raf.read(buf);
            if (intContent == -1) {
                return -1;
            }
            n_buf++;
            for (int id = 0; id < intContent - 3; id++) {
                if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4) { //This should be used unsigned number or others. 
                    longNext += (id + n_buf * buf.length + 1);
                    boolNoFound = false;
                    break;
                }
            }
        } while (boolNoFound);
        return longNext;
    }

    public void readIndex(boolean hasHead, String chrom) throws FileNotFoundException, IOException {
        if (!index.exists()) {
            return;
        }
        BufferedReader br = new BufferedReader(new FileReader(index));
        String strLine;
        if (hasHead) {
            strLine = br.readLine();
        }
        boolean hasChromName = true;
        if (chrom == null || chrom.isEmpty()) {
            hasChromName = false;
        }
        boolean boolFlag = false;
        ialtIndex = new IntArrayList();
        laltIndex = new LongArrayList();
        //It is very often that two or more chromsomes 
        String lastChromPos = null;
        boolean hasNotStart = true;
        if (hasChromName) {
            while ((strLine = br.readLine()) != null) {
                String[] strItem = strLine.split("\t");
                if (!strItem[0].equals(chrom)) {
                    if (boolFlag) {
                        break;
                    }
                    lastChromPos = strItem[2];
                    continue;
                } else {
                    if (hasNotStart) {
                        if (lastChromPos != null) {
                            ialtIndex.add(0);
                            laltIndex.add(Long.valueOf(lastChromPos));
                        }
                        hasNotStart = false;
                    }
                    ialtIndex.add(Integer.valueOf(strItem[1]));
                    laltIndex.add(Long.valueOf(strItem[2]));
                    boolFlag = true;
                }

            }
        } else {
            while ((strLine = br.readLine()) != null) {
                String[] strItem = strLine.split("\t");
                ialtIndex.add(Integer.valueOf(strItem[0]));
                laltIndex.add(Long.valueOf(strItem[1]));
                boolFlag = true;
            }
        }
        br.close();
    }

    public long[] findIndex(int start, int end) {
        if (ialtIndex == null | ialtIndex.isEmpty()) {
            return null;
        }
        if (laltIndex == null | laltIndex.isEmpty()) {
            return null;
        }
        if (start > end) {
            int temp = start;
            start = end;
            end = temp;
        }

        long[] secPos = new long[2];
        int posStart = ialtIndex.binarySearch(start);
        if (posStart >= 0) {
            secPos[0] = laltIndex.get(posStart);
        } else {
            posStart = -posStart - 2;
            if (posStart < 0) {
                secPos[0] = laltIndex.get(0);
            } else {
                secPos[0] = laltIndex.get(posStart);
            }
        }

        int posEnd = ialtIndex.binarySearch(end);
        if (posEnd >= 0) {
            posEnd = posEnd + 1;
            if (posEnd >= laltIndex.size()) {
                secPos[1] = input.length();
            } else {
                secPos[1] = laltIndex.get(posEnd);
            }
        } else {
            posEnd = -posEnd;
            if (posEnd >= laltIndex.size()) {
                secPos[1] = input.length();
            } else {
                secPos[1] = laltIndex.get(posEnd);
            }
        }
        return secPos;
    }

    public boolean checkIndex() throws IOException {
        if (!input.exists()) {
            return false;
        }
        File fle = new File(input.getCanonicalPath() + ".idx");
        if (!fle.exists()) {
            return false;
        } else {
            index = fle;
            return true;
        }
    }

    public static void main(String[] args) throws IOException {
        String strFile = "D:\\kgg3d5\\resources\\LD\\v5\\all\\ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";

        // String strFile = "AF.all.snp.bed.bgz";
        // String strFile="D:\\01WORK\\KGGseq\\testdata\\hg19_funcnote_encode_megamix.bed.gz.Histone.cmp.gz";
        try {
//            InputStream is = new BufferedInputStream(new FileInputStream(strFile));
//            LimitInputStream cis;
//            cis = new LimitInputStream(is, 1411);
//            cis.skip(0);
            int test1 = 2;
            if (test1 == 1) {
                BGZFInputStream bf = new BGZFInputStream(strFile, 4, true);
                if (!bf.checkIndex()) {
                    bf.adjustPos();
                    bf.buildIndex(strFile, 0, 1, false);
                }

                bf.readIndex(true, "X");
                long[] pos = bf.findIndex(20922786, 22192439);
                bf.adjustPos(pos[0], pos[1]);
                bf.creatSpider(pos[0] != 0);
                int[] temp = new int[500];
                byte[] bty;
                for (int i = 0; i < bf.threadNum; i++) {
                    while ((bty = bf.spider[i].readLine(temp)) != null) {
                        String str = new String(bty);
                        System.out.println(temp[0] + " : " + str);
                        // System.out.println(str.substring(0, 20));
                    }
                }
            } else {
                BGZFInputStream bf = new BGZFInputStream(strFile, 1, true);
                bf.adjustPos();
                bf.creatSpider();
                int[] temp = new int[500];
                byte[] bty;
                int t=0;
                while ((bty = bf.spider[0].readLine(temp)) != null) {
                   // String str = new String(bty);
                   // System.out.println(temp[0] + " : " + str);
                    // System.out.println(str.substring(0, 20));
                    t++;
                }
                System.out.println(t);
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }

//        BGZFInputStream bf = new BGZFInputStream(strFile, 1);
//        int[] intDelPos = new int[3000];
//        bf.adjustPos();
//        bf.creatSpider();
//
//        int intLine = 0;
//        for (int i = 0; i < bf.threadNum; i++) {
//            byte[] bytTemp = bf.spider[i].readLine(intDelPos);
//            while (bytTemp != null) {
////            System.out.println(intLine+++": "+new String(bytTemp));               
//                bytTemp = bf.spider[i].readLine(intDelPos);
//                intLine++;
//            }
//        }
//        System.out.println(intLine);
    }
}
