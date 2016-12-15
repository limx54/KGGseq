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
<<<<<<< HEAD
import java.io.FileOutputStream;
=======
>>>>>>> origin/master
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
<<<<<<< HEAD
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
=======
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
>>>>>>> origin/master
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

<<<<<<< HEAD
    File inputFile;
=======
    File input;
>>>>>>> origin/master
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
        this.inputFile = new File(strFile);
        this.threadNum = n_thread;
        pos = new long[n_thread + 1];
        if (strFile.endsWith(".gz")) {
            intFileFormat = 0;
        } else {
            intFileFormat = 1;
        }

    }

    public BGZFInputStream(String strFile, int n_thread, boolean isGz) {
<<<<<<< HEAD
        this.inputFile = new File(strFile);
=======
        this.input = new File(strFile);
>>>>>>> origin/master
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
            pos[i] = (inputFile.length() / threadNum * i);
        }
        pos[pos.length - 1] = inputFile.length();
        if (intFileFormat == 1) {
            //no need split the file when its size is small
<<<<<<< HEAD
            if (inputFile.length() < 1024 * 1024) {
=======
            if (input.length() < 1024 * 1024) {
>>>>>>> origin/master
                threadNum = 1;
                pos[1] = pos[pos.length - 1];
            }
            return;
        }

        raf = new RandomAccessFile(inputFile.getCanonicalFile(), "r");
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
<<<<<<< HEAD
                    break;
                }
                for (int id = 0; id < buf.length - 4; id++) {
                    if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4 && buf[id + 4] == 0) { //This should be used unsigned number or others. 
=======
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
>>>>>>> origin/master
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

    public void adjustPos(long start, long end) throws IOException {
        long dist = end - start;
        for (int i = 0; i < threadNum; i++) {
            pos[i] = start + (dist / threadNum * i);
        }
        pos[pos.length - 1] = end;

        raf = new RandomAccessFile(inputFile.getCanonicalFile(), "r");
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
                for (int id = 0; id < buf.length - 4; id++) {
                    if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4 && buf[id + 4] == 0) { //This should be used unsigned number or others. 
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
<<<<<<< HEAD
            spider[i] = new BZPartReader(inputFile, i, intFileFormat, this.pos[i], this.pos[i + 1], '\t');
            // System.out.println("BZPartReader" + i + ": created!");
            if (i > 0) {
                spider[i - 1].bytEndLine = spider[i].getStartLine();
            }
=======
            spider[i] = new BZPartReader(i, intFileFormat, this.pos[i], this.pos[i + 1], '\t');
            // System.out.println("BZPartReader" + i + ": created!");
>>>>>>> origin/master
        }
    }

    public void creatSpider(boolean boolRFL) throws IOException {
        spider = new BZPartReader[this.threadNum];
        for (int i = 0; i < this.threadNum; i++) {
<<<<<<< HEAD
            spider[i] = new BZPartReader(inputFile, i, intFileFormat, this.pos[i], this.pos[i + 1], '\t', boolRFL);
=======
            spider[i] = new BZPartReader(i, intFileFormat, this.pos[i], this.pos[i + 1], '\t', boolRFL);
>>>>>>> origin/master
            // System.out.println("BZPartReader" + i + ": created!");
        }
    }

    public void creatSpider(long longStart, long longEnd) throws IOException {
        spider = new BZPartReader[1];
<<<<<<< HEAD
        spider[0] = new BZPartReader(inputFile, 0, 0, longStart, longEnd, '\t');
=======
        spider[0] = new BZPartReader(0, 0, longStart, longEnd, '\t');
>>>>>>> origin/master
    }

    public void creatSpider(long longPos[]) throws IOException {
        spider = new BZPartReader[longPos.length - 1];
        for (int i = 0; i < spider.length; i++) {
<<<<<<< HEAD
            spider[i] = new BZPartReader(inputFile, i, 0, this.pos[i], this.pos[i + 1], '\t');
            // System.out.println("BZPartReader" + i + ": created!");
            if (i > 0) {
                spider[i - 1].bytEndLine = spider[i].getStartLine();
=======
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
>>>>>>> origin/master
            }
        }
    }

    private void creatBuilders(int intChrom, int intStart, boolean hasHead, String prefixName) throws IOException {
        builder = new BZPartBuilder[pos.length - 1];
        for (int i = 0; i < (pos.length - 1); i++) {
            builder[i] = new BZPartBuilder(inputFile, i, this.pos[i], this.pos[i + 1], intChrom, intStart, hasHead, prefixName);
            if (i > 0) {
                builder[i - 1].bzpr.bytEndLine = builder[i].bzpr.getStartLine();
            }
        }
    }

    public class ParseTask extends Task implements Callable<String> {

        BZPartReader bzprSpider;
        int intChrom;
        int intStart;
        long longCurr;

<<<<<<< HEAD
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
                    //System.out.println(longCurr);
                } else {
                    break;
=======
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
>>>>>>> origin/master
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

        creatBuilders(intChrom, intStart, hasHead, fleOutput.getCanonicalPath());

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

        String[] files = new String[runningThread];
        for (int j = 0; j < files.length; j++) {
            files[j] = fleOutput.getCanonicalPath() + "." + j;
        }
        concatenateFiles(files, fleOutput.getCanonicalPath());
        for (int j = 0; j < files.length; j++) {
            File f = new File(files[j]);
            f.delete();
        }
        /*
        BufferedWriter bw = new BufferedWriter(new FileWriter(fleOutput));
        for (BZPartBuilder builder1 : builder) {
            for (int j = 0; j < builder1.altIndex.size(); j++) {
                bw.write(builder1.altIndex.get(j));
            }
        }        
        bw.close();
         */

        index = fleOutput;

        return 0;
    }

    private long getNextStart(RandomAccessFile raf, long longCurr) throws IOException {
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
            for (int id = 0; id < intContent - 4; id++) {
                if (buf[id] == 31 && buf[id + 1] == -117 && buf[id + 2] == 8 && buf[id + 3] == 4 && buf[id + 4] == 0) { //This should be used unsigned number or others. 
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
        String lastChrom1 = null;
        String lastChromPos1 = null;
        boolean noMatchedChr = true;
        if (hasChromName) {
            List<String> tmpChrs = new ArrayList<String>();
            LongArrayList tmpIndex = new LongArrayList();
            while ((strLine = br.readLine()) != null) {
                String[] strItem = strLine.split("\t");
                //record the start and end of each chromsome
                if (lastChrom1 == null) {
                    lastChrom1 = strItem[0];
                    lastChromPos1 = strItem[2];
                    tmpChrs.add(lastChrom1);
                    tmpIndex.add(Long.valueOf(lastChromPos1));
                } else if (!lastChrom1.equals(strItem[0])) {
                    tmpChrs.add(lastChrom1);
                    tmpIndex.add(Long.valueOf(lastChromPos1));
                    lastChrom1 = strItem[0];
                    lastChromPos1 = strItem[2];
                    tmpChrs.add(lastChrom1);
                    tmpIndex.add(Long.valueOf(lastChromPos1));
                } else {
                    lastChromPos1 = strItem[2];
                }

                if (!strItem[0].equals(chrom)) {
                    lastChromPos = strItem[2];
                    if (boolFlag) {
                        ialtIndex.add(Integer.MAX_VALUE);
                        laltIndex.add(Long.valueOf(strItem[2]));
                        break;
                    }
                    continue;
                } else {
                    if (noMatchedChr) {
                        if (lastChromPos != null) {
                            ialtIndex.add(0);
                            laltIndex.add(Long.valueOf(lastChromPos));
                        }
                        noMatchedChr = false;
                    }
                    ialtIndex.add(Integer.valueOf(strItem[1]));
                    laltIndex.add(Long.valueOf(strItem[2]));
                    boolFlag = true;
                }
            }

            if (noMatchedChr) {
                String[] STAND_CHROM_NAMES = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                    "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "XY", "M", "Un"};
                int index = 0;
                for (index = 0; index < STAND_CHROM_NAMES.length; index++) {
                    if (STAND_CHROM_NAMES[index].equals(chrom)) {
                        break;
                    }
                }

                int index1 = index - 1;
                boolean hasFound = false;
                while (index1 >= 0) {
                    String testChr = STAND_CHROM_NAMES[index1];
                    for (int i = 0; i < tmpChrs.size(); i++) {
                        if (tmpChrs.get(i).equals(testChr)) {
                            ialtIndex.add(0);
                            laltIndex.add(tmpIndex.getQuick((i + 1) < tmpChrs.size() ? (i + 1) : i));
                            hasFound = true;
                            break;
                        }
                    }
                    if (hasFound) {
                        break;
                    }
                    index1--;
                }
                if (index1 < 0) {
                    ialtIndex.add(0);
                    laltIndex.add(0);
                }

                index1 = index + 1;
                hasFound = false;
                while (index1 < STAND_CHROM_NAMES.length) {
                    String testChr = STAND_CHROM_NAMES[index1];
                    for (int i = 0; i < tmpChrs.size(); i++) {
                        if (tmpChrs.get(i).equals(testChr)) {
                            ialtIndex.add(0);
                            laltIndex.add(tmpIndex.getQuick(i));
                            hasFound = true;
                            break;
                        }
                    }
                    if (hasFound) {
                        break;
                    }
                    index1++;
                }

            } else {
                tmpChrs.clear();
                tmpIndex.clear();
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
        if (ialtIndex == null || ialtIndex.isEmpty()) {
            return null;
        }
        if (laltIndex == null || laltIndex.isEmpty()) {
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

<<<<<<< HEAD
=======
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

>>>>>>> origin/master
        int posEnd = ialtIndex.binarySearch(end);
        if (posEnd >= 0) {
            posEnd = posEnd + 1;
            if (posEnd >= laltIndex.size()) {
<<<<<<< HEAD
                secPos[1] = inputFile.length();
=======
                secPos[1] = input.length();
>>>>>>> origin/master
            } else {
                secPos[1] = laltIndex.get(posEnd);
            }
        } else {
<<<<<<< HEAD
            posEnd = -posEnd - 1;
            if (posEnd >= laltIndex.size()) {
                secPos[1] = inputFile.length();
=======
            posEnd = -posEnd;
            if (posEnd >= laltIndex.size()) {
                secPos[1] = input.length();
>>>>>>> origin/master
            } else {
                secPos[1] = laltIndex.get(posEnd);
            }
        }
        return secPos;
    }

    public boolean checkIndex() throws IOException {
<<<<<<< HEAD
        if (!inputFile.exists()) {
            return false;
        }
        File fle = new File(inputFile.getCanonicalPath() + ".idx");
        if (!fle.exists()) {
            return false;
        } else {                                //one month later 
            Date indexFileDate = new Date(fle.lastModified());
            Date orgFileDate = new Date(inputFile.lastModified());
            //if too small or too early
            if (orgFileDate.after(indexFileDate)) {
                //an incomplete file
                fle.delete();
                return false;
            }
            index = fle;
            return true;
        }
    }

    public boolean concatenateFiles(String[] altFiles, String outFile) {
        try {
            try (FileChannel fclOutput = new FileOutputStream(outFile).getChannel()) {
                for (String strFile : altFiles) {
                    try (FileChannel fclInput = new FileInputStream(strFile).getChannel()) {
                        fclInput.transferTo(0, fclInput.size(), fclOutput);
                    }
                }
                fclOutput.force(true);
            }
            return true;
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
            return false;
        } catch (IOException ex) {
            ex.printStackTrace();
            return false;
        }
    }

    public static void main(String[] args) throws IOException {

        // String strFile = "E:\\home\\mxli\\MyJava\\kggseq3\\resources\\hg19\\hg19_snp141.txt.gz";
        // String strFile = "E:\\home\\mxli\\MyJava\\kggseq3\\resources\\hg19\\hg19_known_SNV_dbNCFP.gz";
        //String strFile = "E:\\home\\mxli\\MyJava\\kggseq3\\resources\\hg19\\hg19_all_SNV_dbNCFP.gz";
        // String strFile = "E:\\home\\mxli\\MyJava\\kggseq3\\resources\\hg19\\hg19_snp141.txt.gz";
        // String strFile = "E:\\home\\mxli\\MyJava\\kggseq3\\resources\\hg19\\hg19_known_SNV_dbNCFP.gz";
        //String strFile = "E:\\home\\mxli\\MyJava\\kggseq3\\resources\\hg19\\hg19_all_SNV_dbNCFP.gz";
         String strFile = "1kgafreur.20150813.flt.vcf.gz";
        //String strFile = "1kgafr.20150813.vcf.gz";
        //    String strFile="D:\\01WORK\\KGGseq\\testdata\\hg19_funcnote_encode_megamix.bed.gz.Histone.cmp.gz";
        // String strFile = "D:\\01WORK\\KGGseq\\Test\\test.txt.gz.aa";
//        String[] altFiles = new String[3];
//        altFiles[0] = ("D:\\01WORK\\KGGseq\\Test\\test.txt.gz.aa");
//        altFiles[1] = ("D:\\01WORK\\KGGseq\\Test\\test.txt.gz.ab");
//        altFiles[2] = ("D:\\01WORK\\KGGseq\\Test\\test.txt.gz.ac");
//        String outFile = "D:\\01WORK\\KGGseq\\Test\\test.txt.gz";
//        String strFile="D:\\01WORK\\KGGseq\\Test\\known_SNV_dbNCFP.new.shrink.txt.gz";
=======
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
>>>>>>> origin/master
        try {
//            InputStream is = new BufferedInputStream(new FileInputStream(strFile));
//            LimitInputStream cis;
//            cis = new LimitInputStream(is, 1411);
//            cis.skip(0);
<<<<<<< HEAD

            int test1 = 2;

            if (test1 == 1) {
//                System.out.println("CP1 Start -----> "+System.currentTimeMillis());
                BGZFInputStream bf = new BGZFInputStream(strFile, 4, true);
                //boolean flag = bf.concatenateFiles(altFiles, outFile);
                // System.out.println(flag);
=======
            int test1 = 2;
            if (test1 == 1) {
                BGZFInputStream bf = new BGZFInputStream(strFile, 4, true);
>>>>>>> origin/master
                if (!bf.checkIndex()) {
                    bf.adjustPos();
                    bf.buildIndex(strFile, 0, 1, false);
                }
<<<<<<< HEAD
//                System.out.println("CP1 Start -----> "+new Date());
                bf.readIndex(true, "1");
=======

                bf.readIndex(true, "X");
>>>>>>> origin/master
                long[] pos = bf.findIndex(20922786, 22192439);
                bf.adjustPos(pos[0], pos[1]);
                bf.creatSpider(pos[0] != 0);
                int[] temp = new int[500];
                byte[] bty;
                for (int i = 0; i < bf.threadNum; i++) {
                    while ((bty = bf.spider[i].readLine(temp)) != null) {
                        String str = new String(bty);
<<<<<<< HEAD
                        // System.out.println(temp[0] + " : " + str);
=======
                        System.out.println(temp[0] + " : " + str);
>>>>>>> origin/master
                        // System.out.println(str.substring(0, 20));
                    }
                }
            } else {
<<<<<<< HEAD
                BGZFInputStream bf = new BGZFInputStream(strFile, 4, true);
=======
                BGZFInputStream bf = new BGZFInputStream(strFile, 1, true);
>>>>>>> origin/master
                bf.adjustPos();
                bf.creatSpider();
                int[] temp = new int[500];
                byte[] bty;
<<<<<<< HEAD
                int t = 0;
                while ((bty = bf.spider[0].readLine(temp)) != null) {
                    // String str = new String(bty);
                    // System.out.println(temp[0] + " : " + str);
=======
                int t=0;
                while ((bty = bf.spider[0].readLine(temp)) != null) {
                   // String str = new String(bty);
                   // System.out.println(temp[0] + " : " + str);
>>>>>>> origin/master
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
