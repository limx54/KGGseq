/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.file;

import cern.colt.list.ByteArrayList;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.io.Reader;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.channels.ReadableByteChannel;
import java.util.ArrayList;
import java.util.List;
import java.util.Observable;
import java.util.zip.GZIPInputStream;
import org.cobi.util.text.Util;

/**
 *
 * @author modified from Internet
 * http://blog.csdn.net/hbyscl/article/details/22923683
 */
public class ReadGZipFile extends Observable {

    private int bufSize = 1024 * 1024 * 10;
    // 换行符  
    private byte key = "\n".getBytes()[0];
    private byte key1 = "\t".getBytes()[0];
    // 当前行数  
    private long lineNum = 0;
    // 文件编码,默认为gb2312  
    private String encode = "utf-8";
    // 具体业务逻辑监听器  
    private ReaderFileListener readerListener;

    static class ProcessDataByPostgisListeners extends ReaderFileListener {

        public ProcessDataByPostgisListeners(String encode) {
            this.encode = encode;
        }

        @Override
        public void output(List<String> stringList) {
            for (String item : stringList) {
                Util.tokenize(item, '\t');
                // System.out.println(item);
            }

        }

    }

    public long getUncompresedSize(File dataFile) throws Exception {
        RandomAccessFile raf = new RandomAccessFile(dataFile, "r");
        raf.seek(raf.length() - 4);
        int b4 = raf.read();
        int b3 = raf.read();
        int b2 = raf.read();
        int b1 = raf.read();
        long availableUncompressedSize = (b1 << 24) | (b2 << 16) + (b3 << 8) + b4;
        raf.close();
        return availableUncompressedSize;
    }

    public static void main(String[] args) throws Exception {
        File file = new File("build.xml.gz");
        // File file = new File("SCZ_raw_redo_snp_recal_indel_recal.vcf.gz");

        FileInputStream fis = null;
        try {
            // Runtime runtime = Runtime.getRuntime();
            // int nrOfProcessors = runtime.availableProcessors();
            //  System.out.println("Number of processors available to the Java Virtual Machine: " + nrOfProcessors);

            ReadGZipFile readFile = new ReadGZipFile();
            fis = new FileInputStream(file);
            int available = fis.available();
            int maxThreadNum = 1;
            // 线程粗略开始位置  
            long time = System.nanoTime();
            int filePos = available / maxThreadNum;
            for (int j = 0; j < maxThreadNum; j++) {
                // 计算精确开始位置  
                long startNum = j == 0 ? 0 : readFile.getStartNum(file, filePos * j);
                long endNum = j + 1 < maxThreadNum ? readFile.getStartNum(file, filePos * (j + 1)) : -2;
                // 具体监听实现  
                //here the utf-8 is much faster than gb2312 
                ProcessDataByPostgisListeners listeners = new ProcessDataByPostgisListeners("utf-8");
                // ProcessDataByPostgisListeners listeners = new ProcessDataByPostgisListeners("gb2312");
                new ReadFileThread(listeners, startNum, endNum, file.getPath()).start();
            }

            time = System.nanoTime() - time;
            time = time / 1000000000;
            long min = time / 60;
            long sec = time % 60;
            System.out.println("Elapsed time: " + min + " min. " + sec + " sec.");
            time = System.nanoTime();
            /*
             BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
             String line;
             while (reader.ready()) {
             line = reader.readLine();
             Util.tokenize(line, '\t');
             }
             reader.close();
             time = System.nanoTime() - time;
             time = time / 1000000000;
             min = time / 60;
             sec = time % 60;
             System.out.println("Elapsed time: " + min + " min. " + sec + " sec.");
             */
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void setEncode(String encode) {
        this.encode = encode;
    }

    public void setReaderListener(ReaderFileListener readerListener) {
        this.readerListener = readerListener;
    }

    /**
     * 从设置的开始位置读取文件，一直到结束为止。如果 end设置为负数,刚读取到文件末尾
     *
     * @param fullPath
     * @param start
     * @param end
     * @throws Exception
     */
    public void readFileByLine(String fullPath, long start, long end) throws Exception {
        File fin = new File(fullPath);
        if (fin.exists()) {
            ReadableByteChannel fcin;
            FileChannel fc = new RandomAccessFile(fullPath, "r").getChannel();
            fcin = Channels.newChannel(new GZIPInputStream(Channels.newInputStream(fc)));
            // fcin = Channels.newChannel(new GZIPInputStream(Channels.newInputStream(Channels.newChannel(new FileInputStream(fin)))));
            // fcin.position(start);
            try {
                // MappedByteBuffer rBuffer = fc.map(FileChannel.MapMode.READ_ONLY, 0L, fc.size());
                ByteBuffer rBuffer = ByteBuffer.allocateDirect(bufSize);
                // 每次读取的内容  
                byte[] bs = new byte[bufSize];
                // 缓存  
                byte[] tempBs = new byte[0];
                String line = "";
                // 当前读取文件位置  
                long nowCur = start;
                int fromIndexCheck = 0;
                int fromIndexGetStr = 0;
                int endIndexCheck = 0;
                boolean isEnd = false;
                byte[] newStrByte;
                byte[] bLine;
                int tL;
                while (fcin.read(rBuffer) != -1) {
                    nowCur += bufSize;

                    int rSize = rBuffer.position();
                    rBuffer.rewind();
                    rBuffer.get(bs);
                    rBuffer.clear();
                    newStrByte = bs;
                    // 如果发现有上次未读完的缓存,则将它加到当前读取的内容前面  
                    if (null != tempBs) {
                        tL = tempBs.length;
                        newStrByte = new byte[rSize + tL];
                        System.arraycopy(tempBs, 0, newStrByte, 0, tL);
                        System.arraycopy(bs, 0, newStrByte, tL, rSize);
                    }
                    // 是否已经读到最后一位                     
                    // 如果当前读取的位数已经比设置的结束位置大的时候，将读取的内容截取到设置的结束位置  
                    if (end > 0 && nowCur > end) {
                        // 缓存长度 - 当前已经读取位数 - 最后位数  
                        int l = newStrByte.length - (int) (nowCur - end);
                        newStrByte = substring(newStrByte, 0, l);
                        isEnd = true;
                    }

                    //This is the most time consuming part. Java has down a lot of of optimization in its BufferedRead
                    // 每次读一行内容，以 key（默认为\n） 作为结束符  
                    while ((endIndexCheck = indexOf(newStrByte, fromIndexCheck)) != -1) {
                        bLine = substring(newStrByte, fromIndexGetStr, endIndexCheck);
                        line = new String(bLine, 0, bLine.length, encode);
                        //System.out.println(line);
                        lineNum++;
                        // 输出一行内容，处理方式由调用方提供  
                        //readerListener.outLine(line.trim(), lineNum, false);

                        fromIndexCheck = endIndexCheck + 1;
                        fromIndexGetStr = fromIndexCheck;
                    }
                    // 将未读取完成的内容放到缓存中  
                    tempBs = substring(newStrByte, fromIndexCheck, newStrByte.length);
                    fromIndexCheck = tempBs.length;
                    fromIndexGetStr = 0;
                    if (isEnd) {
                        break;
                    }
                }
                // 将剩下的最后内容作为一行，输出，并指明这是最后一行  
                String lineStr = new String(tempBs, 0, tempBs.length, encode);
                //System.out.println(lineStr);
                // readerListener.outLine(lineStr.trim(), lineNum, true);
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                fcin.close();
            }

        } else {
            throw new FileNotFoundException("没有找到文件：" + fullPath);
        }
        // 通知观察者,当前工作已经完成  
        setChanged();
        notifyObservers(start + "-" + end);
    }

    /**
     * 查找一个byte[]从指定位置之后的一个换行符位置
     *
     * @param src
     * @param fromIndex
     * @return
     * @throws Exception
     */
    private int indexOf(byte[] src, int fromIndex) throws Exception {
        for (int i = fromIndex; i < src.length; i++) {
            if (src[i] == key) {
                return i;
            }
        }
        return -1;
    }

    /**
     * 从指定开始位置读取一个byte[]直到指定结束位置为止生成一个全新的byte[]
     *
     * @param src
     * @param fromIndex
     * @param endIndex
     * @return
     * @throws Exception
     */
    private byte[] substring(byte[] src, int fromIndex, int endIndex) throws Exception {
        int size = endIndex - fromIndex;
        byte[] ret = new byte[size];
        System.arraycopy(src, fromIndex, ret, 0, size);
        return ret;
    }

    /**
     * 获取准确开始位置
     *
     * @param file
     * @param position
     * @return
     * @throws Exception
     */
    public long getStartNum(File file, long position) throws Exception {
        long startNum = position;
        FileChannel fcin = new RandomAccessFile(file, "r").getChannel();
        //  GZIPInputStream gis = new GZIPInputStream(Channels.newInputStream(fcin));
        Reader reader = new InputStreamReader(Channels.newInputStream(fcin));
        reader.skip(position);

        try {
            int cache = 1024;
            CharBuffer rBuffer = CharBuffer.allocate(cache);
            // 每次读取的内容  
            char[] bs = new char[cache];
            // 缓存  
            char[] tempBs = new char[0];
            String line = "";
            char[] newStrByte;
            int readLen = 0;
            while ((readLen = reader.read(rBuffer)) != -1) {
                int rSize = rBuffer.position();
                rBuffer.rewind();
                rBuffer.get(bs);
                rBuffer.clear();
                newStrByte = bs;
                // 如果发现有上次未读完的缓存,则将它加到当前读取的内容前面  
                if (null != tempBs) {
                    int tL = tempBs.length;
                    newStrByte = new char[rSize + tL];
                    System.arraycopy(tempBs, 0, newStrByte, 0, tL);
                    System.arraycopy(bs, 0, newStrByte, tL, rSize);
                }
                // 获取开始位置之后的第一个换行符  
                int endIndex = -1;
                for (int i = 0; i < newStrByte.length; i++) {
                    if (newStrByte[i] == key) {
                        endIndex = i;
                        break;
                    }
                }
                if (endIndex != -1) {
                    return position + endIndex;
                }

                tempBs = new char[newStrByte.length];
                System.arraycopy(newStrByte, 0, tempBs, 0, newStrByte.length);
                startNum += readLen;
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            fcin.close();
            //  gis.close();
            reader.close();
        }
        return startNum;
    }

    /**
     * 获取准确开始位置
     *
     * @param file
     * @param position
     * @return
     * @throws Exception
     */
    public long getStartNum1(File file, long position) throws Exception {
        long startNum = position;
        FileChannel fcin = new RandomAccessFile(file, "r").getChannel();
        GZIPInputStream gis = new GZIPInputStream(Channels.newInputStream(fcin));
        Reader reader = new InputStreamReader(gis);
        reader.skip(position);

        try {
            int cache = 1024;
            CharBuffer rBuffer = CharBuffer.allocate(cache);
            // 每次读取的内容  
            char[] bs = new char[cache];
            // 缓存  
            char[] tempBs = new char[0];
            String line = "";
            while (reader.read(rBuffer) != -1) {
                int rSize = rBuffer.position();
                rBuffer.rewind();
                rBuffer.get(bs);
                rBuffer.clear();
                char[] newStrByte = bs;
                // 如果发现有上次未读完的缓存,则将它加到当前读取的内容前面  
                if (null != tempBs) {
                    int tL = tempBs.length;
                    newStrByte = new char[rSize + tL];
                    System.arraycopy(tempBs, 0, newStrByte, 0, tL);
                    System.arraycopy(bs, 0, newStrByte, tL, rSize);
                }
                // 获取开始位置之后的第一个换行符  
                int endIndex = indexOf(newStrByte, 0);
                if (endIndex != -1) {
                    return startNum + endIndex;
                }
                tempBs = substring(newStrByte, 0, newStrByte.length);
                startNum += 1024;
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            fcin.close();
        }
        return position;
    }

    /**
     * 从设置的开始位置读取文件，一直到结束为止。如果 end设置为负数,刚读取到文件末尾
     *
     * @param fullPath
     * @param start
     * @param end
     * @throws Exception
     */
    public void readFileByLine1(String fullPath, long start, long end) throws Exception {
        File fin = new File(fullPath);
        if (fin.exists()) {
            FileChannel fcin = new RandomAccessFile(fin, "r").getChannel();
            GZIPInputStream gis = new GZIPInputStream(Channels.newInputStream(fcin));
            Reader reader = new InputStreamReader(gis);
            reader.skip(start);
            StringBuilder sb = new StringBuilder();
            try {
                bufSize = bufSize * bufSize;
                CharBuffer rBuffer = CharBuffer.allocate(bufSize);
                // 每次读取的内容  
                char[] bs = new char[bufSize];
                // 缓存  
                char[] tempBs = new char[0];
                String line = "";
                // 当前读取文件位置  
                long nowCur = start;
                while (reader.read(rBuffer) != -1) {
                    nowCur += bufSize;

                    int rSize = rBuffer.position();
                    rBuffer.rewind();
                    rBuffer.get(bs);
                    rBuffer.clear();
                    char[] newStrByte = bs;
                    // 如果发现有上次未读完的缓存,则将它加到当前读取的内容前面  
                    if (null != tempBs) {
                        int tL = tempBs.length;
                        newStrByte = new char[rSize + tL];
                        System.arraycopy(tempBs, 0, newStrByte, 0, tL);
                        System.arraycopy(bs, 0, newStrByte, tL, rSize);
                    }
                    // 是否已经读到最后一位  
                    boolean isEnd = false;
                    // 如果当前读取的位数已经比设置的结束位置大的时候，将读取的内容截取到设置的结束位置  
                    if (end > 0 && nowCur > end) {
                        // 缓存长度 - 当前已经读取位数 - 最后位数  
                        int l = newStrByte.length - (int) (nowCur - end);
                        newStrByte = substring(newStrByte, 0, l);
                        isEnd = true;
                    }
                    int fromIndex = 0;
                    int endIndex = 0;
                    // 每次读一行内容，以 key（默认为\n） 作为结束符  
                    while ((endIndex = indexOf(newStrByte, fromIndex)) != -1) {
                        char[] bLine = substring(newStrByte, fromIndex, endIndex);
                        sb.append(bLine);
                        line = sb.toString();
                        sb.delete(0, sb.length());
                        lineNum++;
                        // 输出一行内容，处理方式由调用方提供  
                        //  readerListener.outLine(line.trim(), lineNum, false);
                        fromIndex = endIndex + 1;
                    }
                    // 将未读取完成的内容放到缓存中  
                    tempBs = substring(newStrByte, fromIndex, newStrByte.length);
                    if (isEnd) {
                        break;
                    }
                }
                // 将剩下的最后内容作为一行，输出，并指明这是最后一行  
                String lineStr = new String(tempBs, 0, tempBs.length);
                //readerListener.outLine(lineStr.trim(), lineNum, true);
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                fcin.close();
            }

        } else {
            throw new FileNotFoundException("没有找到文件：" + fullPath);
        }
        // 通知观察者,当前工作已经完成  
        setChanged();
        notifyObservers(start + "-" + end);
    }
    private static ThreadLocal<String[]> tempArray = new ThreadLocal<String[]>();

    public void readFileByLine2(String fullPath, long start, long end) throws Exception {
        File fin = new File(fullPath);
        List<String[]> lineBuffer = new ArrayList<String[]>();
        if (fin.exists()) {
            FileChannel fcin = new RandomAccessFile(fin, "r").getChannel();
            GZIPInputStream gis = new GZIPInputStream(Channels.newInputStream(fcin));
            Reader reader = new InputStreamReader(gis);
            reader.skip(start);
            StringBuilder strBuilder = new StringBuilder();
            final char delimiter = '\t';

            int len = 0;
            try {
                CharBuffer rBuffer = CharBuffer.allocate(bufSize);
                // 每次读取的内容  
                char[] bs = new char[bufSize];

                String line = "";
                long nowCur = start;
                int acuLen = 0;
                long offset = 0;
                boolean isEnd = false;
                int fromIndexCheck = 0;
                int endIndexCheck = -1;
                int fromIndexGetStr = 0;
                // 是否已经读到最后一位  
                isEnd = false;
                int wordCount = 0;
                int wordStartI = 0;
                String[] temp = tempArray.get();
                int tempLength = 0;
                while ((acuLen = reader.read(rBuffer)) != -1) {
                    nowCur += acuLen;
                    rBuffer.rewind();
                    rBuffer.get(bs);
                    rBuffer.clear();
                    strBuilder.append(bs, 0, acuLen);
                    len = strBuilder.length();
                    // 如果当前读取的位数已经比设置的结束位置大的时候，将读取的内容截取到设置的结束位置  
                    if (end > 0 && nowCur >= end) {
                        offset = len - nowCur + end;
                        strBuilder.delete((int) offset, len);
                        isEnd = true;
                        len = strBuilder.length();
                    }

                    endIndexCheck = -1;
                    tempLength = (len / 2) + 1;
                    if (temp == null || temp.length < tempLength) {
                        if (temp != null) {
                            String[] temp1 = new String[temp.length];
                            System.arraycopy(temp, 0, temp1, 0, temp.length);
                            temp = new String[tempLength];
                            System.arraycopy(temp1, 0, temp, 0, temp1.length);
                        } else {
                            temp = new String[tempLength];
                        }
                        tempArray.set(temp);
                    }

                    do {
                        endIndexCheck = -1;
                        for (int i = fromIndexCheck; i < len; i++) {
                            if (strBuilder.charAt(i) == delimiter) {
                                temp[wordCount++] = strBuilder.substring(wordStartI, i);
                                wordStartI = i + 1;
                            } else if (strBuilder.charAt(i) == key) {
                                endIndexCheck = i;

                                if (wordStartI < i) {
                                    temp[wordCount++] = strBuilder.substring(wordStartI, i);
                                }
                                String[] result = new String[wordCount];
                                System.arraycopy(temp, 0, result, 0, wordCount);

                                // lineBuffer.add(result);
                                wordCount = 0;

                                wordStartI = i + 1;
                                break;
                            }
                        }
//end of the line
                        if (endIndexCheck >= 0) {
                            // line = strBuilder.substring(fromIndexGetStr, endIndexCheck);
                            // System.out.println(line);
                            lineNum++;

                            // 输出一行内容，处理方式由调用方提供  
                            //  readerListener.outLine(line.trim(), lineNum, false);
                            fromIndexCheck = endIndexCheck + 1;
                            fromIndexGetStr = fromIndexCheck;
                        }
                    } while (endIndexCheck != -1);

                    // 将未读取完成的内容放到缓存中  
                    strBuilder.delete(0, fromIndexGetStr);
                    wordStartI = wordStartI - fromIndexGetStr;
                    fromIndexCheck = strBuilder.length();
                    len = fromIndexCheck;
                    fromIndexGetStr = 0;
                    if (isEnd) {
                        break;
                    }
                }

                // 将剩下的最后内容作为一行，输出，并指明这是最后一行  
                if (wordStartI < len) {
                    temp[wordCount++] = strBuilder.substring(wordStartI);
                }
                if (wordCount > 0) {
                    String[] result = new String[wordCount];
                    System.arraycopy(temp, 0, result, 0, wordCount);
                    // lineBuffer.add(result);
                }

                // line = strBuilder.toString();
                //System.out.println(line);
                //readerListener.outLine(lineStr.trim(), lineNum, true);
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                fcin.close();
                gis.close();
                reader.close();
            }
        } else {
            throw new FileNotFoundException("Failed to find the file" + fullPath);
        }
       
        setChanged();
        notifyObservers(start + "-" + end);
    }

    /**
     * 查找一个byte[]从指定位置之后的一个换行符位置
     *
     * @param src
     * @param fromIndex
     * @return
     * @throws Exception
     */
    private int indexOf(char[] src, int fromIndex) throws Exception {

        for (int i = fromIndex; i < src.length; i++) {
            if (src[i] == key) {
                return i;
            }
        }
        return -1;
    }

    /**
     * 从指定开始位置读取一个byte[]直到指定结束位置为止生成一个全新的byte[]
     *
     * @param src
     * @param fromIndex
     * @param endIndex
     * @return
     * @throws Exception
     */
    private char[] substring(char[] src, int fromIndex, int endIndex) throws Exception {
        int size = endIndex - fromIndex;
        char[] ret = new char[size];
        System.arraycopy(src, fromIndex, ret, 0, size);
        return ret;
    }

}
