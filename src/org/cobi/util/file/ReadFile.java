/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.List;
import java.util.Observable;
import org.cobi.util.text.Util;

/**
 *
 * @author mxli
 */
public class ReadFile extends Observable {

    private int bufSize = 1024 * 1024;
    // 换行符  
    private byte key = "\n".getBytes()[0];
    // 当前行数  
    private long lineNum = 0;
    // 文件编码,默认为gb2312  
    // private String encode = "gb2312";
    private String encode = "utf-8";
    // 具体业务逻辑监听器  
    private ReaderFileListener readerListener;

    public void setEncode(String encode) {
        this.encode = encode;
    }

    public void setReaderListener(ReaderFileListener readerListener) {
        this.readerListener = readerListener;
    }

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

    public static void main(String[] args) throws Exception {
        File file = new File("als.flt.txt");

        FileInputStream fis = null;
        try {
            // Runtime runtime = Runtime.getRuntime();
            // int nrOfProcessors = runtime.availableProcessors();
            //  System.out.println("Number of processors available to the Java Virtual Machine: " + nrOfProcessors);

            ReadFile readFile = new ReadFile();
            fis = new FileInputStream(file);
            int available = fis.available();
            int maxThreadNum = 3;
            // 线程粗略开始位置  
            long time = System.nanoTime();
            int filePos = available / maxThreadNum;

            for (int j = 0; j < maxThreadNum; j++) {
                // 计算精确开始位置  
                long startNum = j == 0 ? 0 : readFile.getStartNum(file, filePos * j);
                long endNum = j + 1 < maxThreadNum ? readFile.getStartNum(file, filePos * (j + 1)) : -2;
                // 具体监听实现  
                //here the utf-8 is much faster than gb2312 
                ReadFile.ProcessDataByPostgisListeners listeners = new ReadFile.ProcessDataByPostgisListeners("utf-8");
                // ProcessDataByPostgisListeners listeners = new ProcessDataByPostgisListeners("gb2312");
                new ReadFileThread(listeners, startNum, endNum, file.getPath()).start();
            }
            time = System.nanoTime() - time;
            time = time / 1000000000;
            long min = time / 60;
            long sec = time % 60;
            System.out.println("Elapsed time: " + min + " min. " + sec + " sec.");
            time = System.nanoTime();
            //It turned out the BufferedReader has already been very fast;
            BufferedReader reader = new BufferedReader(new InputStreamReader((new FileInputStream(file))));
            String line;
            while (reader.ready()) {
                line = reader.readLine();
                // Util.tokenize(line, '\t');
            }
            reader.close();
            time = System.nanoTime() - time;
            time = time / 1000000000;
            min = time / 60;
            sec = time % 60;
            System.out.println("Elapsed time: " + min + " min. " + sec + " sec.");

        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
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
        fcin.position(position);
        try {
            int cache = 1024;
            ByteBuffer rBuffer = ByteBuffer.allocate(cache);
            // 每次读取的内容  
            byte[] bs = new byte[cache];
            // 缓存  
            byte[] tempBs = new byte[0];
            String line = "";
            while (fcin.read(rBuffer) != -1) {
                int rSize = rBuffer.position();
                rBuffer.rewind();
                rBuffer.get(bs);
                rBuffer.clear();
                byte[] newStrByte = bs;
                // 如果发现有上次未读完的缓存,则将它加到当前读取的内容前面  
                if (null != tempBs) {
                    int tL = tempBs.length;
                    newStrByte = new byte[rSize + tL];
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
    public void readFileByLine(String fullPath, long start, long end) throws Exception {
        File fin = new File(fullPath);
        if (fin.exists()) {
            FileChannel fcin = new RandomAccessFile(fin, "r").getChannel();
            fcin.position(start);
            try {
                // MappedByteBuffer rBuffer = fcin.map(FileChannel.MapMode.READ_ONLY,  0L, fcin.size());

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
                while (fcin.read(rBuffer) != -1) {
                    nowCur += bufSize;

                    int rSize = rBuffer.position();
                    rBuffer.rewind();
                    rBuffer.get(bs);
                    rBuffer.clear();
                    byte[] newStrByte = bs;
                    // 如果发现有上次未读完的缓存,则将它加到当前读取的内容前面  
                    if (null != tempBs) {
                        int tL = tempBs.length;
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
                        byte[] bLine = substring(newStrByte, fromIndexGetStr, endIndexCheck);
                        line = new String(bLine, 0, bLine.length, encode);
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

}
