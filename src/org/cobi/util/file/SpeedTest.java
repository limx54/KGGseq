/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.file;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;

import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;
import org.cobi.util.text.BGZFInputStream;

/**
 *
 * @author mxli
 */
public class SpeedTest {

    static String inFile = "SCZ_raw_redo_snp_recal_indel_recal.vcf.gz";
    //  static String inFile = "test.txt.gz";
    static String outFile = inFile + "tmp";

    public void convertPhased2Unphased(String filePath, String outPath) throws Exception {
        FileInputStream fin = new FileInputStream(filePath);
        InputStreamReader isr = null;
        BufferedReader br = null;
        File outFile = new File(outPath);
        BlockCompressedOutputStream bw = new BlockCompressedOutputStream(outFile);

        try {
            String line = null;
            GZIPInputStream gzis = new GZIPInputStream(fin);
            isr = new InputStreamReader(gzis, "utf-8");
            br = new BufferedReader(isr, 1024 * 1024);
            String tmpStr;
            while ((line = br.readLine()) != null) {
                //  String[] cells = line.split("\t");
                bw.write(line.replace('|', '/').getBytes());
                bw.write("\n".getBytes());
            }
            bw.close();
            fin.close();
            isr.close();
            br.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void parsingSpeed(String filePath) throws Exception {
        FileInputStream fin = new FileInputStream(filePath);
        InputStreamReader isr = null;
        BufferedReader br = null;

        try {
            String line = null;
            GZIPInputStream gzis = new GZIPInputStream(fin);
            isr = new InputStreamReader(gzis, "utf-8");
            br = new BufferedReader(isr, 1024 * 1024);
            String tmpStr;
            while ((line = br.readLine()) != null) {
                // String[] cells = line.split("\t");
                // bw.write(line + "\n");
                StringTokenizer st = new StringTokenizer(line, "\t");
                while (st.hasMoreTokens()) {
                    tmpStr = st.nextToken();
                }
            }
            fin.close();
            isr.close();
            br.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void parsingSpeedByte(String filePath) throws Exception {
        try {
            byte[] line = null;
            BGZFInputStream bf = new BGZFInputStream(filePath, 1);
            bf.adjustPos();
            bf.creatSpider();
            int[] cellDelimts = new int[2500];
            while ((line = bf.spider[0].readLine(cellDelimts)) != null) {

                // String[] cells = line.split("\t");
                // bw.write(line + "\n");
            }
            bf.spider[0].closeInputStream();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //Conculstion: I spent one day in comparing performance of io and nio. I found io is much faster than nio on Windows
    //Probaly io has been optimized a lot for text parsing.
    public static void main(String[] args) {
        // multiply2();
        //byteIntTest();
        try {
            SpeedTest st = new SpeedTest();
            Long startTime = System.nanoTime();
            //readWriteByIo();
            String inFilePath = "D:\\kgg3d5\\resources\\LD\\v5\\all\\ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
            String outFilePath = "D:\\kgg3d5\\resources\\LD\\v5\\all\\ALL.chr1.unphase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
            st.convertPhased2Unphased(inFilePath, outFilePath);
            /*
            st.parsingSpeed(filePath);
            Long endTime = System.nanoTime();
            System.out.println("用时：" + (endTime - startTime) / 1000000000.0);
            startTime = System.nanoTime();
            //st.parsingSpeedByte(filePath);
            endTime = System.nanoTime();
            System.out.println("用时：" + (endTime - startTime) / 1000000000.0);
             */
        } catch (Exception e) {
            e.printStackTrace();
        }
        /*
         Long startTime = System.nanoTime();
         //readWriteByIo();
         readByIo();
         Long endTime = System.nanoTime();
         System.out.println("用时：" + (endTime - startTime) / 1000000000.0);
         startTime = System.nanoTime();
         readByNio();
         endTime = System.nanoTime();
         System.out.println("用时：" + (endTime - startTime) / 1000000000.0);

         startTime = System.nanoTime();
         // FileInputStream fis = new FileInputStream(inFile);

         //long available = fis.available();
         ReadFile rf = new ReadFile();
         // rf.getStartNum(new File(inFile), available / 2);
         rf.readFileByLine(inFile, 0, -1);

         ReadGZipFile rgz = new ReadGZipFile();
         //long start = rgz.getStartNum(new File(inFile), available / 2);
         //rgz.readFileByLine(inFile, 0, -1);
         endTime = System.nanoTime();
         System.out.println("用时：" + (endTime - startTime) / 1000000000.0);
         */
    }

    public static void byteIntTest() throws Exception {
        int n = 300000000;
        long start = System.nanoTime();
        int a = 3233323;
        byte b = 124;
        boolean x;
        Map<Integer, boolean[]> unphasedGtyCodingMapInt = new HashMap<Integer, boolean[]>();
        Map<Byte, boolean[]> unphasedGtyCodingMapByte = new HashMap<Byte, boolean[]>();
        for (int i = 0; i < 100; i++) {
            unphasedGtyCodingMapInt.put(i, null);
            unphasedGtyCodingMapByte.put((byte) i, null);
        }
        start = System.nanoTime();
        for (int i = 0; i < n; i++) {
            //   a = a | GlobalManager.intOpers[4];
            //  a = a >> 3;
            // x = (a & GlobalManager.intOpers[16]) == GlobalManager.intOpers[16];
            unphasedGtyCodingMapInt.containsKey(i);
        }

        System.out.println((System.nanoTime() - start) / 1000);
        start = System.nanoTime();
        for (int i = 0; i < n; i++) {
            //  b = (byte) (b | GlobalManager.byteOpers[4]);
            //  b = (byte) (b >> 3);
            // x = (b & GlobalManager.byteOpers[4]) == GlobalManager.byteOpers[4];
            unphasedGtyCodingMapByte.containsKey((byte) i);
        }
        System.out.println((System.nanoTime() - start) / 1000);
    }

    public static void readWriteFileByNio() throws Exception {

        // 获取源文件和目标文件的输入输出流  
        FileInputStream fin = new FileInputStream(inFile);
        FileOutputStream fout = new FileOutputStream(outFile);
        // 获取输入输出通道  
        FileChannel fcin = fin.getChannel();
        FileChannel fcout = fout.getChannel();
        /*
         // 创建缓冲区  
         ByteBuffer ib = fcin.map(FileChannel.MapMode.READ_ONLY, 0, fcin.size()).asReadOnlyBuffer(); 
         while (ib.hasRemaining()) {
         ib.get();
         }
         fcin.close();
         */
        ByteBuffer buffer = ByteBuffer.allocate(1024 * 1024);
        while (true) {
            // clear方法重设缓冲区，使它可以接受读入的数据  
            buffer.clear();
            // 从输入通道中将数据读到缓冲区  
            int r = fcin.read(buffer);
            // read方法返回读取的字节数，可能为零，如果该通道已到达流的末尾，则返回-1  
            if (r == -1) {
                break;
            }
            // flip方法让缓冲区可以将新读入的数据写入另一个通道  
            buffer.flip();
            // 从输出通道中将数据写入缓冲区  
            fcout.write(buffer);
        }

    }

    public static void readWriteByIo() throws FileNotFoundException {
        // 获取源文件和目标文件的输入输出流  
        FileInputStream fin = new FileInputStream(inFile);
        FileOutputStream fout = new FileOutputStream(outFile);
        InputStreamReader isr = null;
        OutputStreamWriter osw = null;
        BufferedReader br = null;
        BufferedWriter bw = null;
        try {
            String line = null;
            GZIPInputStream gzis = new GZIPInputStream(fin);
            isr = new InputStreamReader(gzis, "utf-8");
            br = new BufferedReader(isr, 1024 * 1024);

            osw = new OutputStreamWriter(fout, "utf-8");
            bw = new BufferedWriter(osw, 1024 * 1024);

            while ((line = br.readLine()) != null) {
                // bw.write(line + "\n");
            }
            bw.flush();
            fin.close();
            fout.close();
            isr.close();
            osw.close();
            br.close();
            bw.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void multiply2() throws FileNotFoundException {
        int size = 2000000000;
        int s = 0;
        Long startTime = System.nanoTime();
        for (int i = 0; i < size; i += 2) {
            s = i * 2;
            // s = (i + 1) * 2;
        }

        Long endTime = System.nanoTime();
        System.out.println("用时：" + (endTime - startTime) / 1000000000.0);
        s = 0;
        startTime = System.nanoTime();
        for (int i = 0; i < size; i += 2) {
            s = i << 1;
            //s = (i + 1) << 1;
        }
        endTime = System.nanoTime();
        System.out.println("用时：" + (endTime - startTime) / 1000000000.0);
    }

    public static void readByIo() throws FileNotFoundException {
        // 获取源文件和目标文件的输入输出流  
        FileInputStream fin = new FileInputStream(inFile);
        InputStreamReader isr = null;
        BufferedReader br = null;

        try {
            String line = null;
            if (inFile.endsWith(".gz")) {
                GZIPInputStream gzis = new GZIPInputStream(fin);
                isr = new InputStreamReader(gzis, "utf-8");
                br = new BufferedReader(isr, 1024 * 1024);
            } else {
                isr = new InputStreamReader(fin, "utf-8");
                br = new BufferedReader(isr);
            }
            LineNumberReader reader = new LineNumberReader(br);
            while ((line = br.readLine()) != null) {
                // bw.write(line + "\n");
            }

            fin.close();
            isr.close();
            br.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void readByNio() throws FileNotFoundException {
        // 获取源文件和目标文件的输入输出流  

        LineNumberReader br = null;
        try {
            String line = null;
            if (inFile.endsWith(".gz")) {
                br = new LineNumberReader(Channels.newReader(Channels.newChannel(new GZIPInputStream(Channels.newInputStream(Channels.newChannel(new FileInputStream(inFile))))), Charset.forName("UTF-8").newDecoder(), 1024 * 1024 * 10));
            } else {
                br = new LineNumberReader(Channels.newReader(new FileInputStream(inFile).getChannel(), Charset.forName("UTF-8").newDecoder(), 1024 * 1024));
            }

            while ((line = br.readLine()) != null) {
                // bw.write(line + "\n");
            }

            br.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
