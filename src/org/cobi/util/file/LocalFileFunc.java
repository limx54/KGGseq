// (c) 2008-2009 Miaoxin Li
// This file is distributed as part of the IGG source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
// Permission is granted for you to use this file to compile IGG.
// All computer programs have bugs. Use this file at your own risk.
// Saturday, January 17, 2009
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.file;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.channels.ReadableByteChannel;
import java.nio.channels.WritableByteChannel;
import java.nio.charset.Charset;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.CharsetEncoder;
import java.nio.charset.CodingErrorAction;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;
import org.apache.log4j.Logger;
import org.cobi.util.coding.MD5File;

/**
 *
 * @author mxli
 */
public class LocalFileFunc {

    private final Logger LOG = Logger.getLogger(LocalFileFunc.class);

    static public void gunzipFile(String decryptedImageName, String tarredImageName) throws Exception {
        GZIPInputStream in = new GZIPInputStream(new FileInputStream(new File(decryptedImageName)));
        File outFile = new File(tarredImageName);
        try (ReadableByteChannel inChannel = Channels.newChannel(in); WritableByteChannel outChannel = new FileOutputStream(outFile).getChannel()) {
            ByteBuffer buffer = ByteBuffer.allocate(65536);
            while (inChannel.read(buffer) != -1) {
                buffer.flip();
                outChannel.write(buffer);
                buffer.clear();
            }
        }
    }

    static public String gunzipFile(String decryptedImageName) throws Exception {
        //Expect decryptedImageName ends with .gz
        String tarredImageName = decryptedImageName.substring(0, decryptedImageName.length() - 3);
        GZIPInputStream in = new GZIPInputStream(new FileInputStream(new File(decryptedImageName)));
        File outFile = new File(tarredImageName);
        try (ReadableByteChannel inChannel = Channels.newChannel(in); WritableByteChannel outChannel = new FileOutputStream(outFile).getChannel()) {
            ByteBuffer buffer = ByteBuffer.allocate(65536);
            while (inChannel.read(buffer) != -1) {
                buffer.flip();
                outChannel.write(buffer);
                buffer.clear();
            }
        }
        return tarredImageName;
    }

    public static void main(String[] args) {
        try {
            String[] STAND_CHROM_NAMES = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"};
            // String path = "E:\\home\\mxli\\MyJava\\kggseq3\\resources\\hg19\\hg19_mdbNSFP3.0.chr";
            String path = "E:\\home\\mxli\\MyJava\\kggseq3\\als\\";
            path = args[0];
            File folder = new File(path);
            if (!folder.isDirectory()) {
                System.err.println(path + " is not a director!");
                // return;
            }
            LocalFileFunc.bgzipInFolder(folder.getCanonicalPath(), true);

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    static public void bgzipFile(String decryptedImageName) throws Exception {
        String tarredImageName = decryptedImageName + ".tmp";
        File outFile = new File(tarredImageName);
        File inFile = new File(decryptedImageName);
        int len = 0;
        try (GZIPInputStream in = new GZIPInputStream(new FileInputStream(inFile))) {
            BlockCompressedOutputStream bcos = new BlockCompressedOutputStream(outFile);
            try (ReadableByteChannel inChannel = Channels.newChannel(in)) {
                ByteBuffer buffer = ByteBuffer.allocate(65536);
                while ((len = inChannel.read(buffer)) != -1) {
                    buffer.flip();
                    if (len < 65536) {
                        bcos.write(Arrays.copyOfRange(buffer.array(), 0, len));
                    } else {
                        bcos.write(buffer.array());
                    }
                    buffer.clear();
                }
                bcos.close();
            }
        }

        inFile.delete();
        outFile.renameTo(inFile);
        VCFFileReader vcfReader = new VCFFileReader(inFile, false);
        TabixIndexCreator creator = new TabixIndexCreator(new TabixFormat().VCF);
        int i = 0;
        for (VariantContext context : vcfReader) {
            creator.addFeature(context, i);
            i++;
        }
        vcfReader.close();
        creator.finalizeIndex(i).writeBasedOnFeatureFile(inFile);
    }

    static public void bgzipInFolder(String fileFolder, boolean filterGZ) {
        File folder = new File(fileFolder);
        String orgName = null;
        try {
            if (!folder.isDirectory()) {
                orgName = folder.getCanonicalPath();
                String tarredImageName = orgName + ".tmp";
                File outFile = new File(tarredImageName);
                File inFile = new File(orgName);
                int len = 0;
                System.out.println("Compressing " + orgName);
                try (GZIPInputStream in = new GZIPInputStream(new FileInputStream(inFile))) {
                    BlockCompressedOutputStream bcos = new BlockCompressedOutputStream(outFile);
                    try (ReadableByteChannel inChannel = Channels.newChannel(in)) {
                        ByteBuffer buffer = ByteBuffer.allocate(65536);
                        while ((len = inChannel.read(buffer)) != -1) {
                            buffer.flip();
                            if (len < 65536) {
                                bcos.write(Arrays.copyOfRange(buffer.array(), 0, len));
                            } else {
                                bcos.write(buffer.array());
                            }
                            buffer.clear();
                        }
                        bcos.close();
                    }
                }
                inFile.delete();
                outFile.renameTo(inFile);

                FileWriter wtr = new FileWriter(new File(orgName + ".md5"));
                String md5 = MD5File.getFileMD5StringStable(inFile);
                wtr.write(md5);
                wtr.close();
                return;
            }
            File[] files = folder.listFiles();
            for (File f : files) {

                try {
                    if (f.isDirectory()) {
                        bgzipInFolder(f.getCanonicalPath(), filterGZ);
                    }
                    orgName = f.getCanonicalPath();
                    if (filterGZ && !orgName.endsWith(".gz")) {
                        continue;
                    }
                } catch (Exception ex) {
                    System.err.println(ex.toString());
                }

                String tarredImageName = orgName + ".tmp";
                File outFile = new File(tarredImageName);
                File inFile = new File(orgName);
                int len = 0;
                System.out.println("Compressing " + orgName);
                try (GZIPInputStream in = new GZIPInputStream(new FileInputStream(inFile))) {
                    BlockCompressedOutputStream bcos = new BlockCompressedOutputStream(outFile);
                    try (ReadableByteChannel inChannel = Channels.newChannel(in)) {
                        ByteBuffer buffer = ByteBuffer.allocate(65536);
                        while ((len = inChannel.read(buffer)) != -1) {
                            buffer.flip();
                            if (len < 65536) {
                                bcos.write(Arrays.copyOfRange(buffer.array(), 0, len));
                            } else {
                                bcos.write(buffer.array());
                            }
                            buffer.clear();
                        }
                        bcos.close();
                    }
                }
                inFile.delete();
                outFile.renameTo(inFile);

                FileWriter wtr = new FileWriter(new File(orgName + ".md5"));
                String md5 = MD5File.getFileMD5StringStable(inFile);
                wtr.write(md5);
                wtr.close();

            }
        } catch (Exception ex) {
            System.err.println(ex.toString());
        }
    }

    static public void tabixFile(String decryptedImageName) {
        File inFile = new File(decryptedImageName);
        TabixIndexCreator creator;
        int i;
        try {
            try (VCFFileReader vcfReader = new VCFFileReader(inFile, false)) {
                creator = new TabixIndexCreator(new TabixFormat().VCF);
                i = 0;
                for (VariantContext context : vcfReader) {
                    creator.addFeature(context, i);
                    i++;
                }
            }
            creator.finalizeIndex(i).writeBasedOnFeatureFile(inFile);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    static public BufferedReader getBufferedReader(String filePath) throws Exception {
        BufferedReader br = null;
        File dataFile = new File(filePath);
        if (dataFile.exists()) {
            CharsetDecoder decoder = Charset.forName("UTF-8").newDecoder();
            decoder.onMalformedInput(CodingErrorAction.IGNORE);
            try {
                //br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(dataFile))));
                //This will work slightly better under multiplthread senarios
                br = new BufferedReader(Channels.newReader(Channels.newChannel(new GZIPInputStream(Channels.newInputStream(new RandomAccessFile(dataFile, "r").getChannel()))), decoder, 1024 * 1024));
            } catch (ZipException ex) {
                // br = new BufferedReader(new InputStreamReader(new FileInputStream(dataFile)));
                br = new BufferedReader(Channels.newReader(Channels.newChannel(Channels.newInputStream(new RandomAccessFile(dataFile, "r").getChannel())), decoder, 1024 * 1024));

            }
        } else {
            throw new Exception("No input file: " + dataFile.getCanonicalPath());
        }
        return br;
    }

    static public BufferedWriter getBufferedWriter(String filePath, boolean isGzip) throws Exception {
        BufferedWriter bw = null;
        File dataFile = new File(filePath);
//        if (!dataFile.getParentFile().exists()) {
//            dataFile.getParentFile().mkdirs();
//        }
        CharsetEncoder decoder = Charset.forName("UTF-8").newEncoder();
        decoder.onMalformedInput(CodingErrorAction.IGNORE);
        if (isGzip) {
            //This will work slightly better under multiplthread senarios
            bw = new BufferedWriter(Channels.newWriter(Channels.newChannel(new GZIPOutputStream(Channels.newOutputStream(new RandomAccessFile(dataFile, "rw").getChannel()))), decoder, 1024 * 1024));
        } else {
            // br = new BufferedWriter(new InputStreamWriter(new FileInputStream(dataFile)));
            bw = new BufferedWriter(Channels.newWriter(Channels.newChannel(Channels.newOutputStream(new RandomAccessFile(dataFile, "rw").getChannel())), decoder, 1024 * 1024));
        }

        return bw;
    }

    static public void delAll(File f) throws IOException {
        if (!f.exists()) {
            // throw new IOException("Cannot delete:" + f.getName());
            return;
        }
        boolean rslt = true;

        if (!(rslt = f.delete())) {
            File subs[] = f.listFiles();
            for (int i = 0; i <= subs.length - 1; i++) {
                if (subs[i].isDirectory()) {
                    delAll(subs[i]);
                }
                rslt = subs[i].delete();
            }
            rslt = f.delete();
        }

        if (!rslt) {
            throw new IOException("Cannot delete:" + f.getName());
        }
        return;

    }

    static public boolean delAllInside(File f) throws IOException {
        if (!f.exists()) {
            throw new IOException("Cannot delete:" + f.getName());
        }
        boolean rslt = false;
        if (f.isDirectory()) {
            File subs[] = f.listFiles();
            for (int i = 0; i < subs.length; i++) {
                if (subs[i].isDirectory()) {
                    delAll(subs[i]);
                }
                rslt = subs[i].delete();
            }
        } else {
            rslt = f.delete();
        }
        return rslt;
    }

    static public void copyFile(String srcFilename, String dstFilename) throws IOException {
        File f = new File(srcFilename);
        if (!f.exists()) {
            throw new IOException("File not exist:" + f.getName());
        }

        // Create channel on the source  
        FileChannel srcChannel = new FileInputStream(srcFilename).getChannel();

        // Create channel on the destination  
        FileChannel dstChannel = new FileOutputStream(dstFilename).getChannel();

        // Copy file contents from source to destination  
        dstChannel.transferFrom(srcChannel, 0, srcChannel.size());
        // Close the channels  
        srcChannel.close();
        dstChannel.close();
    }

    public static boolean makeStorageLoc(String filePath) throws Exception {
        if (filePath == null) {
            return false;
        }

        File nwFile = new File(filePath);
        if (!nwFile.exists()) {
            if (!nwFile.mkdirs()) {
                return false;
            }
        }
        return true;
    }
}
