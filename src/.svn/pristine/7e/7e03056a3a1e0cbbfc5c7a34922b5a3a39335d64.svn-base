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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.CharsetEncoder;
import java.nio.charset.CodingErrorAction;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

/**
 *
 * @author mxli
 */
public class LocalFileFunc {

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
