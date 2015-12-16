/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.file;

import com.ice.tar.TarEntry;
import com.ice.tar.TarInputStream;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author mxli
 */
public class Tar {

    public static void untargz(String tarFileName, String untarDir) throws IOException {
        InputStream in;
        //if (tarFileName.substring(tarFileName.lastIndexOf(".") + 1, tarFileName.lastIndexOf(".") + 3).equalsIgnoreCase("gz")) 
        if (tarFileName.endsWith("tar.gz") || tarFileName.endsWith("gz.tgz")) {
            System.out.println("Creating an GZIPInputStream for the file");
            in = new GZIPInputStream(new FileInputStream(new File(tarFileName)));

        } else {
            System.out.println("Creating an InputStream for the file");
            in = new FileInputStream(new File(tarFileName));
        }
        System.out.println("Reading TarInputStream... ");
        TarInputStream tin = new TarInputStream(in);
        TarEntry tarEntry = tin.getNextEntry();
        if (new File(untarDir).exists()) {
            while (tarEntry != null) {
                File destPath = new File(untarDir + File.separatorChar + tarEntry.getName());
                System.out.println("Processing " + destPath.getAbsoluteFile());
                if (!tarEntry.isDirectory()) {
                    FileOutputStream fout = new FileOutputStream(destPath);
                    tin.copyEntryContents(fout);
                    fout.close();
                } else {
                    destPath.mkdir();
                }
                tarEntry = tin.getNextEntry();
            }
            tin.close();
        } else {
            System.out.println("That destination directory doesn't exist! " + untarDir);
        }

    }

    
    public static void untar(String tarFileName, String untarDir) throws IOException {
        InputStream in;
        //if (tarFileName.substring(tarFileName.lastIndexOf(".") + 1, tarFileName.lastIndexOf(".") + 3).equalsIgnoreCase("gz")) 
        if (tarFileName.endsWith(".tar") ) {
            System.out.println("Creating an GZIPInputStream for the file");
            in = new BufferedInputStream(new FileInputStream(new File(tarFileName)));

        } else {
            System.out.println("Creating an InputStream for the file");
            in = new FileInputStream(new File(tarFileName));
        }
        System.out.println("Reading TarInputStream... ");
        TarInputStream tin = new TarInputStream(in);
        TarEntry tarEntry = tin.getNextEntry();
        if (new File(untarDir).exists()) {
            while (tarEntry != null) {
                File destPath = new File(untarDir + File.separatorChar + tarEntry.getName());
                System.out.println("Processing " + destPath.getAbsoluteFile());
                if (!tarEntry.isDirectory()) {
                    FileOutputStream fout = new FileOutputStream(destPath);
                    tin.copyEntryContents(fout);
                    fout.close();
                } else {
                    destPath.mkdir();
                }
                tarEntry = tin.getNextEntry();
            }
            tin.close();
        } else {
            System.out.println("That destination directory doesn't exist! " + untarDir);
        }

    }

    public static void main(String[] strArgs) throws IOException {
        // String strSourceFile = "resources/phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ASN.vcf.gz.tgz";
        String strSourceFile = "resources/hapmap2.r22.chbjpt.hg19.tar";
        // String strSourceFile = "resources/kgg.tar.gz";

        String strDest = "resources/";
        try {
            // untargz(strSourceFile, strDest);
            untar(strSourceFile, strDest);
            //unTar(new File(strSourceFile), new File(strDest));
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }
}
