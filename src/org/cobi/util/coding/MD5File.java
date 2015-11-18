/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.coding;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

/**
 *
 * @author  http://x7700.iteye.com/blog/584576
 */
public class MD5File {

    protected static char hexDigits[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};

    public static void main(String[] args) throws IOException, NoSuchAlgorithmException {
        long begin = System.currentTimeMillis();

        File big = new File("resources/snp.zip");
        String md5 = null;
        //  md5=   getFileMD5String(big);
        //String md5 = getMD5String("a");
        long end = System.currentTimeMillis();
        System.out.println("md5:" + md5 + " time:" + ((end - begin) / 1000) + "s");
        begin = System.currentTimeMillis();
        md5 = getFileMD5StringStable(big);
        //String md5 = getMD5String("a");
        end = System.currentTimeMillis();
        System.out.println("md5:" + md5 + " time:" + ((end - begin) / 1000) + "s");
    }

    public static String getFileMD5String(File file) throws IOException, NoSuchAlgorithmException {
        FileInputStream in = new FileInputStream(file);
        FileChannel ch = in.getChannel();

        //700000000 bytes are about 670M
        int maxSize = 700000000;
        MessageDigest messageDigest = MessageDigest.getInstance("MD5");
        long startPosition = 0L;
        long step = file.length() / maxSize;

        if (step == 0) {
            MappedByteBuffer byteBuffer = ch.map(FileChannel.MapMode.READ_ONLY, 0, file.length());
            messageDigest.update(byteBuffer);
            return bufferToHex(messageDigest.digest());
        }

        for (int i = 0; i < step; i++) {
            MappedByteBuffer byteBuffer = ch.map(FileChannel.MapMode.READ_ONLY, startPosition, maxSize);
            messageDigest.update(byteBuffer);
            startPosition += maxSize;
        }

        if (startPosition == file.length()) {
            return bufferToHex(messageDigest.digest());
        }

        MappedByteBuffer byteBuffer = ch.map(FileChannel.MapMode.READ_ONLY, startPosition, file.length() - startPosition);
        messageDigest.update(byteBuffer);


        return bufferToHex(messageDigest.digest());
    }

    public static String getFileMD5StringStable(File file) throws IOException, NoSuchAlgorithmException {
        int bufferSize = 1024 * 1024;
        byte buffer[] = new byte[bufferSize];
        MessageDigest messageDigest = MessageDigest.getInstance("MD5");
        InputStream fis = new FileInputStream(file);
        BufferedInputStream bis = new BufferedInputStream(fis, 1024 * 1024 * 10);
        while (null != bis && (bis.read(buffer)) != -1) {
            messageDigest.update(buffer);
        }
        bis.close();
        return bufferToHex(messageDigest.digest());
    }

    public static String getMD5String(String s) throws NoSuchAlgorithmException {
        return getMD5String(s.getBytes());
    }

    public static String getMD5String(byte[] bytes) throws NoSuchAlgorithmException {
        MessageDigest messageDigest = MessageDigest.getInstance("MD5");
        messageDigest.update(bytes);
        return bufferToHex(messageDigest.digest());
    }

    private static String bufferToHex(byte bytes[]) {
        return bufferToHex(bytes, 0, bytes.length);
    }

    private static String bufferToHex(byte bytes[], int m, int n) {
        StringBuffer stringbuffer = new StringBuffer(2 * n);
        int k = m + n;
        for (int l = m; l < k; l++) {
            appendHexPair(bytes[l], stringbuffer);
        }
        return stringbuffer.toString();
    }

    private static void appendHexPair(byte bt, StringBuffer stringbuffer) {
        char c0 = hexDigits[(bt & 0xf0) >> 4];
        char c1 = hexDigits[bt & 0xf];
        stringbuffer.append(c0);
        stringbuffer.append(c1);
    }

    public static boolean checkPassword(String password, String md5PwdStr) throws NoSuchAlgorithmException {
        String s = getMD5String(password);
        return s.equals(md5PwdStr);
    }
}
