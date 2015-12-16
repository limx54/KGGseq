/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.download.stable;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author mxli
 */
public class DownloadTask {

    protected String url = "";
    protected int threadNum = 5;
    protected String localPath = "";
    protected long receivedCount = 0;
    protected List<DownloadThread> threads = Collections.synchronizedList(new LinkedList<DownloadThread>());
    protected long lastCount = 0;
    protected long beginTime = 0;
    protected long endTime = 0;
    protected long autoCallbackSleep = 1000;
    protected List<DownloadTaskListener> listeners = Collections.synchronizedList(new LinkedList<DownloadTaskListener>());
    DownloadTaskBean taskBean = null;
    RandomAccessFile taskRandomFile = null;
    private static boolean DEBUG = false;
    protected final Object object = new Object();
    protected boolean done = false;
    protected String dataMd5=null;
    protected static int repeatTime=0;

    public String getDataMd5() {
        return dataMd5;
    }

    public void setDataMd5(String dataMd5) {
        this.dataMd5 = dataMd5;
    }

    
    public boolean isDone() {
        return done;
    }

    public static void setDebug(boolean debug) {
        DEBUG = debug;
    }

    public static boolean getDebug() {
        return DEBUG;
    }

    public void setAutoCallbackSleep(long autoCallbackSleep) {
        this.autoCallbackSleep = autoCallbackSleep;
    }

    public long getAutoCallbackSleep() {
        return this.autoCallbackSleep;
    }

    public void setLocalPath(String localPath) {
        this.localPath = localPath;
    }

    public String getLocalPath() {
        return localPath;
    }

    public void addTaskListener(DownloadTaskListener listener) {
        listeners.add(listener);
    }

    public void cancel()  {
        for (DownloadThread download : threads) {
            download.stop();
        } 
    }

    class Moniter implements Runnable {

        long availalableCount;
        long totalContentLength;

        public Moniter(long availalableCount, long totalContentLength) {
            this.availalableCount = availalableCount;
            this.totalContentLength = totalContentLength;
        }

        @Override
        public void run() {
            if (receivedCount < totalContentLength && !threads.isEmpty()) {
                showInfo(availalableCount, totalContentLength);
            } else {
                showInfo(availalableCount, totalContentLength);
            }
        }
    }

    protected void showInfo(long availalableCount, long totalContentLength) {
        long currentTime = System.currentTimeMillis();
        double realTimeSpeed = (receivedCount - lastCount) * 1.0 / ((currentTime - endTime) / 1000.0);
        double globalSpeed = receivedCount * 1.0 / ((currentTime - beginTime) / 1000.0);
        lastCount = receivedCount;
        endTime = currentTime;
        fireAutoCallback(new DownloadTaskEvent(receivedCount, availalableCount, totalContentLength, formatSpeed(realTimeSpeed), formatSpeed(globalSpeed), done));
    }

    private void fireAutoCallback(DownloadTaskEvent event) {
        if (listeners.isEmpty()) {
            return;
        }
        for (DownloadTaskListener listener : listeners) {
            listener.autoCallback(event);
        }
    }

    private String formatSpeed(double speed) {
        DecimalFormat format = new DecimalFormat("#,##0.##");
        if (speed < 1024) {
            return format.format(speed) + " B/s";
        }

        speed /= 1024;
        if (speed < 1024) {
            return format.format(speed) + " K/s";
        }

        speed /= 1024;
        if (speed < 1024) {
            return format.format(speed) + " M/s";
        }

        speed /= 1024;
        if (speed < 1024) {
            return format.format(speed) + " G/s";
        }

        speed /= 1024;
        if (speed < 1024) {
            return format.format(speed) + " T/s";
        }

        return format.format(speed) + "B/s";
    }

    //read data from the description file
    public synchronized DownloadTaskBean readTaskBean(RandomAccessFile file) throws IOException {
        DownloadTaskBean taskBeanTmp = new DownloadTaskBean();
        byte[] temp = new byte[DownloadTaskBean.HEAD_SIZE];
        file.seek(0);
        int readed = file.read(temp);
        if (readed != temp.length) {
            throw new RuntimeException();
        }

        ByteArrayInputStream bais = new ByteArrayInputStream(temp);
        DataInputStream dis = new DataInputStream(bais);

        taskBeanTmp.setDownURL(dis.readUTF());
        taskBeanTmp.setSaveFile(dis.readUTF());
        taskBeanTmp.setSectionCount(dis.readInt());
        taskBeanTmp.setContentLength(dis.readLong());
        taskBeanTmp.setIsRange(dis.readBoolean());

        bais.close();
        dis.close();

        file.seek(DownloadTaskBean.HEAD_SIZE);
        int sectionCount = taskBeanTmp.getSectionCount();
        long[] sectionsEnd = new long[sectionCount];
        for (int i = 0; i < sectionCount; i++) {
            sectionsEnd[i] = file.readLong();
        }
        taskBeanTmp.setSectionsEnd(sectionsEnd);

        long[] sectionsOffset = new long[sectionCount];
        for (int i = 0; i < sectionCount; i++) {
            sectionsOffset[i] = file.readLong();
        }
        taskBeanTmp.setSectionsOffset(sectionsOffset);

        return taskBeanTmp;
    }

    protected DownloadTaskBean createTaskBean(boolean acceptRanges, long dwonloadContentLength) throws IOException {
        //inconsistent length
        DownloadTaskBean taskBeanTmp = new DownloadTaskBean();
        taskBeanTmp.setDownURL(url);
        taskBeanTmp.setSaveFile(localPath);
        taskBeanTmp.setContentLength(dwonloadContentLength);
        taskBeanTmp.setSectionCount(threadNum);
        taskBeanTmp.setIsRange(acceptRanges);


        long perThreadLength = dwonloadContentLength / threadNum + 1;
        long startPosition = 0;
        long endPosition = perThreadLength;
        long[] sectionsOffset = new long[threadNum];
        long[] sectionsEnd = new long[threadNum];
        int secCount = 0;
       
        while (secCount < threadNum) {
            endPosition = startPosition + perThreadLength;
            if (endPosition > dwonloadContentLength) {
                endPosition = dwonloadContentLength;
            }
            sectionsOffset[secCount] = startPosition;
            sectionsEnd[secCount] = endPosition;
            startPosition = endPosition + 1;//�˴��� 1,�ӽ���λ�õ���һ���ط���ʼ����
            secCount++;
        }
        taskBeanTmp.setSectionsOffset(sectionsOffset);
        taskBeanTmp.setSectionsEnd(sectionsEnd);

        return taskBeanTmp;
    }

    //create description file
    public synchronized void writeTaskBean(RandomAccessFile file, DownloadTaskBean taskBean) throws IOException {
        long len = taskBean.getTaskFileLen();
        file.setLength(len);
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        DataOutputStream dos = new DataOutputStream(baos);
        dos.writeUTF(taskBean.getDownURL());
        dos.writeUTF(taskBean.getSaveFile());
        dos.writeInt(taskBean.getSectionCount());
        dos.writeLong(taskBean.getContentLength());
        dos.writeBoolean(taskBean.isIsRange());

        byte[] src = baos.toByteArray();
        byte[] temp = new byte[DownloadTaskBean.HEAD_SIZE];
        System.arraycopy(src, 0, temp, 0, src.length);
        file.seek(0);
        file.write(temp);
        baos.close();
        dos.close();
        taskRandomFile.seek(DownloadTaskBean.HEAD_SIZE);
        long[] sectionsEnd = taskBean.getSectionsEnd();
        for (int i = 0; i < sectionsEnd.length; i++) {
            file.writeLong(sectionsEnd[i]);
        }
        writeOffsetTaskBean(file, taskBean);

    }
//update the downloading process

    public synchronized void writeOffsetTaskBean(RandomAccessFile taskRandomFile, DownloadTaskBean taskBean) throws IOException {
        long[] sectionsOffset = taskBean.getSectionsOffset();
        taskRandomFile.seek(DownloadTaskBean.HEAD_SIZE + 8 * sectionsOffset.length);
        for (int i = 0; i < sectionsOffset.length; i++) {
            taskRandomFile.writeLong(sectionsOffset[i]);
        }
    }
}
