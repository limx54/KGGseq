/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.download.stable;

import java.io.InputStream;
import java.io.RandomAccessFile;
import java.util.concurrent.Callable;
import org.apache.http.Header;
import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.DefaultHttpClient;

/**
 *
 */
public class HttpClient4DownloadThread extends DownloadThread implements Callable {

    public HttpClient4DownloadThread(int sectionID, DownloadTaskBean taskBean) {
        super(sectionID, taskBean);
    }

    @Override
    public String call() throws Exception {
        if (HttpClient4DownloadTask.getDebug()) {
            System.out.println("Thread:" + sectionID + " Start:" + taskBean.getSectionsOffset()[sectionID] + "-" + taskBean.getSectionsEnd()[sectionID]);
        }
        HttpClient httpClient = new DefaultHttpClient();
        try {
            HttpGet httpGet = new HttpGet(taskBean.getDownURL());
            //  System.out.println(httpGet.getURI().toString());
            if (taskBean.isIsRange()) {//���߳�����
                httpGet.addHeader("Range", "bytes=" + taskBean.getSectionsOffset()[sectionID] + "-" + taskBean.getSectionsEnd()[sectionID]);
            }
            HttpResponse response = httpClient.execute(httpGet);
            int statusCode = response.getStatusLine().getStatusCode();
            if (HttpClient4DownloadTask.getDebug()) {
                for (Header header : response.getAllHeaders()) {
                    System.out.println(header.getName() + ":" + header.getValue());
                }
                System.out.println("statusCode:" + statusCode);
            }
            if (statusCode == 206 || (statusCode == 200 && !taskBean.isIsRange())) {
                InputStream inputStream = response.getEntity().getContent();
                // BufferedInputStream bis = new BufferedInputStream(is, temp.length);
                //��������д��
                RandomAccessFile outputStream = new RandomAccessFile(taskBean.getSaveFile() + ".save", "rw");
                long offset = taskBean.getSectionsOffset()[sectionID];
                /*
                //it does not work. i do not know why
                 FileChannel channel = outputStream.getChannel();
                 long size = taskBean.getSectionsEnd()[sectionID] - offset;
                 final FileLock lock = channel.tryLock(offset, size, false);
                 // final MappedByteBuffer buffer1 = channel.map(MapMode.READ_WRITE, offset, size);
                 if (lock != null) {
                 channel.transferFrom(Channels.newChannel(inputStream), offset, size);
                 lock.release();
                 }
                 channel.close()
                 */;

                outputStream.seek(taskBean.getSectionsOffset()[sectionID]);
                int count = 0;
                byte[] buffer = new byte[10 * 1024];
                while (notStop && (count = inputStream.read(buffer, 0, buffer.length)) > 0) {
                    outputStream.write(buffer, 0, count);
                    offset += count;
                    taskBean.getSectionsOffset()[sectionID] = offset;
                    //���������¼�
                    fireAfterPerDown(new DownloadThreadEvent(this, count));
                }

                outputStream.close();
                if (notStop) {
                    response.getEntity().consumeContent();
                }
            }
            httpGet.abort();
        } finally {
            fireDownCompleted(new DownloadThreadEvent(this, taskBean.getSectionsEnd()[sectionID]));
            if (HttpClient4DownloadTask.getDebug()) {
                System.out.println("End:" + taskBean.getSectionsOffset()[sectionID] + "-" + taskBean.getSectionsEnd()[sectionID]);
            }
            httpClient.getConnectionManager().shutdown();
            return "Download thread " + sectionID + " over";
        }
    }
}
