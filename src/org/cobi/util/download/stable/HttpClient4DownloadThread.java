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
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;

/**
 *
 */
public class HttpClient4DownloadThread extends DownloadThread implements Callable {

    public HttpClient4DownloadThread(int sectionID, DownloadTaskBean taskBean) {
        super(sectionID, taskBean);
    }

    @Override
    public String call() {
        if (HttpClient4DownloadTask.getDebug()) {
            System.out.println("Thread:" + sectionID + " Start:" + taskBean.getSectionsOffset()[sectionID] + "-" + taskBean.getSectionsEnd()[sectionID]);
        }
        //HttpClient httpClient = HttpClientBuilder.create().build();
        CloseableHttpClient httpClient = HttpClients.createDefault();
        // RequestConfig requestConfig = RequestConfig.custom().setSocketTimeout(5000).setConnectTimeout(5000).setConnectionRequestTimeout(5000).build();

        long downloadLen = 0;
        long startSite = taskBean.getSectionsOffset()[sectionID];
        try {

            HttpGet httpGet = new HttpGet(taskBean.getDownURL());
            // httpGet.setConfig(requestConfig);
            //httpGet.setHeader("User-Agent", "Mozilla/5.0 (X11; U; Linux x86_64; en-US; rv:1.9.2.13) Gecko/20101206 Firefox/3.6.13");
            //httpGet.setHeader("Content-Type", "application/x-gzip");

            //  System.out.println(httpGet.getURI().toString());
            if (taskBean.isIsRange()) {//���߳�����
                httpGet.addHeader("Range", "bytes=" + taskBean.getSectionsOffset()[sectionID] + "-" + taskBean.getSectionsEnd()[sectionID]);
            }

// Execution context can be customized locally.
            // HttpClientContext context = HttpClientContext.create();
            // Contextual attributes set the local context level will take
            // precedence over those set at the client level.
            //context.setAttribute("http.protocol.version", HttpVersion.HTTP_1_1);
            HttpResponse response = httpClient.execute(httpGet);
            int statusCode = response.getStatusLine().getStatusCode();
            if (HttpClient4DownloadTask.getDebug()) {
                for (Header header : response.getAllHeaders()) {
                    System.out.println(header.getName() + ":" + header.getValue());
                }
                System.out.println("statusCode:" + statusCode + "\n");
            }
            if (statusCode == 206 || (statusCode == 200 && !taskBean.isIsRange())) {
                InputStream inputStream = response.getEntity().getContent();
                // BufferedInputStream bis = new BufferedInputStream(is, temp.length);

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
                 channel.close();
                 */

                outputStream.seek(taskBean.getSectionsOffset()[sectionID]);
                int count = 0;
                byte[] buffer = new byte[10 * 1024];

                do {
                    count = inputStream.read(buffer, 0, buffer.length);
                    if (count <= 0) {
                        break;
                    }
                    outputStream.write(buffer, 0, count);
                    offset += count;
                    downloadLen += count;
                    taskBean.getSectionsOffset()[sectionID] = offset;
                    //���������¼�
                    fireAfterPerDown(new DownloadThreadEvent(this, count));
                } while (notStop);

                outputStream.close();
                if (notStop) {
                    //response.getEntity().consumeContent();
                    EntityUtils.consume(response.getEntity());
                }
            }
            httpGet.abort();
            httpClient.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        } finally {

            fireDownCompleted(new DownloadThreadEvent(this, taskBean.getSectionsEnd()[sectionID]));
            if (HttpClient4DownloadTask.getDebug()) {
                System.out.println("Task " + sectionID + ": " + startSite + "-" + taskBean.getSectionsEnd()[sectionID] + " download bytes: " + downloadLen);
            }
            // httpClient.getConnectionManager().shutdown();

            return "Download thread " + sectionID + " over";
        }
    }
}
