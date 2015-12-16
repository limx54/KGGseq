/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.download.stable;

/**
 *
 * @author mxli
 */
public class DownloadTaskEvent {

    private long receivedCount = 0;
    private long totalCount = 0;
    private long availableCount = 0;
    private String realTimeSpeed = "";
    private String globalSpeed = "";
    private boolean complete = false;


    public DownloadTaskEvent(long ReceivedCount, long avaibleCount, long totalCount,
            String realTimeSpeed, String globalSpeed) {
        this(totalCount, avaibleCount, totalCount, globalSpeed, globalSpeed, false);
    }

    public DownloadTaskEvent(long ReceivedCount, long avaibleCount, long totalCount,
            String realTimeSpeed, String globalSpeed, boolean complete) {
        this.receivedCount = ReceivedCount;
        this.availableCount = avaibleCount;
        this.totalCount = totalCount;
        this.realTimeSpeed = realTimeSpeed;
        this.globalSpeed = globalSpeed;
        this.complete = complete;
    }

    public long getAvailableCount() {
        return availableCount;
    }

    public long getTotalDownloadedCount(){
        return (availableCount+receivedCount);
    }
    /**
     * ��ȡ�ļ��ѽ��մ�С(�ֽ���)
     * @return
     */
    public long getReceivedCount() {
        return receivedCount;
    }

    /**
     * ��ȡ�ļ��ܴ�С(�ֽ���)
     * @return
     */
    public long getTotalCount() {
        return totalCount;
    }

    /**
     * ��ȡʵʱ�ٶ�
     * @return
     */
    public String getRealTimeSpeed() {
        return realTimeSpeed;
    }

    /**
     * ��ȡȫ���ٶ�
     * @return
     */
    public String getGlobalSpeed() {
        return globalSpeed;
    }

    public boolean isComplete() {
        return complete;
    }
}
