/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.download.stable;

/**
 *
 * @author mxli
 */
public class DownloadThreadEvent {

    private DownloadThread target = null;
    private long count = 0;

    public DownloadThreadEvent(DownloadThread target, long count) {
        this.target = target;
        this.count = count;
    }

    /**
    
     * @return
     */
    public DownloadThread getTarget() {
        return target;
    }

    /**
    
     * @return
     */
    public long getCount() {
        return count;
    }
}
