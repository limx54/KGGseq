/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.download.stable;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author mxli
 */
public class DownloadThread {

    protected int sectionID;
    protected List<DownloadThreadListener> listeners = Collections.synchronizedList(new LinkedList<DownloadThreadListener>());
    protected DownloadTaskBean taskBean = null;
    protected boolean notStop = true;

    public void  stop() {
        notStop = false;
    }

    public DownloadThread(int sectionID, DownloadTaskBean taskBean) {
        this.sectionID = sectionID;
        this.taskBean = taskBean;
    }

    protected void fireAfterPerDown(DownloadThreadEvent event) {
        if (listeners.isEmpty()) {
            return;
        }
        for (DownloadThreadListener listener : listeners) {
            listener.afterPerDown(event);
        }
    }

    protected void fireDownCompleted(DownloadThreadEvent event) {
        if (listeners.isEmpty()) {
            return;
        }

        for (DownloadThreadListener listener : listeners) {
            listener.downCompleted(event);
        }
    }

    protected void addDownloadListener(DownloadThreadListener listener) {
        listeners.add(listener);
    }
}
