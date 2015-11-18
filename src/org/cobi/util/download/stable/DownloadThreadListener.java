/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.cobi.util.download.stable;

/**
 *
 * @author mxli
 */
public interface DownloadThreadListener {
   /**
     
     * @param event
     */
    public void afterPerDown( DownloadThreadEvent event);

    /**
     
     * @param event
     */
    public void downCompleted( DownloadThreadEvent event);
}
