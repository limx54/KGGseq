/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.download.stable;

/**
 *
 * @author Miaoxin Li
 */
public class DownloadTaskBean {

    private String downURL;
    private String saveFile;
    private int sectionCount;
    private long contentLength;
    private long[] sectionsOffset;
    private long[] sectionsEnd;
    private boolean isRange = true;
    public static final int HEAD_SIZE = 4096;
    private long taskFileLen = 0;

    public DownloadTaskBean() {
    }



    public long getTaskFileLen() {
        taskFileLen = HEAD_SIZE + 8 * 2 * sectionCount;
        return taskFileLen;
    }

    public boolean isIsRange() {
        return isRange;
    }

    public void setIsRange(boolean isRange) {
        this.isRange = isRange;
    }

    public void setSectionsEnd(long[] sectionsEnd) {
        this.sectionsEnd = sectionsEnd;
    }

    public long getContentLength() {
        return contentLength;
    }

    public long[] getSectionsEnd() {
        return sectionsEnd;
    }

    public void setContentLength(long contentLength) {
        this.contentLength = contentLength;
    }

    public String getDownURL() {
        return downURL;
    }

    public void setDownURL(String downURL) {
        this.downURL = downURL;
    }

    public String getSaveFile() {
        return saveFile;
    }

    public void setSaveFile(String saveFile) {
        this.saveFile = saveFile;
    }

    public int getSectionCount() {
        return sectionCount;
    }

    public void setSectionCount(int sectionCount) {
        this.sectionCount = sectionCount;
    }

    public long[] getSectionsOffset() {
        return sectionsOffset;
    }

    public void setSectionsOffset(long[] sectionsOffset) {
        this.sectionsOffset = sectionsOffset;
    }
}
