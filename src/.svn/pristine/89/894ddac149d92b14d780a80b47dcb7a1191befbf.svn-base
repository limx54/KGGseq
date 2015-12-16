/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.file;

/**
 *
 * @author mxli
 */
public class ReadFileThread extends Thread {

    private ReaderFileListener processPoiDataListeners;
    private String filePath;
    private long start;
    private long end;

    public ReadFileThread(ReaderFileListener processPoiDataListeners, long start, long end, String file) {
        this.setName(this.getName() + "-ReadFileThread");
        this.start = start;
        this.end = end;
        this.filePath = file;
        this.processPoiDataListeners = processPoiDataListeners;
    }

    @Override
    public void run() {
        //ReadFile readFile = new ReadFile();  
        ReadGZipFile readFile = new ReadGZipFile();
        readFile.setReaderListener(processPoiDataListeners);
        readFile.setEncode(processPoiDataListeners.getEncode());
//        readFile.addObserver();  
        try {
            readFile.readFileByLine(filePath, start, end + 1);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
