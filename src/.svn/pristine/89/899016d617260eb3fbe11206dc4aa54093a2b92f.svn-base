/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.net;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPReply;



/**
 *
 * @author JiangLi
 */
public class FTP {
    FTPClient ftp=null;
    String strIP;
    int intPort;
    String strUser;
    String strPassword;

    public FTP(String strIP, int intPort, String strUser, String strPassword) {
        this.strIP = strIP;
        this.intPort = intPort;
        this.strUser = strUser;
        this.strPassword = strPassword;
    }
    
    public boolean openFTP() throws IOException{
        if(ftp!=null && ftp.isConnected())  return true;
        ftp=new FTPClient();
        ftp.connect(this.strIP, this.intPort);
        ftp.login(this.strUser, this.strPassword);
        int intResponse=ftp.getReplyCode();
        if(!FTPReply.isPositiveCompletion(intResponse)){
            this.closeFTP();
            System.err.println("FTP visiting is not successful!");
            return false;
        }
        return true;
    }
    
    public void closeFTP() throws IOException{
        if(ftp!=null && ftp.isConnected())  ftp.disconnect();
    }
     
    public boolean downloadFTP(String strLink,File fleOutput) throws FileNotFoundException, IOException{      
        OutputStream os=new BufferedOutputStream(new FileOutputStream(fleOutput));
        boolean boolFlag=false;
        if(!ftp.isConnected())  return false;
        try{
            boolFlag=ftp.retrieveFile(strLink, os);
           // System.out.print(strLink+" ---> "+ftp.getReplyString());
        }catch(Exception e){
            e.printStackTrace();
        }finally{
            os.close();
            return boolFlag;
        }  
    }    
    
}
