/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.thread;

/**
 *
 * @author mxli
 */
public interface TaskListener {

 //   public void autoCallback(String infor);

    public void taskCompleted() throws Exception;
}
