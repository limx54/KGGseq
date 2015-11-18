/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.thread;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author mxli
 */
public class Task {

    protected List<TaskListener> listeners = Collections.synchronizedList(new LinkedList<TaskListener>());

    public void addTaskListener(TaskListener listener) {
        listeners.add(listener);
    }

    /*
    protected void fireAutoCallback(String infor) {
        if (listeners.isEmpty()) {
            return;
        }
        for (TaskListener listener : listeners) {
            listener.autoCallback(infor);
        }
    }
*/
    protected void fireTaskComplete() throws Exception {
        if (listeners.isEmpty()) {
            return;
        }
        for (TaskListener listener : listeners) {
            listener.taskCompleted();
        }
    }
}
