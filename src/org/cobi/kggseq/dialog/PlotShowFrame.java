/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * PlotShowFrame.java
 *
 * Created on 2011-5-9, 9:05:34
 */
package org.cobi.kggseq.dialog;

import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.text.Document;

/**
 *
 * @author mxli
 */
public class PlotShowFrame extends javax.swing.JFrame {

    /** Creates new form PlotShowFrame */
    public PlotShowFrame() {
        JFrame.setDefaultLookAndFeelDecorated(true);
        initComponents();
       
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jScrollPane1 = new javax.swing.JScrollPane();
        breifTextPane = new javax.swing.JTextPane();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);

        breifTextPane.setEditable(false);
        jScrollPane1.setViewportView(breifTextPane);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 635, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 590, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        java.awt.EventQueue.invokeLater(new Runnable() {

            @Override
            public void run() {
                new PlotShowFrame().setVisible(true);
            }
        });
    }

    public void insertImage2PlottingPane(File path) {
        try {
            //  URL pageUrl = new URL("file:///" + path.getCanonicalPath());
            //  String infor = "<img src=\"" + pageUrl + "\"/>";
            this.setTitle(path.getCanonicalPath());
            BufferedImage image = javax.imageio.ImageIO.read(path);
            this.setSize(image.getWidth() + 40, image.getHeight() + 40);
            Dimension displaySize = Toolkit.getDefaultToolkit().getScreenSize();
            Dimension frameSize = this.getSize();
            if (frameSize.width > displaySize.width) {
                frameSize.width = displaySize.width;
            }
            if (frameSize.height > displaySize.height) {
                frameSize.height = displaySize.height;
            }
            this.setLocation((displaySize.width - frameSize.width) / 2, (displaySize.height - frameSize.height) / 2);
            breifTextPane.setSelectionStart(breifTextPane.getDocument().getLength());
            breifTextPane.insertIcon(new ImageIcon(image));
            Document docs = breifTextPane.getDocument();
            docs.insertString(docs.getLength(), "\n Source: "+path.getCanonicalPath()+"\n",null);
            breifTextPane.setSelectionStart(breifTextPane.getDocument().getLength());
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTextPane breifTextPane;
    private javax.swing.JScrollPane jScrollPane1;
    // End of variables declaration//GEN-END:variables
}