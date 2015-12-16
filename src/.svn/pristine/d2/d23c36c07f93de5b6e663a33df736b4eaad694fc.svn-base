/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.download.stable;


import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.border.EmptyBorder;
 

public class DownloadFrame extends JFrame {

   
    private JPanel contentPane;
    private JTextField txtURL;
    private JTextField txtLocal;
    private JLabel tSpeedl;
    private JLabel rtSpeedl;
    private JButton btnStart;
    private JProgressBar progressBar;
    private JLabel label_4;
    private JTextField txtThreadCount;
    private JCheckBox cbAutoNamed;

    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        try {
            UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
        } catch (Exception e2) {
            e2.printStackTrace();
        }
        EventQueue.invokeLater(new Runnable() {

            public void run() {
                try {
                    DownloadFrame frame = new DownloadFrame();
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }

    /**
     * Create the frame.
     */
    public DownloadFrame() {
        setTitle("\u591A\u7EBF\u7A0B\u4E0B\u8F7D\u6D4B\u8BD5\u7A97\u4F53--ZYWANG");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setBounds(100, 100, 450, 242);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        setContentPane(contentPane);
        contentPane.setLayout(null);

        btnStart = new JButton("\u5F00\u59CB\u4E0B\u8F7D");
        btnStart.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    btnStart.setEnabled(false);
                    progressBar.setValue(0);
                    rtSpeedl.setText("");
                    tSpeedl.setText("");
                    HttpClient4DownloadTask.setDebug(true);  
                    HttpClient4DownloadTask task = new HttpClient4DownloadTask(txtURL.getText(), txtLocal.getText(), Integer.valueOf(txtThreadCount.getText()));

                              
                    task.addTaskListener(new DownloadTaskListener() {

                        @Override
                        public void autoCallback(DownloadTaskEvent event) {
                            int progess = (int) (event.getReceivedCount() * 100.0 / event.getTotalCount());
                            progressBar.setValue(progess);
                            rtSpeedl.setText(event.getRealTimeSpeed());
                            tSpeedl.setText(event.getGlobalSpeed());
                            if (event.isComplete()) {
                                btnStart.setEnabled(true);
                            }
                            System.out.println(progess + " " + event.getRealTimeSpeed() + " " + event.getGlobalSpeed());
                        }
                        @Override
                          public void taskCompleted() throws Exception {
                              
                          }
                    });
                    try {
                        if (cbAutoNamed.isSelected()) {
                            File file = new File(txtLocal.getText());
                            if (file.isDirectory()) {
                                txtLocal.setText(new File(file, task.guessFileName()).getPath());
                            } else {
                                txtLocal.setText(new File(file.getParent(), task.guessFileName()).getPath());
                            }
                            task.setLocalPath(txtLocal.getText());
                        }

                        task.call();

                    } catch (Exception e1) {
                        btnStart.setEnabled(true);
                        JOptionPane.showMessageDialog(null, "����ʧ��\n" + e1.getMessage(), "����", JOptionPane.ERROR_MESSAGE);
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }
        });
        btnStart.setBounds(331, 165, 93, 23);
        contentPane.add(btnStart);

        JLabel label = new JLabel("\u4E0B\u8F7D\u76EE\u6807\uFF1A");
        label.setBounds(10, 13, 68, 15);
        contentPane.add(label);

        txtURL = new JTextField();
        txtURL.setText("http://localhost:8080/eu/abc.mp3");
        txtURL.setBounds(85, 10, 339, 21);
        contentPane.add(txtURL);
        txtURL.setColumns(10);

        JLabel label_1 = new JLabel("\u672C\u5730\u5730\u5740\uFF1A");
        label_1.setBounds(10, 44, 68, 15);
        contentPane.add(label_1);

        txtLocal = new JTextField();
        txtLocal.setText("d:\\");
        txtLocal.setBounds(85, 41, 250, 21);
        contentPane.add(txtLocal);
        txtLocal.setColumns(10);

        progressBar = new JProgressBar();
        progressBar.setStringPainted(true);
        progressBar.setBounds(10, 116, 414, 14);
        contentPane.add(progressBar);

        JLabel label_2 = new JLabel("\u5168\u5C40\u901F\u5EA6\uFF1A");
        label_2.setBounds(293, 140, 97, 15);
        contentPane.add(label_2);

        JLabel label_3 = new JLabel("\u5B9E\u65F6\u901F\u5EA6\uFF1A");
        label_3.setBounds(164, 140, 93, 15);
        contentPane.add(label_3);

        rtSpeedl = new JLabel("");
        rtSpeedl.setBounds(226, 140, 68, 15);
        contentPane.add(rtSpeedl);

        tSpeedl = new JLabel("");
        tSpeedl.setBounds(350, 140, 74, 15);
        contentPane.add(tSpeedl);

        label_4 = new JLabel("\u7EBF\u7A0B\u4E2A\u6570\uFF1A");
        label_4.setBounds(10, 75, 68, 15);
        contentPane.add(label_4);

        txtThreadCount = new JTextField();
        txtThreadCount.setHorizontalAlignment(SwingConstants.RIGHT);
        txtThreadCount.setText("5");
        txtThreadCount.setBounds(85, 72, 40, 21);
        contentPane.add(txtThreadCount);
        txtThreadCount.setColumns(10);

        cbAutoNamed = new JCheckBox("\u667A\u80FD\u547D\u540D");
        cbAutoNamed.setSelected(true);
        cbAutoNamed.setBounds(341, 40, 83, 23);
        contentPane.add(cbAutoNamed);
    }
}
