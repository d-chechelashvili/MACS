// JCount.java

/*
 Basic GUI/Threading exercise.
*/

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class JCount extends JPanel {

    private JTextField destinationNum;
    private JLabel currCount;
    private JButton startButton;
    private JButton stopButton;
    private WorkerThread worker;
    private int count;

    private static final int INITIAL_COUNT = 0;

    public JCount() {
        // Set the JCount to use Box layout
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        worker = null;
        count = INITIAL_COUNT;

        destinationNum = new JTextField();
        add(destinationNum);

        currCount = new JLabel("" + count);
        add(currCount);

        startButton = new JButton("Start");
        startButton.addActionListener(e -> start());
        add(startButton);

        stopButton = new JButton("Stop");
        stopButton.addActionListener(e -> stop());
        add(stopButton);

        add(Box.createRigidArea(new Dimension(0, 40)));
    }

    private void stop() {
        if (worker != null) worker.interrupt();
    }

    private void start() {
        if (worker != null) worker.interrupt();
        count = INITIAL_COUNT;
        currCount.setText("" + count);
        startWorker();
    }

    private void startWorker() {
        int targ = Integer.parseInt(destinationNum.getText());
        worker = new WorkerThread(targ);
        worker.start();
    }


    private class WorkerThread extends Thread {

        private int targ;

        public WorkerThread(int targ) {
            this.targ = targ;
        }

        @Override
        public void run() {
            int count = 0;
            while (count < targ) {
                if (isInterrupted()) break;

                if (count % 10000 == 0) {
                    try {
                        sleep(100);
                    } catch (InterruptedException e) {
                        break;
                    }
                    updateDisplay(count);
                }
                count++;
            }
            currCount.setText("Finished Working");
        }

        private void updateDisplay(int count) {
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    currCount.setText("" + count);
                }
            });
        }
    }


    static public void main(String[] args) {
        // Creates a frame with 4 JCounts in it.
        // (provided)

        JFrame frame = new JFrame("The Count");
        frame.setLayout(new BoxLayout(frame.getContentPane(), BoxLayout.Y_AXIS));

        frame.add(new JCount());
        frame.add(new JCount());
        frame.add(new JCount());
        frame.add(new JCount());

        frame.pack();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }

}

