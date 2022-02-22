import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Semaphore;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

public class WebFrame extends JFrame {

    private DefaultTableModel model;
    private JPanel panel;

    private JButton singleThreadFetchButton;
    private JButton concurrentFetchButton;
    private JButton stopButton;

    private JTextField threadNumLimiter;

    private JLabel running;
    private JLabel completed;
    private JLabel elapsedTimeCount;

    private JProgressBar progressBar;

    private long startTime;

    private List<String> urls;
    private List<Thread> threads;

    private ThreadLauncher launcher;

    private Semaphore limiter;
    private Lock interruptionLock;


    private static final int TABLE_STATUS_INDEX = 1;


    public WebFrame(List<String> urls) {
        super("WebLoader");
        this.urls = urls;

        panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        addTable();

        addFetchButtons();
        addTextField();
        addLabels();
        addProgressBar();
        addStopButton();
        add(panel);

        startTime = 0;

        threads = new ArrayList<>();
        interruptionLock = new ReentrantLock();

        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        pack();
        setVisible(true);
    }

    /////////////////////////////////////////////
    //GUI INITIALIZATION FUNCTIONS
    /////////////////////////////////////////////

    private void addTable() {
        model = new DefaultTableModel(new String[]{"url", "status"}, 0);
        JTable table = new JTable(model);
        table.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
        JScrollPane scrollpane = new JScrollPane(table);
        scrollpane.setPreferredSize(new Dimension(600, 300));

        for (String url : urls)
            model.addRow(new String[]{url, ""});

        table.setShowGrid(false);
        panel.add(scrollpane);
        panel.add(Box.createRigidArea(new Dimension(0, 5)));
    }

    private void addFetchButtons() {
        singleThreadFetchButton = new JButton("Single Thread Fetch");
        singleThreadFetchButton.addActionListener(e -> {
            int limit = 1;
            fetch(limit);
        });
        panel.add(singleThreadFetchButton);
        panel.add(Box.createRigidArea(new Dimension(0, 5)));

        concurrentFetchButton = new JButton("Concurrent Fetch");
        concurrentFetchButton.addActionListener(e -> {
            int limit = Integer.parseInt(threadNumLimiter.getText());
            fetch(limit);
        });
        panel.add(concurrentFetchButton);
        panel.add(Box.createRigidArea(new Dimension(0, 8)));
    }

    private void addTextField() {
        threadNumLimiter = new JTextField();
        threadNumLimiter.setMaximumSize(new Dimension(50, 30));
        panel.add(threadNumLimiter);
        panel.add(Box.createRigidArea(new Dimension(0, 1)));
    }

    private void addLabels() {
        running = new JLabel();
        running.setText("Running:0");
        panel.add(running);

        completed = new JLabel();
        completed.setText("Completed:0");
        panel.add(completed);

        elapsedTimeCount = new JLabel();
        elapsedTimeCount.setText("Elapsed:0.0");
        panel.add(elapsedTimeCount);
        panel.add(Box.createRigidArea(new Dimension(0, 10)));
    }


    private void addProgressBar() {
        progressBar = new JProgressBar();
        panel.add(progressBar);
        panel.add(Box.createRigidArea(new Dimension(0, 10)));
    }

    private void addStopButton() {
        stopButton = new JButton("Stop");
        stopButton.setEnabled(false);
        stopButton.addActionListener(e -> {
            interruptionLock.lock();
            launcher.interrupt();
            for (Thread thread : threads) thread.interrupt();
            interruptionLock.unlock();
            SwingUtilities.invokeLater(() -> stopButton.setEnabled(false));
        });
        panel.add(stopButton);
        panel.add(Box.createRigidArea(new Dimension(0, 5)));
    }


    /////////////////////////////////////////////
    //GUI HELPER FUNCTIONS
    /////////////////////////////////////////////


    private void fetch(int limit) {
        changeStateToRunning();
        launcher = new ThreadLauncher(limit);
        launcher.start();
    }


    private void changeStateToRunning() {
        resetFrame();
        startTime = System.currentTimeMillis();
        SwingUtilities.invokeLater(() -> {
            stopButton.setEnabled(true);
            concurrentFetchButton.setEnabled(false);
            singleThreadFetchButton.setEnabled(false);
            progressBar.setMaximum(urls.size());
        });
    }

    private void resetFrame() {
        startTime = 0;
        SwingUtilities.invokeLater(() -> {
            running.setText("Running:0");
            completed.setText("Completed:0");
            progressBar.setValue(0);
        });
    }


    private void changeStateToReady() {
        double timeElapsed = (double) (System.currentTimeMillis() - startTime) / 1000;
        SwingUtilities.invokeLater(() -> {
            elapsedTimeCount.setText("Elapsed:" + timeElapsed);
            stopButton.setEnabled(false);
            concurrentFetchButton.setEnabled(true);
            singleThreadFetchButton.setEnabled(true);
        });
    }

    /////////////////////////////////////////////
    //PROTECTED FUNCTIONS FOR WEBWORKERS
    /////////////////////////////////////////////

    protected void incrementActiveThreads() {
        synchronized (running) {
            SwingUtilities.invokeLater(() -> {
                int activeThreadCount = Integer.parseInt(running.getText().substring(8)) + 1; //length of "running:"
                running.setText("Running:" + activeThreadCount);
            });
        }
    }

    protected void decrementActiveThreads() {
        synchronized (running) {

            SwingUtilities.invokeLater(() -> {
                int activeThreadCount = Integer.parseInt(running.getText().substring(8)) - 1; //length of "running:"
                running.setText("Running:" + activeThreadCount);
            });
        }

        synchronized (completed) {
            SwingUtilities.invokeLater(() -> {
                int completedThreadCount = Integer.parseInt(completed.getText().substring(10)) + 1;//length of "completed:"
                completed.setText("Completed:" + completedThreadCount);
                progressBar.setValue(progressBar.getValue() + 1);
            });
            limiter.release();
        }
    }

    protected void updateStatus(String status, int row) {
        synchronized (model) {
            SwingUtilities.invokeLater(() -> model.setValueAt(status, row, TABLE_STATUS_INDEX));
        }
    }


    /////////////////////////////////////////////
    //THREADLAUNCHER CLASS
    /////////////////////////////////////////////


    private class ThreadLauncher extends Thread {
        private int threadNumLimit;

        public ThreadLauncher(int limit) {
            threadNumLimit = limit;
        }

        @Override
        public void run() {
            incrementActiveThreads();
            limiter = new Semaphore(threadNumLimit);
            for (int i = 0; i < urls.size(); i++) {
                try {
                    limiter.acquire();
                } catch (InterruptedException e) {
                    break;
                }
                interruptionLock.lock();
                if (isInterrupted()) break;
                WebWorker worker = new WebWorker(urls.get(i), i, WebFrame.this);
                threads.add(worker);
                worker.start();
                interruptionLock.unlock();
            }
            endThreadLauncher();
        }

        private void threadLauncherInterrupted() {
            for (int i = 0; i < threadNumLimit - 1; i++) {
                try {
                    limiter.acquire();
                } catch (InterruptedException ignored) {
                }
            }
        }

        private void endThreadLauncher() {
            for (Thread thread : threads) {
                try {
                    thread.join();
                } catch (InterruptedException e) {
                    threadLauncherInterrupted();
                }
            }
            updateGUI();
        }

        private void updateGUI() {
            SwingUtilities.invokeLater(() -> {
                int activeThreadCount = Integer.parseInt(running.getText().substring(8)) - 1; //length of "running:"
                running.setText("Running:" + activeThreadCount);
            });
            changeStateToReady();
        }
    }
}
