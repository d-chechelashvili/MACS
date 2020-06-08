import java.io.*;
import java.net.*;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.CountDownLatch;
import javax.swing.*;

public class WebWorker extends Thread {

    private String urlString;
    private int rowNum;
    private WebFrame webFrame;
    private CountDownLatch latch;

    public WebWorker(String url, int rowNumber, WebFrame wFrame) {
        urlString = url;
        rowNum = rowNumber;
        webFrame = wFrame;
    }

    @Override
    public void run() {
        webFrame.incrementActiveThreads();
        download();
        webFrame.decrementActiveThreads();
    }


    //This is the core web/download i/o code...
    private void download() {
        long startTime = System.currentTimeMillis();
        InputStream input = null;
        StringBuilder contents = null;
        try {
            URL url = new URL(urlString);
            URLConnection connection = url.openConnection();

            // Set connect() to throw an IOException
            // if connection does not succeed in this many msecs.
            connection.setConnectTimeout(5000);

            connection.connect();
            input = connection.getInputStream();

            BufferedReader reader = new BufferedReader(new InputStreamReader(input));

            char[] array = new char[1000];
            int len;
            contents = new StringBuilder(1000);
            while ((len = reader.read(array, 0, array.length)) > 0) {
                if(isInterrupted()) {
                    updateGUIFail(true);
                    return;
                }
                contents.append(array, 0, len);
                Thread.sleep(100);
            }

            // Successful download if we get here
            long timeElapsed = System.currentTimeMillis() - startTime;
            updateGUISuccess(timeElapsed, contents.length());

        }
        // Otherwise control jumps to a catch...
        catch(IOException ignored) {
            updateGUIFail(false);
        }
        catch(InterruptedException exception) {
            updateGUIFail(true);
        }

        // "finally" clause, to close the input stream
        // in any case
        finally {
            try {
                if (input != null) input.close();

            } catch (IOException ignored) {
            }
        }
    }

    private void updateGUISuccess(long timeElapsed, int bytes) {
        String completionTime = new SimpleDateFormat("HH:mm:ss").format(new Date(System.currentTimeMillis()));
        SimpleDateFormat s = new SimpleDateFormat();
        String status = completionTime + "  " + timeElapsed + "ms  " + bytes + "bytes";
        webFrame.updateStatus(status, rowNum);
    }

    private void updateGUIFail(boolean interrupted) {
        String status;
        if (interrupted) status = "interrupted";
        else status = "err";
        webFrame.updateStatus(status, rowNum);
    }
}
