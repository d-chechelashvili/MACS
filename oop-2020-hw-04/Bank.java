// Bank.java

/*
 Creates a bunch of accounts and uses threads
 to post transactions to the accounts concurrently.
*/

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

public class Bank {
    public static final int ACCOUNTS = 20;     // number of accounts
    public static final int ACCOUNT_INITIAL_BALANCE = 1000;
    private final Transaction nullTrans = new Transaction(-1, 0, 0);
    private CountDownLatch latch;
    private List<Account> accounts;
    private BlockingQueue<Transaction> transactionQueue;
    private int numWorkers;


    public Bank(int numWorkers) {
        this.numWorkers = numWorkers;
        latch = new CountDownLatch(numWorkers);
        transactionQueue = new ArrayBlockingQueue<Transaction>(numWorkers);
        initAccounts();
    }


    private void initAccounts() {
        accounts = new ArrayList<Account>();
        for (int i = 0; i < ACCOUNTS; i++) {
            Account acc = new Account(this, i, ACCOUNT_INITIAL_BALANCE);
            accounts.add(acc);
        }
    }

    private void startWorkers() {
        for (int i = 0; i < numWorkers; i++) {
            Worker worker = new Worker();
            worker.start();
        }
    }


    /*
     Reads transaction data (from/to/amt) from a file for processing.
     (provided code)
     */
    private void readFile(String file) throws IOException, InterruptedException {
        BufferedReader reader = new BufferedReader(new FileReader(file));

        // Use stream tokenizer to get successive words from file
        StreamTokenizer tokenizer = new StreamTokenizer(reader);

        while (true) {
            int read = tokenizer.nextToken();
            if (read == StreamTokenizer.TT_EOF) break;  // detect EOF
            int from = (int) tokenizer.nval;

            tokenizer.nextToken();
            int to = (int) tokenizer.nval;

            tokenizer.nextToken();
            int amount = (int) tokenizer.nval;

            Transaction tr = new Transaction(from, to, amount);
            transactionQueue.put(tr);
        }
    }

    /*
     Processes one file of transaction data
     -fork off workers
     -read file into the buffer
     -wait for the workers to finish
    */
    public void processFile(String file) throws InterruptedException, IOException {
        startWorkers();
        readFile(file);

        for (int i = 0; i < numWorkers; i++) {
            transactionQueue.put(nullTrans);
        }
        latch.await();
        printInfo();
    }

    private void printInfo() {
        for (Account account : accounts) {
            System.out.println(account);
        }
    }

    protected Account getAccount(int id) {
        return accounts.get(id);
    }


    /*
     Looks at commandline args and calls Bank processing.
    */
    public static void main(String[] args) throws InterruptedException, IOException {
        // deal with command-lines args
        if (args.length == 0) {
            System.out.println("Args: transaction-file [num-workers [limit]]");
            return;
        }

        String file = args[0];

        int numWorkers = 1;
        if (args.length >= 2) {
            numWorkers = Integer.parseInt(args[1]);
        }
        System.out.println(args[0] + " " + args[1]);

        Bank bank = new Bank(numWorkers);
        bank.processFile(file);
    }


    private class Worker extends Thread {
        @Override
        public void run() {
            while (true) {
                try {
                    Transaction tr = transactionQueue.take();
                    if (tr.equals(nullTrans)) break;
                    Account.makeTransaction(Bank.this, tr);
                    latch.countDown();
                } catch (InterruptedException ignored) {}
            }
        }
    }
}