import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import static org.junit.jupiter.api.Assertions.*;

class BankTest {

    @Test
    void testConstructor() {
        int numWorkers = 5;
        Bank b = new Bank(numWorkers);
    }

    @Test
    void testProcessFile() throws IOException, InterruptedException {
        int numWorkers = 5;
        Bank b = new Bank(numWorkers);
        b.processFile("5k.txt");
        for (int i = 0; i < Bank.ACCOUNTS; i++) {
            assertEquals(b.getAccount(i).getBalance(), Bank.ACCOUNT_INITIAL_BALANCE);
        }

    }

    @Test
    void testBank() throws IOException, InterruptedException {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(os);
        System.setOut(ps);
        int numWorkers = 5;
        String[] args = new String[2];
        args[0] = "small.txt";
        args[1] = String.valueOf(numWorkers);
        Bank.main(args);
        assertEquals(os.toString(),
                args[0] + " " + args[1] + "\nacct:0 bal:999 trans:1\n" +
                "acct:1 bal:1001 trans:1\n" +
                "acct:2 bal:999 trans:1\n" +
                "acct:3 bal:1001 trans:1\n" +
                "acct:4 bal:999 trans:1\n" +
                "acct:5 bal:1001 trans:1\n" +
                "acct:6 bal:999 trans:1\n" +
                "acct:7 bal:1001 trans:1\n" +
                "acct:8 bal:999 trans:1\n" +
                "acct:9 bal:1001 trans:1\n" +
                "acct:10 bal:999 trans:1\n" +
                "acct:11 bal:1001 trans:1\n" +
                "acct:12 bal:999 trans:1\n" +
                "acct:13 bal:1001 trans:1\n" +
                "acct:14 bal:999 trans:1\n" +
                "acct:15 bal:1001 trans:1\n" +
                "acct:16 bal:999 trans:1\n" +
                "acct:17 bal:1001 trans:1\n" +
                "acct:18 bal:999 trans:1\n" +
                "acct:19 bal:1001 trans:1\n");
    }

    @Test
    void testBank2() throws IOException, InterruptedException {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(os);
        System.setOut(ps);
        int numWorkers = 5;
        String[] args = new String[0];
        Bank.main(args);
        assertEquals(os.toString(),"Args: transaction-file [num-workers [limit]]\n");
    }
    @Test
    void testBank3() throws IOException, InterruptedException {
        String[] args = new String[1];
        args[0] = "small.txt";
        assertThrows(ArrayIndexOutOfBoundsException.class,()->{
            Bank.main(args);
        } );
    }

}