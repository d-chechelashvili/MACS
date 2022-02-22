import org.junit.jupiter.api.Test;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.security.NoSuchAlgorithmException;

import static org.junit.jupiter.api.Assertions.*;

class CrackerTest {

    @Test
    void testCracker1() throws InterruptedException, NoSuchAlgorithmException {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(os);
        System.setOut(ps);
        String[] args = new String[0];
        Cracker.main(args);
        assertEquals(os.toString(), "Args: target length [workers]\n");
    }

    @Test
    void testCracker2() throws InterruptedException, NoSuchAlgorithmException {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(os);
        System.setOut(ps);
        String[] args = new String[1];
        args[0] = "a!";
        Cracker.main(args);
        assertEquals("34800e15707fae815d7c90d49de44aca97e2d759\n", os.toString());
    }

    @Test
    void testCracker3() throws InterruptedException, NoSuchAlgorithmException {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(os);
        System.setOut(ps);
        String[] args = new String[1];
        args[0] = "xyz";
        Cracker.main(args);
        assertEquals("66b27417d37e024c46526c2f6d358a754fc552f3\n", os.toString());
    }

    @Test
    void testCracker4() throws InterruptedException, NoSuchAlgorithmException {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(os);
        System.setOut(ps);
        String[] args = new String[3];
        args[0] = "66b27417d37e024c46526c2f6d358a754fc552f3";
        args[1] = "3";
        args[2] = "5";
        Cracker.main(args);
        assertEquals("xyz", os.toString().substring(0, os.toString().indexOf('\n') - 1));
    }

    @Test
    void testCracker5() throws InterruptedException, NoSuchAlgorithmException {
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(os);
        System.setOut(ps);
        String[] args = new String[3];
        args[0] = "34800e15707fae815d7c90d49de44aca97e2d759";
        args[1] = "3";
        args[2] = "4";
        Cracker.main(args);
        assertEquals("a!", os.toString().substring(0, os.toString().indexOf('\n') - 1));
    }
}