// Cracker.java
/*
 Generates SHA hashes of short strings in parallel.
*/

import java.security.*;
import java.util.Arrays;
import java.util.concurrent.CountDownLatch;

public class Cracker {
    // Array of chars used to produce strings
    public static final char[] CHARS = "abcdefghijklmnopqrstuvwxyz0123456789.,-!".toCharArray();

    private int numWorkers;
    private CountDownLatch latch;
    private byte[] target;


    public Cracker(int numWorkers, String target) {
        this.numWorkers = numWorkers;
        latch = new CountDownLatch(numWorkers);
        this.target = hexToArray(target);
    }


    public void crack(int maxLength) throws NoSuchAlgorithmException, InterruptedException {
        long initialTime = System.currentTimeMillis();
        startWorkers(maxLength);

        latch.await();

        long timeElapsed = System.currentTimeMillis() - initialTime;
        System.out.println("all done");
        System.out.println("Time Elapsed: " + timeElapsed);

    }

    private void startWorkers(int maxLength) throws NoSuchAlgorithmException {
        boolean inbounds = true;
        int startIndex = 0;
        int endIndex = -1;
        while (inbounds) {
            startIndex = endIndex + 1;
            endIndex += CHARS.length / numWorkers;
            if (endIndex >= CHARS.length - 1) {
                endIndex = CHARS.length - 1;
                inbounds = false;
            }
            Worker worker = new Worker(startIndex, endIndex, maxLength);
            worker.start();
        }
    }

    /*
     Given a byte[] array, produces a hex String,
     such as "234a6f". with 2 chars for each byte in the array.
     (provided code)
    */
    public static String hexToString(byte[] bytes) {
        StringBuffer buff = new StringBuffer();
        for (int i = 0; i < bytes.length; i++) {
            int val = bytes[i];
            val = val & 0xff;  // remove higher bits, sign
            if (val < 16) buff.append('0'); // leading 0
            buff.append(Integer.toString(val, 16));
        }
        return buff.toString();
    }

    /*
     Given a string of hex byte values such as "24a26f", creates
     a byte[] array of those values, one byte value -128..127
     for each 2 chars.
     (provided code)
    */
    public static byte[] hexToArray(String hex) {
        byte[] result = new byte[hex.length() / 2];
        for (int i = 0; i < hex.length(); i += 2) {
            result[i / 2] = (byte) Integer.parseInt(hex.substring(i, i + 2), 16);
        }
        return result;
    }


    public static void main(String[] args) throws NoSuchAlgorithmException, InterruptedException {
        if (args.length < 1) {
            System.out.println("Args: target length [workers]");
            return;
        }
        // args: targ len [num]
        String targ = args[0];

        if (args.length > 2) {
            int len = Integer.parseInt(args[1]);
            int num = Integer.parseInt(args[2]);
            Cracker cracker = new Cracker(num, targ);
            cracker.crack(len);
        } else generateHash(targ);


        // a! 34800e15707fae815d7c90d49de44aca97e2d759
        // xyz 66b27417d37e024c46526c2f6d358a754fc552f3

    }

    private static void generateHash(String targ) throws NoSuchAlgorithmException {

        MessageDigest md = MessageDigest.getInstance("SHA-1");
        md.update(targ.getBytes());
        String hash = hexToString(md.digest());
        System.out.println(hash);

    }


    private class Worker extends Thread {

        private int rangeStart;
        private int rangeEnd;
        private int maxLength;
        private MessageDigest md;

        public Worker(int rangeStart, int rangeEnd, int maxLength) throws NoSuchAlgorithmException {
            this.rangeStart = rangeStart;
            this.rangeEnd = rangeEnd;
            this.maxLength = maxLength;
            md = MessageDigest.getInstance("SHA-1");
        }

        @Override
        public void run() {
            for (int i = rangeStart; i <= rangeEnd; i++) {
                StringBuilder soFar = new StringBuilder(String.valueOf(CHARS[i]));
                generatePassword(soFar);
            }
            latch.countDown();
        }

        private void generatePassword(StringBuilder soFar) {
            if (soFar.length() > maxLength) return;

            md.update(soFar.toString().getBytes());
            if (Arrays.equals(md.digest(), target)) System.out.println(soFar);

            for (char ch : CHARS) {
                soFar.append(ch);
                generatePassword(soFar);
                soFar.deleteCharAt(soFar.length() - 1);
            }
        }
    }
}
