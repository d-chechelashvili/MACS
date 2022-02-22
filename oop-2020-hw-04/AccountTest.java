import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class AccountTest {


    @Test
    void testAccount() {
        Bank b = new Bank(5);
        Account acc = new Account(b, 0, 1000);
        assertEquals(1000, acc.getBalance());
    }

    @Test
    void testMakeTransaction() {
        Bank b = new Bank(5);
        Transaction tr = new Transaction(0, 1, 500);
        Account.makeTransaction(b, tr);
        assertEquals(500, b.getAccount(0).getBalance());
        assertEquals(1500, b.getAccount(1).getBalance());
    }

    @Test
    void changeBalance() {
        Bank b = new Bank(5);
        Account acc = new Account(b, 0, 1000);
        acc.changeBalance(500);
        assertEquals(1500, acc.getBalance());
    }

    @Test
    void testToString() {
        Bank b = new Bank(5);
        Account acc = new Account(b, 0, 1000);
        acc.changeBalance(500);
        assertEquals("acct:0 bal:1500 trans:1", acc.toString());
    }
}