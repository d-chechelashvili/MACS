import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class TransactionTest {


    @Test
    void testTransaction(){
        Transaction tr = new Transaction(0, 1, 500);
        assertEquals(0, tr.getSourceId());
        assertEquals(1, tr.getDestinationId());
        assertEquals(500, tr.getAmount());
    }

    @Test
    void testEquals() {
        Transaction tr = new Transaction(0, 1, 500);
        Transaction tr1 = new Transaction(0, 1, 500);
        assertTrue(tr1.equals(tr));
    }

    @Test
    void testToString() {
        Transaction tr = new Transaction(0, 1, 500);
        assertEquals("from:0 to:1 amt:500", tr.toString());
    }
}