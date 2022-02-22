import java.util.Objects;

// Transaction.java
/*
 (provided code)
 Transaction is just a dumb struct to hold
 one transaction. Supports toString.
*/
public class Transaction {

    private int from;
    private int to;
    private int amount;

    public Transaction(int from, int to, int amount) {
        this.from = from;
        this.to = to;
        this.amount = amount;
    }

    public int getSourceId() {
        return from;
    }

    public int getDestinationId() {
        return to;
    }

    public int getAmount() {
        return amount;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Transaction that = (Transaction) o;
        return from == that.from &&
                to == that.to &&
                amount == that.amount;
    }

    public String toString() {
        return ("from:" + from + " to:" + to + " amt:" + amount);
    }
}
