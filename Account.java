// Account.java

/*
 Simple, thread-safe Account class encapsulates
 a balance and a transaction count.
*/
public class Account {
    private int id;
    private int balance;
    private int transactions;

    // It may work out to be handy for the account to
    // have a pointer to its Bank.
    // (a suggestion, not a requirement)
    private Bank bank;

    public Account(Bank bank, int id, int balance) {
        this.bank = bank;
        this.id = id;
        this.balance = balance;
        transactions = 0;
    }

    public int getBalance() {
        return balance;
    }

    //makes transaction
    public static void makeTransaction(Bank bank, Transaction tr) {
        Account from = bank.getAccount(tr.getSourceId());
        Account to = bank.getAccount(tr.getDestinationId());
        int amount = tr.getAmount();
        from.changeBalance(-amount);
        to.changeBalance(amount);
    }

    public synchronized void changeBalance(int amount) {
        balance += amount;
        transactions++;
    }

    @Override
    public String toString() {
        return "acct:" + id + " bal:" + balance + " trans:" + transactions;
    }
}
