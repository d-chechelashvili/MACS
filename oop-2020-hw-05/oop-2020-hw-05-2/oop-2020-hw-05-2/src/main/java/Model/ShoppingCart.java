package Model;

import java.util.HashMap;
import java.util.Set;

public class ShoppingCart {
    public static final String ATTRIBUTE = "shopping cart";
    private HashMap<String, Integer> cart;

    public ShoppingCart() {
        cart = new HashMap<>();
    }

    public void add(String productID) {
        cart.merge(productID, 1, Integer::sum);
    }

    public Set<String> getItems() {
        return cart.keySet();
    }

    public int getCount(String productID) {
        return cart.get(productID);
    }

    public void setCount(String productID, int count) {
        if (count <= 0) cart.remove(productID);
        else cart.replace(productID, count);
    }
}
