package model;

import java.util.HashMap;

public class AccountManager {
    public static final String ATTRIBUTE = "acc_manager";
    private HashMap<String, String> data;

    public AccountManager(){
        data = new HashMap<>();
        data.put("Patrick", "1234");
        data.put("Molly", "FloPup");
    }

    public boolean accountExists(String username){
        return data.containsKey(username);
    }

    public boolean validPassword(String username, String password){
        return data.containsKey(username) && data.get(username).equals(password);
    }

    public void createAccount(String username, String password){
        data.putIfAbsent(username, password);
    }
}
