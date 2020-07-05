package Model;

import java.sql.*;
import java.util.ArrayList;

public class StudentStoreDatabase {
    public static final String ATTRIBUTE = "database";
    private static final String DB_NAME = "my_database";
    private Connection dbCon;

    public StudentStoreDatabase() {
        try {
            Class.forName("com.mysql.cj.jdbc.Driver");
            dbCon = DriverManager.getConnection(
                    "jdbc:mysql://localhost",
                    "root",
                    "my_password");
            Statement useDbStm = dbCon.createStatement();
            useDbStm.execute("USE " + DB_NAME);
        } catch (SQLException | ClassNotFoundException e) {
            throw new RuntimeException("");
        }

    }

    public ArrayList<Product> getProductList() throws SQLException {
        Statement searchStatement = dbCon.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE, ResultSet.CONCUR_READ_ONLY);
        String query = "SELECT * FROM products";
        ResultSet rs = searchStatement.executeQuery(query);
        ArrayList<Product> products = new ArrayList<>();
        while (rs.next()) {
            Product p = getProduct(rs);
            products.add(p);
        }
        return products;
    }

    public Product findProduct(String productId) {
        try {
            Statement searchStatement = dbCon.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE, ResultSet.CONCUR_READ_ONLY);
            String query = "SELECT * FROM products WHERE productid = \"" + productId + "\"";
            ResultSet rs = searchStatement.executeQuery(query);
            rs.next();
            return getProduct(rs);
        } catch (SQLException throwables) {
            throwables.printStackTrace();
        }
        return null;
    }

    private Product getProduct(ResultSet rs) throws SQLException {
        String id = rs.getString("productid");
        String name = rs.getString("name");
        String imageFile = rs.getString("imagefile");
        double price = rs.getDouble("price");
        Product p = new Product(id, name, imageFile, price);
        return p;
    }

}
