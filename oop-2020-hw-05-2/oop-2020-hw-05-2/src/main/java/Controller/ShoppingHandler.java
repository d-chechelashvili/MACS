package Controller;

import Model.ShoppingCart;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.util.Enumeration;

public class ShoppingHandler extends HttpServlet {

    @Override
    protected void doPost(HttpServletRequest req, HttpServletResponse resp) throws ServletException, IOException {
        String productID = req.getParameter("productID");
        ShoppingCart cart = (ShoppingCart) req.getSession().getAttribute(ShoppingCart.ATTRIBUTE);
        if (productID != null) {
            cart.add(productID);
        } else {
            Enumeration<String> items = req.getParameterNames();
            while (items.hasMoreElements()) {
                String id = items.nextElement();
                int count = Integer.parseInt(req.getParameter(id));
                cart.setCount(id, count);
            }
        }
        req.getRequestDispatcher("/WEB-INF/shopping cart.jsp").forward(req, resp);
    }
}
