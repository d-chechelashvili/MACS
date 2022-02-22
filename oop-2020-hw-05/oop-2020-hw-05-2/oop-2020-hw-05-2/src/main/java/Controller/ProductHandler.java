package Controller;

import Model.Product;
import Model.StudentStoreDatabase;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;

public class ProductHandler extends HttpServlet {
    @Override
    protected void doGet(HttpServletRequest req, HttpServletResponse resp) throws ServletException, IOException {
        String productId = (String) req.getParameter("productId");
        StudentStoreDatabase db = (StudentStoreDatabase) getServletContext().getAttribute(StudentStoreDatabase.ATTRIBUTE);
        Product p = db.findProduct(productId);
        req.setAttribute("productName", p.getName());
        req.setAttribute("imageFile", p.getImageFile());
        req.setAttribute("price", "$" + p.getPrice() + " ");
        req.setAttribute("productID", p.getId());
        req.getRequestDispatcher("/WEB-INF/product page.jsp").forward(req, resp);
    }
}
