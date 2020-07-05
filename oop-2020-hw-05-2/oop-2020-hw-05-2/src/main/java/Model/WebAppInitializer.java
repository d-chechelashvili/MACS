package Model;

import javax.servlet.ServletContext;
import javax.servlet.ServletContextEvent;
import javax.servlet.ServletContextListener;
import javax.servlet.annotation.WebListener;
import javax.servlet.http.HttpSession;
import javax.servlet.http.HttpSessionEvent;
import javax.servlet.http.HttpSessionListener;

@WebListener
public class WebAppInitializer implements ServletContextListener, HttpSessionListener {

    @Override
    public void contextInitialized(ServletContextEvent sce) {
        StudentStoreDatabase db = new StudentStoreDatabase();
        ServletContext sc = sce.getServletContext();
        sc.setAttribute(StudentStoreDatabase.ATTRIBUTE, db);
    }

    @Override
    public void sessionCreated(HttpSessionEvent se) {
        ShoppingCart cart = new ShoppingCart();
        HttpSession session = se.getSession();
        session.setAttribute(ShoppingCart.ATTRIBUTE, cart);
    }
}
