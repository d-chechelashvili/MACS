package model;

import javax.servlet.ServletContext;
import javax.servlet.ServletContextEvent;
import javax.servlet.ServletContextListener;
import javax.servlet.annotation.WebListener;

@WebListener
public class WebAppInitializer implements ServletContextListener {

    @Override
    public void contextInitialized(ServletContextEvent sce) {
        AccountManager am = new AccountManager();
        ServletContext sc  = sce.getServletContext();
        sc.setAttribute(AccountManager.ATTRIBUTE, am);
    }
}
