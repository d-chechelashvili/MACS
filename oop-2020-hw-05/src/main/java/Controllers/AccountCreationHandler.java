package Controllers;

import model.AccountManager;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;

public class AccountCreationHandler extends HttpServlet {
    @Override
    protected void doGet(HttpServletRequest req, HttpServletResponse resp) throws ServletException, IOException {
        req.getRequestDispatcher("/WEB-INF/create account.jsp").forward(req, resp);
    }

    @Override
    protected void doPost(HttpServletRequest req, HttpServletResponse resp) throws ServletException, IOException {
        AccountManager am = (AccountManager) getServletContext().getAttribute(AccountManager.ATTRIBUTE);

        String userName = req.getParameter("username");
        String password = req.getParameter("password");

        req.setAttribute("username", userName);
        if(!am.accountExists(userName)){
            am.createAccount(userName, password);
            req.getRequestDispatcher("/WEB-INF/user welcome.jsp").forward(req, resp);
        }else{
            req.getRequestDispatcher("WEB-INF/account name in use.jsp").forward(req, resp);
        }
    }
}
