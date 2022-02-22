<%@ page import="Model.Product" %>
<%@ page import="Model.StudentStoreDatabase" %>
<%@ page import="java.sql.SQLException" %>
<%@ page import="java.util.ArrayList" %><%--
  Created by IntelliJ IDEA.
  User: Whiskeyjack
  Date: 04-Jul-20
  Time: 22:41
  To change this template use File | Settings | File Templates.
--%>
<%@ page contentType="text/html;charset=UTF-8" language="java" %>
<html>
<head>
    <title>Student Store</title>
</head>
<body>
<h1>Student Store</h1>

<p>Items available:</p>
<ul>
    <% try {
        ArrayList<Product> products = ((StudentStoreDatabase) request.getServletContext().
                getAttribute(StudentStoreDatabase.ATTRIBUTE)).getProductList();
        for (Product p : products) {
    %>
    <li><a href="showProduct?productId=<%= p.getId()%>"><%= p.getName() %>
    </a></li>
    <%
            }
        } catch (SQLException throwables) {
            throwables.printStackTrace();
        }
    %>


</ul>
</body>
</html>
