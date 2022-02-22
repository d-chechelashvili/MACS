<%@ page import="Model.Product" %>
<%@ page import="Model.ShoppingCart" %>
<%@ page import="Model.StudentStoreDatabase" %>
<%@ page import="java.util.Set" %>
<%--
  Created by IntelliJ IDEA.
  User: Whiskeyjack
  Date: 05-Jul-20
  Time: 16:46
  To change this template use File | Settings | File Templates.
--%>
<%@ page contentType="text/html;charset=UTF-8" language="java" %>
<html>
<head>
    <title>Shopping Cart</title>
</head>
<body>

<h1>Shopping Cart</h1>
<form action="shopping-cart" method="post">
    <ul>
        <%
            ShoppingCart cart = (ShoppingCart) request.getSession()
                    .getAttribute(ShoppingCart.ATTRIBUTE);
            Set<String> items = cart.getItems();
            StudentStoreDatabase db = (StudentStoreDatabase) request.getServletContext().
                    getAttribute(StudentStoreDatabase.ATTRIBUTE);

            double totalCost = 0;
            for (String productID : items) {
                Product p = db.findProduct(productID);
        %>
        <li><label>
            <input name="<%=p.getId()%>" type="number" value="<%=cart.getCount(productID)%>">
            <%=p.getName() + ", " + p.getPrice()%>
        </label></li>
        <%
                totalCost += cart.getCount(productID) * p.getPrice();
            }

        %>
    </ul>
    <label>Total: $<%=totalCost%> <input type="submit" value="Update Cart"></label>

</form>
<a href="/" style="text-align:center">Continue Shopping</a>
</body>
</html>
