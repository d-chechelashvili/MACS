<%--
  Created by IntelliJ IDEA.
  User: Whiskeyjack
  Date: 05-Jul-20
  Time: 13:05
  To change this template use File | Settings | File Templates.
--%>
<%@ page contentType="text/html;charset=UTF-8" language="java" %>
<html>
<head>
    <title>${productName}</title>
</head>
<body>

<h1>${productName}</h1>

<img src="store-images/${imageFile}" alt="${productName}">

<form action="shopping-cart" method="post">
    <input name="productID" type="hidden" value="${productID}"/>
    <label>${price}<input type="submit" value="Add to Cart"></label>
</form>

</body>
</html>
