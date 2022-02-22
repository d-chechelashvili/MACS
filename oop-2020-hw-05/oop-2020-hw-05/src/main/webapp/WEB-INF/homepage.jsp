<%--
  Created by IntelliJ IDEA.
  User: Whiskeyjack
  Date: 04-Jul-20
  Time: 14:06
  To change this template use File | Settings | File Templates.
--%>
<%@ page contentType="text/html;charset=UTF-8" language="java" %>
<html>
<head>
    <title>Welcome</title>
</head>
<body>

<h2>Welcome to Homework 5</h2>

<p>Please log in.</p>


<form action="" method="post">
    <label for="user-name-label">User name: </label>
    <input type="text" id="user-name-label" name="username"><br><br>
    <label for="password-label">Password: </label>
    <input type="text" id="password-label" name="password">
    <input type="submit" value="Login"><br><br>
</form>

<a href="signup">Create New Account</a>

</body>
</html>
