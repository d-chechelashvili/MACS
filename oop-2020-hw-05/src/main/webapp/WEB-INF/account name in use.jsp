<%--
  Created by IntelliJ IDEA.
  User: Whiskeyjack
  Date: 04-Jul-20
  Time: 15:05
  To change this template use File | Settings | File Templates.
--%>
<%@ page contentType="text/html;charset=UTF-8" language="java" %>
<html>
<head>
    <title>Create Account</title>
</head>
<body>

<h2>The Name ${username} is Already In Use</h2>

<p>Please enter another name and password.</p>

<form action="signup" method="post">
    <label for="user-name-label">User name: </label>
    <input type="text" id="user-name-label" name="username"><br><br>
    <label for="password-label">Password: </label>
    <input type="text" id="password-label" name="password">
    <input type="submit" value="Login">
</form>

</body>
</html>
