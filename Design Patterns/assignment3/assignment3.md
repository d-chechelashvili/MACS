# Assignment 3

Implement simple REPL (Read Evaluate Print Loop) expecting the following commands

 - `subscribe <username> to <channel>`: from now on given user should get a notification when the video is published on the given channel
 - `publish video on <channel>`: print usernames of the subscribers that will get notified

It is expected that you will show the understanding of Observer and Repository patterns as well as your ability to test on the "Edge"

## Persistence

You can use [SQLite](https://docs.python.org/3/library/sqlite3.html)) but a simple text file with your own made-up format will also suffice.

## Example

```
>>> subscribe <Alice> to <Continuous Delivery>
Alice subscribed to Continuous Delivery
>>> subscribe <Bob> to <Continuous Delivery>
Bob subscribed to Continuous Delivery
>>> publish video on <Continuous Delivery>
Notifying subscribers of Continuous Delivery:
  Alice
  Bob
```

## Unit testing

Provide unit tests that prove the correctness of your software artifacts

## Linting/formatting

Format your code using `black` auto formatter

Sort your imports with `isort` using the following configuration:

```
[settings]
profile = black
```

Check your static types with `mypy` using the following configuration:

```
[mypy]
python_version = 3.9
ignore_missing_imports = True
strict = True
```

Check your code with `flake8` using the following configuration:

```
[flake8]
max-line-length = 88
select = C,E,F,W,B,B950
ignore = E501,W503
```

## Grading

- 40%: Usage of design patterns
- 40%: Usage of S.O.L.I.D. principals
- 10%: Unit testing
- 10%: Linting/formatting

## Disclaimer

We reserve the right to penalize you for violating well-known software principles that you covered in previous courses such as decomposition or DRY. We sincerely ask you to not make a mess of your code.