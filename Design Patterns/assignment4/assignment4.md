# Assignment 4

## Intro

Within the scope of this assignment, we are going to design a service with HTTP API for Point of Sales (POS) system.

## User Stories

> As a *cashier*, I would like to open a receipt so that, I can serve a customer.

> As a *cashier*, I would like to add items to an open receipt so that, I can calculate how much the customer needs to pay.

> As a *customer*, I would like to see a receipt with all my items so that, I know how much I have to pay.

> As a *cashier*, I would like to close the paid receipt so that, I can start serving the next customer.

> As a *store manager*, I would like to make X reports so that, I can see the state of the store.


## Technical Details

- Store can sell items as singles
- Store can sell items as batches/packs. (think 6-pack of beer cans :D)
- Closed receipt cannot be modified.
- Receipt contains list of items with units sold, price and total price, as well as a grand total price of all items.
- X reports contain revenue, count of each item sold, and number of closed receipts (on a given day).
- Registration/Authorization is out of scope, all endpoints should be available for everyone.
- UI and Concurancy are out of scope.
- Use [SQLite](https://docs.python.org/3/library/sqlite3.html)) for persistence.
- Use [FastAPI](https://fastapi.tiangolo.com/) as a web framework.


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

- 20%: Architecture
- 25%: Usage of design patterns
- 25%: Usage of S.O.L.I.D. principals
- 20%: Unit testing
- 10%: Linting/formatting

## Disclaimer

We reserve the right to penalize you for violating well-known software principles that you covered in previous courses such as decomposition or DRY. We sincerely ask you to not make a mess of your code.
