# Bonus Project

## Intro

The assignment is inspired by the 2008 life simulation RTS game [Spore](https://www.spore.com/). Spore allows you to evolve creatures with different characteristics such as feathers, fur, tails, and many many more.

## Simulation

On the Tribal Stage of Spore, your tribe can socialize with others in order to charm and befriend them. Your task is to model a simple simulation of such socialization between the guest and host tribes. Provide 100 simulations of the interaction between tribes as follows:

- Evolve a random guest tribe (log characteristics of each member)
- Evolve a random host tribe (log characteristics of each member)
- Give the random number of each instrument (Horn, Maraca, Didgeridoo) to each tribe
- Guest tribe performs for hosts until:
  * Either all hosts run out of composure. (log message: "Guests have charmed hosts.")
  * Or all guests run out of energy. (log message: "Hosts have disappointed guests.")

## Technical Details

- Each tribal member can evolve `tails` (+3 charm per tail)
- Each tribal member can evolve `fur`
  * `Regular` x2 charm
  * `Stripped` x4 charm
  * `Dotted` x6 charm
- Each tribe can have up to 5 members
- Each tribe can own up to 15 instruments
- Naturally, to perform with an instrument the tribe must own the said instrument first
- Each guest tribe member can decide on which performance to choose at any particular moment
  * (Assuming they have enough energy and instruments)
- Guests emit charm during the performance that shakes the composure of the hosts
- You may distribute total emitted charm across the host tribe however you like. A couple of ideas:
  * Equally among all host tribe members
  * Pick a random member and apply all charms to them
- Each performance has its requirements and effects shown in the following table:

| Performance   | Requires Energy   | Uses Energy   | Charm   |
|-------------  |-----------------  |-------------  |-------  |
| Sing          | 0+                | 1             | 1       |
| Dance         | 20+               | 2             | 3       |
| Horn          | 40+               | 2             | 4       |
| Maraca        | 60+               | 4             | 6       |
| Didgeridoo    | 80+               | 4             | 8       |

## Expected Patterns

You should not need more than, Decorator, Strategy, and Composite. Please do not go overboard with multiple Chains of Responsibilities or Singletons.

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

- 80%: Design
- 20%: Testing

Note: we will not grade your assignment if your code fails any of the Linting/formatting rules and you will automatically get zero.

## Disclaimer

We reserve the right to penalize you for violating well-known software principles that you covered in previous courses such as decomposition or DRY. We sincerely ask you to not make a mess of your code.
