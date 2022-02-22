# Assignment 1

## Intro

The assignment is inspired by the 2008 life simulation RTS game [Spore](https://www.spore.com/). Spore allows you to evolve creatures with different characteristics such as claws, sharp teeth, wings, and many many more.

## Simulation

Your task is to model a simple simulation of the interaction between two creatures Predator and Pray. For now, let's imaging that the shape of the world we are modeling is an infinite ray illustrated below. It starts at position 0 and goes on and on as shown below.

```
012345678...
------------
```

Provide the simulation of the interaction between two creatures as follows

- Evolve a random predator at a random location (log characteristics)
- Evolve a random pray at a random location (log characteristics)
- Predator chases pray, pray runs away until:
  * Predator runs out of stamina. (log message: "Pray ran into infinity")
  * Predator catches pray. In this case, they enter the fight
- In the fight, both creatures attack until:
  * Predator runs out of health. (log message: "Pray ran into infinity")
  * Pray runs out of health. (log message: "Some R rated things have happened") :D

Repeat 100 times

Note: During the chase phase it is up to you to decide how smart or desperate you want your creatures to be.

## Functional requirements

- Ability to evolve `Legs`. Creature needs
  * At least one leg to `Hop`
  * At least two legs to `Walk` or `Run`
- Ability to evolve `Wings`. Creature needs
  * At least two wings to `Fly`
- Ability to evolve `Claws`
  * Small claws multiply creatures attacking power x2
  * Medium claws multiply creatures attacking power x3
  * Big claws multiply creatures attacking power x4
- Ability to evolve mouth with `Teeth`
  * Different degree of sharpness with +3, +6, +9 boost to attacking power
- Ability to `move` the creature in the world
  * Creatures with no limb or wings can Crawl
- Ability to `attack` other creatures

Movement | requires stamina | uses stamina | speed |
---------|------------------|--------------|-------|
Crawling | 0+               | 1            | 1     |
Hopping  | 20+              | 2            | 3     |
Walking  | 40+              | 2            | 4     |
Running  | 60+              | 4            | 6     |
Flying   | 80+              | 4            | 8     |


## Technical Details

Initial `Creature` class is provided.

```python
from dataclasses import dataclass

@dataclass
class Creature:
    position: int
    power: int = 1
    health: int = 100
    stamina: int = 100

```

- *Position:* Non-negative integer. The location of the creature in the world
- *Power:* Non-negative integer. Each creature can deal base damage equal to their power
- *Stamina:* Non-negative integer. Movement consumes stamina
- *health:* Non-negative integer. Each attack experienced by the creature takes away the health

You need to design the "API" of this class and any other artifacts necessary for the simulation.

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

## Disclamer

We reserve the right to penalize you for violating well-known software principles that you covered in previous courses such as decomposition or DRY. We sincerely ask you to not make a mess of your code.
