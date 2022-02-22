from abc import ABC, abstractmethod
from typing import Callable


def invoke_movement(
    creature: "Movable", position_change: int, stamina_change: int
) -> int:
    creature.stamina += stamina_change
    creature.position += position_change
    return position_change


def default_move_strategy(creature: "Movable") -> None:
    if (
        creature.wing_count >= FLYING_REQUIRED_WINGS
        and creature.stamina >= FLYING_REQUIRED_STAMINA
    ):
        invoke_movement(creature, FLYING_POSITION_CHANGE, FLYING_STAMINA_CHANGE)
    elif (
        creature.leg_count >= RUNNING_REQUIRED_LEGS
        and creature.stamina >= RUNNING_REQUIRED_STAMINA
    ):
        invoke_movement(creature, RUNNING_POSITION_CHANGE, RUNNING_STAMINA_CHANGE)
    elif (
        creature.leg_count >= WALKING_REQUIRED_LEGS
        and creature.stamina >= WALKING_REQUIRED_STAMINA
    ):
        invoke_movement(creature, WALKING_POSITION_CHANGE, WALKING_STAMINA_CHANGE)
    elif (
        creature.leg_count >= HOPPING_REQUIRED_LEGS
        and creature.stamina >= HOPPING_REQUIRED_STAMINA
    ):
        invoke_movement(creature, HOPPING_POSITION_CHANGE, HOPPING_STAMINA_CHANGE)
    elif (
        creature.leg_count >= CRAWLING_REQUIRED_LEGS
        and creature.stamina >= CRAWLING_REQUIRED_STAMINA
    ):
        invoke_movement(creature, CRAWLING_POSITION_CHANGE, CRAWLING_STAMINA_CHANGE)


def default_damage_calculation(creature: "Fighter") -> int:
    return int(
        (creature.power + (creature.teeth_sharpness + 1) * 3) * (creature.claw_type + 2)
    )


class Movable(ABC):
    position: int
    stamina: int
    wing_count: int
    leg_count: int

    @abstractmethod
    def move(self, move_strategy: Callable[["Movable"], None]) -> None:
        pass


class Evolvable(ABC):
    wing_count: int
    leg_count: int
    claw_type: int
    teeth_sharpness: int

    @abstractmethod
    def evolve_legs(self, count: int) -> None:
        pass

    @abstractmethod
    def evolve_claws(self, type: int) -> None:
        pass

    @abstractmethod
    def evolve_wings(self, count: int) -> None:
        pass

    @abstractmethod
    def evolve_teeth(self, sharpness: int) -> None:
        pass


class Fighter(ABC):
    power: int
    health: int
    claw_type: int
    teeth_sharpness: int

    @abstractmethod
    def attack(
        self, enemy: "Fighter", damage_calculation_strategy: Callable[["Fighter"], int]
    ) -> None:
        pass

    @abstractmethod
    def take_damage(self, damage: int) -> None:
        pass


class ICreature(Evolvable, Movable, Fighter, ABC):
    position: int = 0


class Creature(ICreature):
    power: int = 1
    health: int = 100
    stamina: int = 100
    wing_count: int = 0
    leg_count: int = 0
    claw_type: int = 0
    teeth_sharpness: int = 0

    def evolve_legs(self, count: int) -> None:
        self.leg_count = count

    def evolve_claws(self, type: int) -> None:
        self.claw_type = type

    def evolve_wings(self, count: int) -> None:
        self.wing_count = count

    def evolve_teeth(self, sharpness: int) -> None:
        self.teeth_sharpness = sharpness

    def move(
        self, move_strategy: Callable[["Movable"], None] = default_move_strategy
    ) -> None:
        move_strategy(self)

    def attack(
        self,
        enemy: Fighter,
        damage_calculation_strategy: Callable[
            [Fighter], int
        ] = default_damage_calculation,
    ) -> None:
        damage = damage_calculation_strategy(self)
        enemy.take_damage(damage)

    def take_damage(self, damage: int) -> None:
        self.health -= damage


FLYING_REQUIRED_WINGS = 2
RUNNING_REQUIRED_LEGS = 2
WALKING_REQUIRED_LEGS = 2
HOPPING_REQUIRED_LEGS = 1
CRAWLING_REQUIRED_LEGS = 0

FLYING_POSITION_CHANGE = 8
RUNNING_POSITION_CHANGE = 6
WALKING_POSITION_CHANGE = 4
HOPPING_POSITION_CHANGE = 3
CRAWLING_POSITION_CHANGE = 1

FLYING_REQUIRED_STAMINA = 81
RUNNING_REQUIRED_STAMINA = 61
WALKING_REQUIRED_STAMINA = 41
HOPPING_REQUIRED_STAMINA = 21
CRAWLING_REQUIRED_STAMINA = 1

FLYING_STAMINA_CHANGE = -4
RUNNING_STAMINA_CHANGE = -4
WALKING_STAMINA_CHANGE = -2
HOPPING_STAMINA_CHANGE = -2
CRAWLING_STAMINA_CHANGE = -1
