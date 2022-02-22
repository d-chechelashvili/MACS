import random
from dataclasses import dataclass
from typing import Callable, List, Optional, Protocol

from app.creature import (
    Creature,
    Evolvable,
    Fighter,
    ICreature,
    Movable,
    default_damage_calculation,
    default_move_strategy,
)

CLAW_NAMES = {0: "Small Claws", 1: "Medium Claws", 2: "Big claws"}


def creature_sort_func(x: ICreature) -> int:
    return int(x.position)


def random_init_creature() -> ICreature:
    creature = Creature()
    creature.position = random.randint(POSITION_RANGE_MIN, POSITION_RANGE_MAX)
    return creature


def evolve_random_features(evolvable: Evolvable) -> None:
    leg_count = random.randint(LEG_COUNT_MIN, LEG_COUNT_MAX)
    evolvable.evolve_legs(leg_count)

    wing_count = random.randint(WING_COUNT_MIN, WING_COUNT_MAX)
    evolvable.evolve_wings(wing_count)

    claw_type = random.randint(CLAW_TYPE_MIN, CLAW_TYPE_MAX)
    evolvable.evolve_claws(claw_type)

    teeth_sharpness = random.randint(TEETH_SHARPNESS_MIN, TEETH_SHARPNESS_MAX)
    evolvable.evolve_teeth(teeth_sharpness)


def default_fight_strategy(predator: Fighter, prey: Fighter) -> None:
    while True:
        predator.attack(prey, default_damage_calculation)
        if prey.health <= 0:
            return

        prey.attack(predator, default_damage_calculation)
        if predator.health <= 0:
            return


class SporeSimulator:
    game_on: bool = False

    def simulate(
        self, creature_count: int, creatures: Optional[List[ICreature]] = None
    ) -> None:
        self.game_on = True
        if creatures is None:
            creatures = self.__initialize_creatures(creature_count)
            self.__evolve_creatures(creatures)
        while self.game_on:
            self.__simulate_turn(creatures)
            self.__update_state(creatures)

    def __initialize_creatures(self, creature_count: int) -> List[ICreature]:
        creatures = []
        for _ in range(creature_count):
            creature = CreatureInitializer().init_creature()
            creatures.append(creature)
        return creatures

    def __evolve_creatures(self, creatures: List[ICreature]) -> None:
        CreatureEvolverLogger(CreatureEvolver()).evolve_creatures(creatures)

    # simulates one turn by sorting creatures by some logic then moving them and simulating fights if any happened
    def __simulate_turn(self, creatures: List[ICreature]) -> None:
        creatures.sort(key=creature_sort_func)
        self.__simulate_move(creatures)
        self.__simulate_fights(creatures)

    def __simulate_move(self, creatures: List[Movable]) -> None:
        for i, creature in enumerate(creatures):
            creature.move(move_strategy=default_move_strategy)

    # simulates fights if there was any intersection
    def __simulate_fights(self, creatures: List[ICreature]) -> None:
        for i in range(len(creatures) - 1):
            for j in range(i + 1, len(creatures)):
                if creatures[i].position >= creatures[j].position:
                    self.__simulate_fight(creatures[i], creatures[j])

    # simulates one fight
    def __simulate_fight(
        self,
        predator: Fighter,
        prey: Fighter,
        fight: Callable[[Fighter, Fighter], None] = default_fight_strategy,
    ) -> None:
        fight(predator, prey)  # FightSimulator for testing

    def __update_state(self, creatures: List[ICreature]) -> None:
        if creatures[-1].health <= 0:
            print("Some R rated things have happened")
            self.game_on = False
        else:
            for creature in creatures[:-1]:
                if creature.health <= 0 or creature.stamina <= 0:
                    creatures.remove(creature)
            if len(creatures) <= 1:
                print("Prey ran into infinity")
                self.game_on = False


# შემეძლო მქონოდა class Predator(Protocol), class Prey(Protocol),
# spore_simulator-ში მქონოდა predator: Predator, prey: Prey
# და interface segregation მექნა რამე განსხვავებული რომ ჰქონდეს რომელსაც დავმალავდით ან რამე
# მაგრამ ორივეს ზუსტად ერთიდაიგივე ქცევები აქვს და არ დამჭირდა.


class CreatureInitializer:
    def __init__(
        self, initialize_strategy: Callable[[], ICreature] = random_init_creature
    ) -> None:
        self.initialize_creature = initialize_strategy

    def init_creature(self) -> ICreature:
        return self.initialize_creature()


class ICreatureEvolver(Protocol):
    def evolve_creatures(self, creatures: List[ICreature]) -> None:
        pass


class CreatureEvolver:
    def __init__(
        self, evolve_strategy: Callable[[Evolvable], None] = evolve_random_features
    ) -> None:
        self.evolve_creature = evolve_strategy

    def evolve_creatures(self, creatures: List[ICreature]) -> None:
        for creature in creatures:
            self.evolve_creature(creature)


@dataclass
class CreatureEvolverDecorator:
    inner_creature_evolver: ICreatureEvolver

    def evolve_creatures(self, creatures: List[ICreature]) -> None:
        self.inner_creature_evolver.evolve_creatures(creatures)


# just for demonstration and training purposes
class CreatureEvolverLogger(CreatureEvolverDecorator):
    def evolve_creatures(self, creatures: List[ICreature]) -> None:
        super().evolve_creatures(creatures)
        for creature in creatures:
            print(
                "Pos: {}, legs: {}, wings: {}, claws: {}, teeth_sharpness: {}".format(
                    creature.position,
                    creature.leg_count,
                    creature.wing_count,
                    CLAW_NAMES[creature.claw_type],
                    creature.teeth_sharpness,
                )
            )


POSITION_RANGE_MIN = 0
POSITION_RANGE_MAX = 100

LEG_COUNT_MIN = 0
LEG_COUNT_MAX = 5

WING_COUNT_MIN = 0
WING_COUNT_MAX = 3

CLAW_TYPE_MIN = 0
CLAW_TYPE_MAX = 2

TEETH_SHARPNESS_MIN = 0
TEETH_SHARPNESS_MAX = 2
