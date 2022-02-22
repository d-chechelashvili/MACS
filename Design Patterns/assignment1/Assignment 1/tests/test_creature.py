import random
from typing import Optional

from app.creature import (
    CRAWLING_POSITION_CHANGE,
    CRAWLING_REQUIRED_LEGS,
    CRAWLING_STAMINA_CHANGE,
    FLYING_POSITION_CHANGE,
    FLYING_REQUIRED_STAMINA,
    FLYING_REQUIRED_WINGS,
    FLYING_STAMINA_CHANGE,
    HOPPING_POSITION_CHANGE,
    HOPPING_REQUIRED_LEGS,
    HOPPING_STAMINA_CHANGE,
    RUNNING_POSITION_CHANGE,
    RUNNING_REQUIRED_LEGS,
    RUNNING_REQUIRED_STAMINA,
    RUNNING_STAMINA_CHANGE,
    WALKING_POSITION_CHANGE,
    WALKING_REQUIRED_LEGS,
    WALKING_REQUIRED_STAMINA,
    WALKING_STAMINA_CHANGE,
    Creature,
)
from app.spore_simulator import (
    CLAW_TYPE_MAX,
    CLAW_TYPE_MIN,
    LEG_COUNT_MAX,
    LEG_COUNT_MIN,
    TEETH_SHARPNESS_MAX,
    TEETH_SHARPNESS_MIN,
    WING_COUNT_MAX,
    WING_COUNT_MIN,
)


def test_creature_creation_health() -> None:
    creature = Creature()
    assert creature.health == 100


def test_creature_creation_stamina() -> None:
    creature = Creature()
    assert creature.stamina == 100


def test_creature_creation_power() -> None:
    creature = Creature()
    assert creature.power == 1


def test_creature_evolution_legs() -> None:
    creature = Creature()

    for i in range(100):
        leg_count = random.randint(LEG_COUNT_MIN, LEG_COUNT_MAX)
        creature.evolve_legs(leg_count)
        assert creature.leg_count == leg_count


def test_creature_evolution_wings() -> None:
    creature = Creature()
    for i in range(100):
        wing_count = random.randint(WING_COUNT_MIN, WING_COUNT_MAX)
        creature.evolve_wings(wing_count)
        assert creature.wing_count == wing_count


def test_creature_evolution_claws() -> None:
    creature = Creature()
    for i in range(100):
        claw_type = random.randint(CLAW_TYPE_MIN, CLAW_TYPE_MAX)
        creature.evolve_claws(claw_type)
        assert creature.claw_type == claw_type


def test_creature_evolution_teeth() -> None:
    creature = Creature()
    for i in range(100):
        teeth_sharpness = random.randint(TEETH_SHARPNESS_MIN, TEETH_SHARPNESS_MAX)
        creature.evolve_teeth(teeth_sharpness)
        assert creature.teeth_sharpness == teeth_sharpness


class TestableCreature(Creature):  # type: ignore
    def __init__(
        self,
        position: Optional[int] = None,
        power: Optional[int] = None,
        health: Optional[int] = None,
        stamina: Optional[int] = None,
        wing_count: Optional[int] = None,
        leg_count: Optional[int] = None,
        claw_type: Optional[int] = None,
        teeth_sharpness: Optional[int] = None,
    ) -> None:
        if position is not None:
            self.position = position
        if power is not None:
            self.power = power
        if health is not None:
            self.health = health
        if stamina is not None:
            self.stamina = stamina
        if wing_count is not None:
            self.wing_count = wing_count
        if leg_count is not None:
            self.leg_count = leg_count
        if claw_type is not None:
            self.claw_type = claw_type
        if teeth_sharpness is not None:
            self.teeth_sharpness = teeth_sharpness


def test_move_no_stamina() -> None:
    for i in range(LEG_COUNT_MIN, LEG_COUNT_MAX + 1):
        creature = TestableCreature(leg_count=i, stamina=0)
        pos = creature.position
        creature.move()
        assert creature.position == pos


def move_test(
    num: int, creature: Creature, position_change: int, stamina_change: int
) -> None:
    for _ in range(num):
        pos = creature.position
        stamina = creature.stamina
        creature.move()
        assert creature.position == pos + position_change
        assert creature.stamina == stamina + stamina_change


def test_move_crawl() -> None:
    creature = TestableCreature(leg_count=CRAWLING_REQUIRED_LEGS, stamina=100)
    move_test(80, creature, CRAWLING_POSITION_CHANGE, CRAWLING_STAMINA_CHANGE)
    creature.evolve_legs(random.randint(CRAWLING_REQUIRED_LEGS + 1, LEG_COUNT_MAX))
    move_test(20, creature, CRAWLING_POSITION_CHANGE, CRAWLING_STAMINA_CHANGE)


def test_move_hop() -> None:
    creature = TestableCreature(leg_count=1, stamina=80)
    move_test(20, creature, HOPPING_POSITION_CHANGE, HOPPING_STAMINA_CHANGE)
    creature.evolve_legs(random.randint(HOPPING_REQUIRED_LEGS + 1, LEG_COUNT_MAX))
    move_test(10, creature, HOPPING_POSITION_CHANGE, HOPPING_STAMINA_CHANGE)


def test_move_fly() -> None:
    creature = TestableCreature(
        stamina=FLYING_REQUIRED_STAMINA,
        wing_count=random.randint(FLYING_REQUIRED_WINGS, WING_COUNT_MAX),
    )
    move_test(1, creature, FLYING_POSITION_CHANGE, FLYING_STAMINA_CHANGE)


def test_move_run() -> None:
    creature = TestableCreature(
        stamina=RUNNING_REQUIRED_STAMINA,
        leg_count=random.randint(RUNNING_REQUIRED_LEGS, LEG_COUNT_MAX),
    )
    move_test(1, creature, RUNNING_POSITION_CHANGE, RUNNING_STAMINA_CHANGE)


def test_move_walk() -> None:
    creature = TestableCreature(
        stamina=WALKING_REQUIRED_STAMINA,
        leg_count=random.randint(WALKING_REQUIRED_LEGS, LEG_COUNT_MAX),
    )
    move_test(1, creature, WALKING_POSITION_CHANGE, WALKING_STAMINA_CHANGE)


def test_fight() -> None:
    creature1 = TestableCreature(claw_type=0, teeth_sharpness=2)
    creature2 = TestableCreature(health=50)
    creature1.attack(creature2)
    assert creature2.health == 30
