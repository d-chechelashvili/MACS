from collections import defaultdict
from copy import copy
from typing import DefaultDict

import pytest

from app.creature import Creature
from app.spore_simulator import (
    CLAW_TYPE_MAX,
    CLAW_TYPE_MIN,
    LEG_COUNT_MAX,
    LEG_COUNT_MIN,
    POSITION_RANGE_MAX,
    POSITION_RANGE_MIN,
    TEETH_SHARPNESS_MAX,
    TEETH_SHARPNESS_MIN,
    WING_COUNT_MAX,
    WING_COUNT_MIN,
    CreatureEvolver,
    CreatureInitializer,
    SporeSimulator,
    default_fight_strategy,
)
from tests.test_creature import TestableCreature


@pytest.fixture
def spore_simulator() -> SporeSimulator:
    return SporeSimulator()


def test_initialize_creature() -> None:
    for _ in range(1000):
        creature = CreatureInitializer().init_creature()
        assert isinstance(creature, Creature)
        assert POSITION_RANGE_MIN <= creature.position <= POSITION_RANGE_MAX


def test_evolve_creatures() -> None:
    leg_values: DefaultDict[int, int] = defaultdict(int)
    wing_values: DefaultDict[int, int] = defaultdict(int)
    claw_values: DefaultDict[int, int] = defaultdict(int)
    teeth_values: DefaultDict[int, int] = defaultdict(int)
    creatures = [
        CreatureInitializer().init_creature(),
        CreatureInitializer().init_creature(),
    ]

    for _ in range(1000):
        CreatureEvolver().evolve_creatures(creatures)
        for creature in creatures:
            assert LEG_COUNT_MIN <= creature.leg_count <= LEG_COUNT_MAX
            leg_values[creature.leg_count] += 1
            assert WING_COUNT_MIN <= creature.wing_count <= WING_COUNT_MAX
            wing_values[creature.wing_count] += 1
            assert CLAW_TYPE_MIN <= creature.claw_type <= CLAW_TYPE_MAX
            claw_values[creature.claw_type] += 1
            assert (
                TEETH_SHARPNESS_MIN <= creature.teeth_sharpness <= TEETH_SHARPNESS_MAX
            )
            teeth_values[creature.teeth_sharpness] += 1

    assert len(leg_values.keys()) > 1
    assert len(wing_values.keys()) > 1
    assert len(claw_values.keys()) > 1
    assert len(teeth_values.keys()) > 1


def test_simulate_fight() -> None:
    creature1 = TestableCreature(claw_type=2, teeth_sharpness=1)
    creature2 = TestableCreature(claw_type=1, teeth_sharpness=2)
    default_fight_strategy(creature1, creature2)
    assert creature1.health == 10
    assert creature2.health == -12
    creature1 = TestableCreature(claw_type=2, teeth_sharpness=1)
    creature2 = TestableCreature(claw_type=2, teeth_sharpness=2)
    default_fight_strategy(creature1, creature2)
    assert creature1.health == -20
    assert creature2.health == 16
    creature1 = TestableCreature(claw_type=0, teeth_sharpness=0)
    creature2 = TestableCreature(claw_type=0, teeth_sharpness=0)
    default_fight_strategy(creature1, creature2)
    assert creature1.health == 4
    assert creature2.health == -4
    creature1 = TestableCreature(claw_type=0, teeth_sharpness=1)
    creature2 = TestableCreature(claw_type=1, teeth_sharpness=0)
    default_fight_strategy(creature1, creature2)
    assert creature1.health == 16
    assert creature2.health == -12


def test_simulate_catch(spore_simulator: SporeSimulator) -> None:
    creatures = [
        TestableCreature(
            position=1, claw_type=2, teeth_sharpness=1, leg_count=2, wing_count=3
        ),
        TestableCreature(
            position=13, claw_type=1, teeth_sharpness=2, leg_count=0, wing_count=1
        ),
    ]
    spore_simulator.simulate(len(creatures), creatures=copy(creatures))
    assert creatures[0].health == 10
    assert creatures[1].health == -12
    assert creatures[0].stamina == 92
    assert creatures[1].stamina == 98
    assert creatures[0].position == 17
    assert creatures[1].position == 15


def test_simulate_escape(spore_simulator: SporeSimulator) -> None:
    creatures = [
        TestableCreature(
            position=1, claw_type=2, teeth_sharpness=1, leg_count=2, wing_count=2
        ),
        TestableCreature(
            position=61, claw_type=1, teeth_sharpness=2, leg_count=1, wing_count=1
        ),
    ]
    spore_simulator.simulate(len(creatures), creatures=copy(creatures))
    assert creatures[0].health == 100
    assert creatures[1].health == 100
    assert creatures[0].stamina == 0
    assert creatures[1].stamina == 10
    assert creatures[0].position == 161
    assert creatures[1].position == 191
