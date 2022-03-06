from typing import List

import pytest

from app.spore_simulator import SporeSimulator
from app.tribe import EqualCharmDistributer, SingPerformanceChooser, Tribe
from app.tribe_member import (
    ITribeMember,
    TribeMember,
    TribeMemberWithFur,
    TribeMemberWithTails,
)


@pytest.fixture
def spore_simulator() -> SporeSimulator:
    return SporeSimulator()


class TestableTribe(Tribe):  # type: ignore
    def set_members(self, members: List[ITribeMember]) -> None:
        self._members = members


def test_should_charm(spore_simulator: SporeSimulator) -> None:
    guest = TestableTribe(
        possible_instruments=[],
        charm_distributer=EqualCharmDistributer(),
        performance_chooser=SingPerformanceChooser(),
    )
    guest_member = TribeMemberWithFur(
        tribe_member=TribeMemberWithTails(tribe_member=TribeMember(), modifier=0),
        modifier=1,
    )
    guest.set_members([guest_member])
    host = TestableTribe(
        possible_instruments=[],
        charm_distributer=EqualCharmDistributer(),
        performance_chooser=SingPerformanceChooser(),
    )
    host_member = TribeMember()
    host_member.composure = 100
    host.set_members([host_member])
    result = spore_simulator.simulate(guest=guest, host=host, evolution=False)
    assert result == 1


def test_should_not_charm(spore_simulator: SporeSimulator) -> None:
    host = TestableTribe(
        possible_instruments=[],
        charm_distributer=EqualCharmDistributer(),
        performance_chooser=SingPerformanceChooser(),
    )
    host_member = TribeMember()
    host_member.composure = 101
    host.set_members([host_member])
    guest = TestableTribe(
        possible_instruments=[],
        charm_distributer=EqualCharmDistributer(),
        performance_chooser=SingPerformanceChooser(),
    )
    guest_member = TribeMemberWithFur(
        tribe_member=TribeMemberWithTails(tribe_member=TribeMember(), modifier=0),
        modifier=1,
    )
    guest.set_members([guest_member])
    result = spore_simulator.simulate(guest=guest, host=host, evolution=False)
    assert result == 0
