import random
from random import randint
from typing import List, Protocol

from app.performance import Dance, Instrument, Sing
from app.tribe_member import (
    ITribeMember,
    TribeMember,
    TribeMemberWithFur,
    TribeMemberWithTails,
)


class Evolvable(Protocol):
    def evolve(self) -> None:
        pass


class Performer(Protocol):
    def perform(self) -> int:
        pass


class Charmable(Protocol):
    def distribute_charm(self, charm: int) -> None:
        pass


class Composed(Protocol):
    def has_composure(self) -> bool:
        pass


class Energetic(Protocol):
    def has_energy(self) -> bool:
        pass


class IHost(Protocol):
    def evolve(self) -> None:
        pass

    def distribute_charm(self, charm: int) -> None:
        pass

    def has_composure(self) -> bool:
        pass


class IGuest(Protocol):
    def evolve(self) -> None:
        pass

    def perform(self) -> int:
        pass

    def has_energy(self) -> bool:
        pass


class ITribe(Protocol):
    def evolve(self) -> None:
        pass

    def perform(self) -> int:
        pass

    def distribute_charm(self, charm: int) -> None:
        pass

    def has_energy(self) -> bool:
        pass

    def has_composure(self) -> bool:
        pass


class Tribe:
    _members: List[ITribeMember]
    _instruments: List[Instrument]
    performance_chooser: "IChoosePerformance"
    charm_distributer: "IDistributeCharm"

    def __init__(
        self,
        possible_instruments: List[Instrument],
        performance_chooser: "IChoosePerformance",
        charm_distributer: "IDistributeCharm",
    ) -> None:
        member_count = randint(1, 5)
        self._members = [TribeMember() for _ in range(member_count)]
        instrument_count = randint(1, 15)
        if len(possible_instruments) > 0:
            self._instruments = [
                random.choice(possible_instruments) for _ in range(instrument_count)
            ]
        else:
            self._instruments = []
        self.performance_chooser = performance_chooser
        self.charm_distributer = charm_distributer

    def evolve(self) -> None:
        for i in range(len(self._members)):
            tail_count = randint(0, 2)
            self._members[i] = TribeMemberWithTails(
                tribe_member=self._members[i], modifier=tail_count * 3
            )
            fur_type = random.choice([1, 2, 4, 6])
            self._members[i] = TribeMemberWithFur(
                tribe_member=self._members[i], modifier=fur_type
            )
            print(f"Tribe Member: {self._members[i]}")

    def perform(self) -> int:
        total_charm = 0
        used_instruments = []
        random.shuffle(self._members)
        for member in self._members:
            index = self.performance_chooser.choose_performance(
                self._instruments, member
            )
            total_charm += member.perform()
            if index >= 0:
                used_instrument = self._instruments.pop(index)
                used_instruments.append(used_instrument)
        self._instruments.extend(used_instruments)
        return total_charm

    def distribute_charm(self, charm: int) -> None:
        self.charm_distributer.distribute_charm(charm, self._members)

    def has_energy(self) -> bool:
        for member in self._members:
            if member.has_energy():
                return True
        return False

    def has_composure(self) -> bool:
        for member in self._members:
            if member.has_composure():
                return True
        return False


class IChoosePerformance(Protocol):
    def choose_performance(
        self, available_instruments: List[Instrument], member: ITribeMember
    ) -> int:  # returns index of chosen instrument, negative means that no instrument was chosen
        pass


class SingPerformanceChooser:
    def choose_performance(
        self, available_instruments: List[Instrument], member: ITribeMember
    ) -> int:
        member.set_next_performance(next_performance=Sing())
        return -1


class RandomPerformanceChooser:
    def choose_performance(
        self, available_instruments: List[Instrument], member: ITribeMember
    ) -> int:
        choice_list = [Sing(), Dance()]
        choice_list.extend(available_instruments)
        index = random.randint(0, len(choice_list) - 1)
        member.set_next_performance(next_performance=choice_list[index])
        return index - 2


class IDistributeCharm(Protocol):
    def distribute_charm(self, charm: int, tribe_members: List[ITribeMember]) -> None:
        pass


class EqualCharmDistributer:
    def distribute_charm(self, charm: int, tribe_members: List[ITribeMember]) -> None:
        equal_charm = charm // len(tribe_members)
        remainder = charm % len(tribe_members)
        for member in tribe_members:
            member_charm = equal_charm
            if remainder > 0:
                member_charm += 1
                remainder -= 1
            member.distribute_charm(member_charm)
