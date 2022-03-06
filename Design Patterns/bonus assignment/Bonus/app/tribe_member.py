from dataclasses import dataclass, field
from typing import Callable, Protocol

from app.performance import Performance, Sing


class ITribeMember(Protocol):
    def perform(self) -> int:
        pass

    def has_composure(self) -> bool:
        pass

    def has_energy(self) -> bool:
        pass

    def distribute_charm(self, charm: int) -> None:
        pass

    def get_energy(self) -> int:
        pass

    def get_composure(self) -> int:
        pass

    def set_next_performance(
        self, next_performance: Callable[[int], Performance]
    ) -> None:
        pass

    def __str__(self) -> str:
        pass


@dataclass
class TribeMember:
    energy: int = 100
    composure: int = 3700

    next_performance: Callable[[int], Performance] = field(default_factory=Sing)

    def perform(self) -> int:
        performance = self.next_performance(self.energy)
        charm: int = performance.calculate_charm()
        self.energy = performance.calculate_energy(self.energy)
        return charm

    def has_composure(self) -> bool:
        return self.composure > 0

    def has_energy(self) -> bool:
        return self.energy > 0

    def distribute_charm(self, charm: int) -> None:
        self.composure -= charm

    def get_energy(self) -> int:
        return self.energy

    def get_composure(self) -> int:
        return self.composure

    def set_next_performance(
        self, next_performance: Callable[[int], Performance]
    ) -> None:
        self.next_performance = next_performance

    def __str__(self) -> str:
        return f"Energy: {self.get_energy()}, Composure: {self.get_composure()}"


@dataclass
class TribeMemberDecorator:
    tribe_member: ITribeMember

    def perform(self) -> int:
        return self.tribe_member.perform()

    def has_composure(self) -> bool:
        return self.tribe_member.has_composure()

    def has_energy(self) -> bool:
        return self.tribe_member.has_energy()

    def distribute_charm(self, charm: int) -> None:
        self.tribe_member.distribute_charm(charm)

    def get_energy(self) -> int:
        return self.tribe_member.get_energy()

    def get_composure(self) -> int:
        return self.tribe_member.get_composure()

    def set_next_performance(
        self, next_performance: Callable[[int], Performance]
    ) -> None:
        self.tribe_member.set_next_performance(next_performance)

    def __str__(self) -> str:
        return self.tribe_member.__str__()


@dataclass
class TribeMemberWithTails(TribeMemberDecorator):
    modifier: int

    def perform(self) -> int:
        return super().perform() + self.modifier

    def __str__(self) -> str:
        return self.tribe_member.__str__() + f", Tail_modifier: {self.modifier}"


@dataclass
class TribeMemberWithFur(TribeMemberDecorator):
    modifier: int

    def perform(self) -> int:
        return super().perform() * self.modifier

    def __str__(self) -> str:
        return self.tribe_member.__str__() + f", Fur_modifier: {self.modifier}"
