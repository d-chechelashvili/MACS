from dataclasses import dataclass
from typing import Protocol


@dataclass
class Performance:
    energy: int = 0
    charm: int = 0

    def calculate_energy(self, old_energy: int) -> int:
        return old_energy - self.energy

    def calculate_charm(self) -> int:
        return self.charm


class Sing:
    def __call__(self, energy: int) -> Performance:
        if energy > 0:
            return Performance(energy=1, charm=1)
        return Performance()


class Dance:
    def __call__(self, energy: int) -> Performance:
        if energy > 20:
            return Performance(energy=2, charm=3)
        return Performance()


class Instrument(Protocol):
    def __call__(self, energy: int) -> Performance:
        pass


class Horn:
    def __call__(self, energy: int) -> Performance:
        if energy > 40:
            return Performance(energy=2, charm=4)
        return Performance()


class Maraca:
    def __call__(self, energy: int) -> Performance:
        if energy > 60:
            return Performance(energy=4, charm=6)
        return Performance()


class Didgeridoo:
    def __call__(self, energy: int) -> Performance:
        if energy > 80:
            return Performance(energy=4, charm=8)
        return Performance()
