from typing import Protocol


class ICalculatePrice(Protocol):
    def calculate_price(self) -> float:
        pass


class IItem(Protocol):
    def calculate_price(self) -> float:
        pass

    def get_name(self) -> str:
        pass

    def get_count(self) -> int:
        pass

    def get_single_price(self) -> float:
        pass

    def __eq__(self, other: object) -> bool:
        pass

    def __hash__(self) -> int:
        pass


class Item:
    _name: str
    _count: int
    _single_price: float

    def __init__(self, name: str, single_price: float, count: int = 1):
        self._name = name
        self._count = count
        self._single_price = single_price

    def calculate_price(self) -> float:
        return self._single_price * self._count

    def get_name(self) -> str:
        return self._name

    def get_count(self) -> int:
        return self._count

    def get_single_price(self) -> float:
        return self._single_price

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Item):
            return NotImplemented
        return (
            self.get_name() == other.get_name()
            and self.get_count() == other.get_count()
        )

    def __hash__(self) -> int:
        return hash((self.get_name(), self.get_count()))
