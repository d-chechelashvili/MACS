from abc import ABC, abstractmethod
from typing import Dict, List, Protocol, Tuple


class ICalculatePrice(Protocol):
    def calculate_price(self, catalogue: Dict["IItem", float]) -> float:
        pass


class IItem(Protocol):
    def calculate_price(
        self, catalogue: Dict["IItem", float], discounts: Dict[Tuple["IItem"], float]
    ) -> float:
        pass

    def get_name(self) -> str:
        pass

    def get_count(self) -> int:
        pass

    def __eq__(self, other: object) -> bool:
        pass

    def __hash__(self) -> int:
        pass


class Item:
    _name: str
    _count: int = 1

    def __init__(self, name: str):
        self._name = name

    def calculate_price(
        self, catalogue: Dict[IItem, float], discounts: Dict[Tuple[IItem], float]
    ) -> float:
        price = 0.0
        if self in catalogue:
            price = catalogue[self]
        else:
            for item in catalogue:
                if item.get_name() == self.get_name():
                    price = catalogue[item]

        if (self,) in discounts:
            price -= price * discounts[(self,)]

        return price

    def get_name(self) -> str:
        return self._name

    def get_count(self) -> int:
        return self._count

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Item):
            return NotImplemented
        return (
            self.get_name() == other.get_name()
            and self.get_count() == other.get_count()
        )

    def __hash__(self) -> int:
        return hash((self.get_name(), self.get_count()))


# using composite pattern because pack is composite of items
# using factory method for demonstration
class Pack(ABC):
    _items: List["IItem"]
    _name: str

    def __init__(self, number: int, name: str) -> None:
        self._name = name
        self._items = []
        self._create_items(number)

    def _create_items(self, number: int) -> None:
        for _ in range(number):
            item = self._make_item()
            self._items.append(item)

    def calculate_price(
        self, catalogue: Dict[IItem, float], discounts: Dict[Tuple[IItem], float]
    ) -> float:
        price = sum(i.calculate_price(catalogue, discounts) for i in self._items)
        if (self,) in discounts:
            price -= price * discounts[(self,)]
        return price

    def get_name(self) -> str:
        return self._name

    def get_count(self) -> int:
        return len(self._items)

    @abstractmethod
    def _make_item(self) -> "IItem":
        pass

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Pack):
            return NotImplemented
        return (
            self.get_name() == other.get_name()
            and self.get_count() == other.get_count()
        )

    def __hash__(self) -> int:
        return hash((self.get_name(), self.get_count()))


class OrdinaryPack(Pack):
    def _make_item(self) -> "IItem":
        return Item(self._name)
