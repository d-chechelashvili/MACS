from collections import defaultdict
from typing import DefaultDict, Dict, Protocol, Tuple

from app.item import IItem


class StoreDatabase:
    # Item - price
    _items: Dict[IItem, float]
    _sold_items: DefaultDict[IItem, int]
    # Item combination - discount value(percentage)
    _discount_list: Dict[Tuple[IItem], float]
    _revenue: float

    def __init__(self) -> None:
        self._items = {}
        self._sold_items = defaultdict(int)
        self._discount_list = {}
        self._revenue = 0.0

    def add_item(self, item: Tuple[IItem, float]) -> None:
        self._items[item[0]] = item[1]

    def add_discount(self, discount: Tuple[Tuple[IItem], float]) -> None:
        self._discount_list[discount[0]] = discount[1]

    def add_sold_item(self, item: IItem) -> None:
        self._sold_items[item] += 1

    def get_price(self, item: IItem) -> float:
        return self._items[item]

    def get_discounts(self) -> Dict[Tuple[IItem], float]:
        return self._discount_list.copy()

    def get_items(self) -> Dict[IItem, float]:
        return self._items.copy()

    def get_sold_items(self) -> DefaultDict[IItem, int]:
        return self._sold_items.copy()

    def set_revenue(self, amount: float) -> None:
        self._revenue = amount

    def get_revenue(self) -> float:
        return self._revenue

    def reset(self) -> None:
        self._revenue = 0
        self._sold_items = defaultdict(int)


class IStoreDatabase(Protocol):
    def add_item(self, item: Tuple[IItem, float]) -> None:
        pass

    def add_discount(self, discount: Tuple[Tuple[IItem], float]) -> None:
        pass

    def add_sold_item(self, item: IItem) -> None:
        pass

    def get_price(self, item: IItem) -> float:
        pass

    def get_discounts(self) -> Dict[Tuple[IItem], float]:
        pass

    def get_items(self) -> Dict[IItem, float]:
        pass

    def get_sold_items(self) -> DefaultDict[IItem, int]:
        pass

    def set_revenue(self, amount: float) -> None:
        pass

    def get_revenue(self) -> float:
        pass

    def reset(self) -> None:
        pass
