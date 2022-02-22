import random
from typing import Dict, Iterator, Protocol, Tuple

from app.item import IItem


class IXReport(Protocol):
    def calculate_price(self) -> float:
        pass

    def __iter__(self) -> Iterator[Tuple[IItem, int]]:
        pass


class XReport:
    _sold_items: Dict[IItem, int]
    _revenue: float

    def __init__(self, sold_items: Dict[IItem, int], revenue: float):
        self._sold_items = sold_items
        self._revenue = revenue

    def calculate_price(self) -> float:
        return self._revenue

    def __iter__(self) -> Iterator[Tuple[IItem, int]]:
        return iter(list(self._sold_items.items()))


class IStoreManager(Protocol):
    def make_X_report(self, sold_items: Dict[IItem, int], revenue: float) -> IXReport:
        pass

    def answer_y_n_question(self) -> bool:
        pass


class StoreManager:
    def make_X_report(self, sold_items: Dict[IItem, int], revenue: float) -> IXReport:
        return XReport(sold_items, revenue)

    def answer_y_n_question(self) -> bool:
        return bool(random.getrandbits(1))
