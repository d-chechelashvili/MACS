from typing import Dict, Iterator, List, Optional, Protocol, Tuple

from app.item import IItem


# Using Builder, Iterator and composite patterns
class IReceipt(Protocol):
    def calculate_price(
        self, catalogue: Dict["IItem", float], discounts: Dict[Tuple[IItem], float]
    ) -> float:
        pass

    def get_items(self) -> List[IItem]:
        pass

    def __iter__(self) -> Iterator[IItem]:
        pass


class IReceiptBuilder(Protocol):
    def with_item(self, item: IItem) -> None:
        pass

    def build(self) -> IReceipt:
        pass

    def clear(self) -> None:
        pass


class ReceiptBuilder:
    def __init__(self, kwargs: Optional[List[IItem]] = None):
        self._kwargs = kwargs or []

    def with_item(self, item: IItem) -> None:
        self._kwargs.append(item)

    def build(self) -> IReceipt:
        return Receipt(self._kwargs)

    def clear(self) -> None:
        self._kwargs = []


class Receipt:
    def __init__(self, items: List[IItem]):
        self._items = items

    def calculate_price(
        self, catalogue: Dict["IItem", float], discounts: Dict[Tuple[IItem], float]
    ) -> float:
        from itertools import chain, combinations

        powerset = list(
            chain.from_iterable(
                combinations(self._items, size)
                for size in range(2, len(self._items) + 1)
            )
        )
        price = sum(i.calculate_price(catalogue, discounts) for i in self._items)
        for set in powerset:
            if set in discounts:
                price -= (
                    sum(i.calculate_price(catalogue, discounts) for i in set)
                    * discounts[set]  # type: ignore
                )
        return price

    def get_items(self) -> List[IItem]:
        return self._items.copy()

    def __iter__(self) -> Iterator[IItem]:
        return iter(self._items)
