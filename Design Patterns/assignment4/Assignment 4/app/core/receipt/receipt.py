from typing import Iterator, List, Optional, Protocol

from app.core.store.item import IItem


class IReceipt(Protocol):
    def calculate_price(self) -> float:
        pass

    def get_items(self) -> List[IItem]:
        pass

    def __iter__(self) -> Iterator[IItem]:
        pass

    def __str__(self) -> str:
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
    _items: List[IItem]

    def __init__(self, items: List[IItem]):
        self._items = items

    def calculate_price(self) -> float:
        price = 0.0
        for item in self:
            price += item.calculate_price()
        return price

    def get_items(self) -> List[IItem]:
        return self._items.copy()

    def __iter__(self) -> Iterator[IItem]:
        return iter(self._items)

    def __str__(self) -> str:
        receipt_str = ""
        receipt_str += "Name             | Units  | Price  | Total  |\n"
        receipt_str += "-----------------|--------|--------|--------|\n"
        for item in self:
            name = item.get_name()
            price = item.get_single_price()
            count = item.get_count()
            total = item.calculate_price()
            receipt_str += name + " "
            receipt_str += (16 - len(name)) * " " + "| "
            receipt_str += str(count) + " "
            receipt_str += (6 - len(str(count))) * " " + "| "
            receipt_str += str(price) + " "
            receipt_str += (6 - len(str(price))) * " " + "| "
            receipt_str += str(round(total, 2)) + " "
            receipt_str += (6 - len(str(round(total, 2)))) * " " + "| "
            receipt_str += "\n"
        receipt_str += "\n\nSum: " + str(round(self.calculate_price(), 2)) + "\n"
        return receipt_str
