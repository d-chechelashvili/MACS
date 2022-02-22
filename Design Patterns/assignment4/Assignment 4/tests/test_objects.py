from typing import List

from app.core.receipt.receipt import Receipt
from app.core.store.item import IItem, Item


def test_item_name() -> None:
    item = Item("Beer", 1.0, 1)
    assert item.get_name() == "Beer"


def test_item_count() -> None:
    item = Item("Beer", 1.0, 1)
    assert item.get_count() == 1


def test_single_item_price() -> None:
    item = Item(name="Beer", count=1, single_price=1.00)
    assert item.calculate_price() == 1.00


def test_multiple_item_price() -> None:
    item = Item(name="Beer", count=3, single_price=1.00)
    assert item.calculate_price() == 3.00


def test_receipt_get() -> None:
    items: List[IItem] = [
        Item(name="Milk", single_price=1.1, count=1),
        Item(name="Lemonade", single_price=1.2, count=1),
        Item(name="Beer", single_price=1.2, count=7),
    ]
    receipt = Receipt(items)
    assert items == receipt.get_items()


def test_receipt_calculate_price() -> None:
    items: List[IItem] = [
        Item(name="Milk", single_price=8.4, count=1),
        Item(name="Lemonade", single_price=1.2, count=1),
        Item(name="Beer", single_price=1.2, count=7),
    ]
    receipt = Receipt(items)
    price = receipt.calculate_price()
    assert price == 18.00
