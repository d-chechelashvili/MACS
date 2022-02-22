from typing import List

from app.item import IItem, Item, OrdinaryPack
from app.receipt import Receipt


def test_item_name() -> None:
    item = Item("Beer")
    assert item.get_name() == "Beer"


def test_item_count() -> None:
    item = Item("Beer")
    assert item.get_count() == 1


def test_item_price() -> None:
    item = Item("Beer")
    catalogue = {Item("Beer"): 1.00}
    assert item.calculate_price(catalogue, {}) == 1.00


def test_item_price_with_discount() -> None:
    item = Item("Beer")
    catalogue = {Item("Beer"): 1.00}
    discounts = {(Item("Beer"),): 0.5}
    assert item.calculate_price(catalogue, discounts) == 0.5


def test_pack_name() -> None:
    pack = OrdinaryPack(5, "Tissue")
    assert pack.get_name() == "Tissue"


def test_pack_count() -> None:
    pack = OrdinaryPack(5, "Tissue")
    assert pack.get_count() == 5


def test_pack_price() -> None:
    pack = OrdinaryPack(5, "Tissue")
    catalogue = {OrdinaryPack(5, "Tissue"): 1.00}
    assert pack.calculate_price(catalogue, {}) == 5.00


def test_pack_price_with_discount() -> None:
    pack = OrdinaryPack(5, "Tissue")
    catalogue = {OrdinaryPack(5, "Tissue"): 1.00}
    discounts = {(OrdinaryPack(5, "Tissue"),): 0.5}
    assert pack.calculate_price(catalogue, discounts) == 2.5


def test_receipt() -> None:
    items: List[IItem] = [Item("Milk"), Item("Lemonade"), OrdinaryPack(7, "Beer")]
    receipt = Receipt(items)
    assert items == receipt.get_items()


def test_receipt_calculate_price() -> None:
    items: List[IItem] = [Item("Milk"), Item("Lemonade"), OrdinaryPack(7, "Beer")]
    receipt = Receipt(items)
    catalogue = {
        Item("Milk"): 3.00,
        Item("Lemonade"): 2.00,
        OrdinaryPack(2, "Beer"): 1.00,
    }
    price = receipt.calculate_price(catalogue, {})
    assert price == 12.00


def test_receipt_calculate_price_with_discount() -> None:
    items: List[IItem] = [Item("Milk"), Item("Lemonade"), OrdinaryPack(7, "Beer")]
    receipt = Receipt(items)
    catalogue = {
        Item("Milk"): 3.00,
        Item("Lemonade"): 2.00,
        OrdinaryPack(2, "Beer"): 1.00,
    }
    discounts = {
        (Item("Milk"),): 1,
        (Item("Lemonade"), OrdinaryPack(7, "Beer")): 0.5,
    }
    price = receipt.calculate_price(catalogue, discounts)
    assert price == (3 * 0 + 2.00 * 0.5 + 7 * 1 * 0.5)
