from typing import Tuple

import pytest

from app.database import IStoreDatabase, StoreDatabase
from app.item import IItem, Item, OrdinaryPack


@pytest.fixture
def db() -> IStoreDatabase:
    return StoreDatabase()


def test_create(db: IStoreDatabase) -> None:
    assert db.get_items() == {}
    assert db.get_revenue() == 0
    assert db.get_discounts() == {}
    assert db.get_sold_items() == {}


def test_item(db: IStoreDatabase) -> None:
    item = Item("Sour Milk")
    item_with_price = (item, 1.5)
    db.add_item(item_with_price)
    assert item in db.get_items()
    assert db.get_items()[item] == 1.5
    pack = OrdinaryPack(5, "Cigarette")
    pack_with_price = (pack, 1.0)
    db.add_item(pack_with_price)
    assert pack in db.get_items()
    assert db.get_items()[pack] == 1.0


def test_sold_item(db: IStoreDatabase) -> None:
    item = Item("Sour Milk")
    db.add_sold_item(item)
    assert item in db.get_sold_items()
    assert db.get_sold_items()[item] == 1
    item1 = Item("Sour Milk")
    assert item1 in db.get_sold_items()
    db.add_sold_item(item1)
    assert db.get_sold_items()[item] == 2
    pack = OrdinaryPack(5, "Cigarette")
    db.add_sold_item(pack)
    assert pack in db.get_sold_items()
    assert db.get_sold_items()[pack] == 1


def test_discount(db: IStoreDatabase) -> None:
    discount = ((Item("Sour Milk"),), 0.1)
    db.add_discount(discount)
    assert discount[0] in db.get_discounts()
    assert db.get_discounts()[discount[0]] == discount[1]
    discount1: Tuple[Tuple[IItem, ...], float] = (
        (Item("Sour Milk"), OrdinaryPack(5, "Beer")),
        0.2,
    )
    db.add_discount(discount1)
    assert discount1[0] in db.get_discounts()
    assert db.get_discounts()[discount1[0]] == discount1[1]


def test_revenue(db: IStoreDatabase) -> None:
    revenue = 5.0
    db.set_revenue(revenue)
    assert db.get_revenue() == revenue
    revenue_change = 3.0
    db.set_revenue(db.get_revenue() + revenue_change)
    assert db.get_revenue() == revenue + revenue_change


def test_reset(db: IStoreDatabase) -> None:
    test_item(db)
    test_discount(db)
    test_revenue(db)
    test_sold_item(db)
    db.reset()
    assert db.get_revenue() == 0
    assert db.get_sold_items() == {}
