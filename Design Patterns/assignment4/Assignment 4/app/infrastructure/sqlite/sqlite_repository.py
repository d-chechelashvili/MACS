import sqlite3
from collections import defaultdict
from sqlite3 import Cursor
from typing import DefaultDict

from app.core.receipt.receipt import IReceipt
from app.core.store.item import IItem, Item


class SQLiteStoreRepository:
    def __init__(self) -> None:
        self._database_name = "store.db"
        with sqlite3.connect(self._database_name) as conn:
            cursor = conn.cursor()
            cursor.execute(
                """CREATE TABLE IF NOT EXISTS sold_items
                    (id INTEGER PRIMARY KEY NOT NULL,
                    name TEXT NOT NULL,
                    count INTEGER NOT NULL,
                    single_price FLOAT NOT NULL,
                    receipt_id INTEGER NOT NULL);"""
            )

    def _add_sold_item(self, item: IItem, receipt_id: int, cursor: Cursor) -> None:
        command = """ INSERT INTO sold_items
                        (name,count,single_price,receipt_id)
                        VALUES(?,?,?,?) """
        data = (
            item.get_name(),
            item.get_count(),
            item.get_single_price(),
            receipt_id,
        )
        cursor.execute(command, data)

    def get_sold_items(self) -> DefaultDict[IItem, int]:
        sold_items: DefaultDict[IItem, int] = defaultdict(int)
        with sqlite3.connect(self._database_name) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT * FROM sold_items")
            rows = cursor.fetchall()
            for item_tuple in rows:
                item = Item(
                    name=item_tuple[1], count=item_tuple[2], single_price=item_tuple[3]
                )
                sold_items[item] += 1
        return sold_items

    def get_revenue(self) -> float:
        revenue = 0.0
        sold_items = self.get_sold_items()
        for item, count in sold_items.items():
            revenue += item.calculate_price() * count
        return revenue

    def _get_max_receipt_id(self) -> int:
        max_id = -1
        with sqlite3.connect(self._database_name) as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT max(receipt_id) FROM sold_items")
            row = cursor.fetchone()
            if row[0] is not None:
                max_id = row[0]
        return max_id

    def add_closed_receipt(self, receipt: IReceipt) -> None:
        new_receipt_id = self._get_max_receipt_id() + 1
        with sqlite3.connect(self._database_name) as conn:
            cursor = conn.cursor()
            for item in receipt:
                self._add_sold_item(item, new_receipt_id, cursor)
            conn.commit()

    def get_receipt_count(self) -> int:
        return self._get_max_receipt_id() + 1
