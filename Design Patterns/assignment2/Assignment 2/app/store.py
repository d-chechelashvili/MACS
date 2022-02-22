from typing import Dict, List, Protocol, Tuple

from app.cashier import ICashier
from app.database import IStoreDatabase
from app.item import IItem
from app.receipt import IReceipt
from app.store_manager import IStoreManager


class Store:
    _manager: IStoreManager
    _cashiers: List[Tuple[ICashier, bool]]
    _database: IStoreDatabase

    def __init__(
        self,
        manager: IStoreManager,
        cashiers: List[ICashier],
        database: IStoreDatabase,
    ) -> None:
        self._manager = manager
        self._cashiers = [(cashier, True) for cashier in cashiers]
        self._database = database

    def get_manager(self) -> IStoreManager:
        return self._manager

    def get_available_cashier(self) -> ICashier:
        for cashierInfo in self._cashiers:
            if cashierInfo[1]:
                return cashierInfo[0]
        return self._cashiers[0][0]

    def get_items(self) -> List[IItem]:
        return list(self._database.get_items().keys())

    def get_discounts(self) -> Dict[Tuple[IItem], float]:
        discounts: Dict[Tuple[IItem], float] = self._database.get_discounts()
        return discounts

    def get_catalogue(self) -> Dict[IItem, float]:
        catalogue: Dict[IItem, float] = self._database.get_items()
        return catalogue

    def add_sold_items(self, closed_receipt: IReceipt) -> None:
        db = self._database
        for item in closed_receipt:
            db.add_sold_item(item)
        db.set_revenue(
            db.get_revenue()
            + closed_receipt.calculate_price(db.get_items(), db.get_discounts())
        )

    def close_shift(self) -> None:
        for cashierInfo in self._cashiers:
            cashierInfo[0].make_Z_report(self._database)

    def get_sold_items(self) -> Dict[IItem, int]:
        sold_items: Dict[IItem, int] = self._database.get_sold_items()
        return sold_items

    def get_revenue(self) -> float:
        revenue: float = self._database.get_revenue()
        return revenue


class IStore(Protocol):
    def get_manager(self) -> IStoreManager:
        pass

    def get_available_cashier(self) -> ICashier:
        pass

    def get_catalogue(self) -> Dict[IItem, float]:
        pass

    def get_items(self) -> List[IItem]:
        pass

    def get_discounts(self) -> Dict[Tuple[IItem], float]:
        pass

    def add_sold_items(self, closed_receipt: IReceipt) -> None:
        pass

    def get_sold_items(self) -> Dict[IItem, int]:
        pass

    def get_revenue(self) -> float:
        pass

    def close_shift(self) -> None:
        pass
