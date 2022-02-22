from app.database import IStoreDatabase, StoreDatabase
from app.factory import DefaultActorFactory, IActorFactory
from app.item import Item, OrdinaryPack
from app.printer import IPrinter, Printer
from app.receipt import ReceiptBuilder
from app.store import IStore, Store

PER_ZREPORT_CUSTOMER_COUNT = 100
PER_XREPORT_CUSTOMER_COUNT = 20


class StoreSimulator:
    factory: IActorFactory
    _store: IStore
    _printer: IPrinter

    def __init__(self) -> None:
        self.factory = DefaultActorFactory()
        manager = self.factory.create_manager()
        cashiers = [self.factory.create_cashier()]
        database = self.__init_database()
        self._store = Store(manager=manager, cashiers=cashiers, database=database)
        self._printer = Printer()

    def simulate_shift(self) -> None:
        customers_served = 0
        while True:
            self.__simulate_serve_customer()
            customers_served += 1
            shift_ended = self.__check_report(customers_served)
            if shift_ended:
                break

    def __simulate_serve_customer(self) -> None:
        cashier = self._store.get_available_cashier()
        customer = self.factory.create_customer(
            cashier=cashier, store_items=self._store.get_items()
        )
        cashier.open_receipt(receipt_builder=ReceiptBuilder())
        customer_items = customer.show_items()
        for item in customer_items:
            cashier.add_item_to_receipt(item)
        receipt = cashier.give_receipt()
        self._printer.print_receipt(
            receipt=receipt,
            catalogue=self._store.get_catalogue(),
            discounts=self._store.get_discounts(),
        )
        customer.pay_for_items()
        cashier.close_receipt()
        self._store.add_sold_items(receipt)

    def __check_report(self, customers_served: int) -> bool:
        if customers_served % PER_ZREPORT_CUSTOMER_COUNT == 0:
            if self._store.get_manager().answer_y_n_question():
                self._store.close_shift()
                return True
        elif customers_served % PER_XREPORT_CUSTOMER_COUNT == 0:
            if self._store.get_manager().answer_y_n_question():
                x_report = self._store.get_manager().make_X_report(
                    self._store.get_sold_items(), self._store.get_revenue()
                )
                self._printer.print_X_report(x_report)
        return False

    def __init_database(self) -> IStoreDatabase:
        database = StoreDatabase()
        items = {
            Item("Milk"): 4.99,
            Item("Mineral Water"): 3.00,
            Item("Bread"): 0.80,
            Item("Diapers"): 1.39,
            OrdinaryPack(6, "Beer"): 1.00,
            Item("Cheese"): 4.00,
            OrdinaryPack(10, "Tissues"): 1.50,
        }
        discounts = {
            (OrdinaryPack(10, "Tissues"),): 0.1,
            (Item("Bread"), Item("Cheese")): 0.05,
            (Item("Mineral Water"),): 0.5,
        }
        for item in items.items():
            database.add_item(item)

        for discount in discounts.items():
            database.add_discount(discount)

        return database
