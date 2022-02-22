from typing import Dict, Protocol, Tuple

from app.item import IItem
from app.receipt import IReceipt
from app.store_manager import IXReport


class IPrinter(Protocol):
    def print_receipt(
        self,
        receipt: IReceipt,
        catalogue: Dict[IItem, float],
        discounts: Dict[Tuple[IItem], float],
    ) -> None:
        pass

    def print_X_report(self, x_report: IXReport) -> None:
        pass


class Printer:
    def print_receipt(
        self,
        receipt: IReceipt,
        catalogue: Dict[IItem, float],
        discounts: Dict[Tuple[IItem], float],
    ) -> None:
        print("------------RECEIPT------------------")
        for item in receipt:
            name = item.get_name()
            count = item.get_count()
            price = catalogue[item]
            total = item.calculate_price(catalogue, discounts)
            print("Name:", name, "| Units", count, "| Price", price, "| Total:", total)
        print("\nSum:", "{:.2f}".format(receipt.calculate_price(catalogue, discounts)))

    def print_X_report(self, x_report: IXReport) -> None:
        print("------------X Report------------------")
        for item in x_report:
            print("Name:", item[0].get_name(), "| Sold:", item[0].get_count() * item[1])
        print("\nTotal Revenue: ", "{:.2f}".format(x_report.calculate_price()))
