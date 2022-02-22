from fastapi import FastAPI

from app.core.facade import StoreService
from app.infrastructure.fastapi.receipt import receipt_api
from app.infrastructure.fastapi.report import report_api
from app.infrastructure.sqlite.sqlite_repository import SQLiteStoreRepository


def setup() -> FastAPI:
    app = FastAPI()
    app.include_router(receipt_api)
    app.include_router(report_api)
    store_repository = SQLiteStoreRepository()
    app.state.core = StoreService.create(store_repository=store_repository)
    return app
