from fastapi import APIRouter, Depends

from app.core.facade import StoreService
from app.core.store.store_interactor import XReportResponse
from app.infrastructure.fastapi.dependables import get_core

report_api = APIRouter()


@report_api.get("/report/x_report")
def make_x_report(core: StoreService = Depends(get_core)) -> XReportResponse:
    return core.make_x_report()
