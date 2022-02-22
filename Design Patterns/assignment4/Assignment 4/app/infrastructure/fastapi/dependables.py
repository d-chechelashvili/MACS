from starlette.requests import Request

from app.core.facade import StoreService


def get_core(request: Request) -> StoreService:
    core: StoreService = request.app.state.core
    return core
