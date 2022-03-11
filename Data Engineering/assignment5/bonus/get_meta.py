from airflow.utils.db import create_session
from airflow.models.taskinstance import TaskInstance

def get_meta(dag_id, execution_date, operator_type):
    with create_session() as session:
        res = session.query(TaskInstance).filter(TaskInstance.dag_id == dag_id,
                                                 TaskInstance.execution_date == execution_date,
                                                 TaskInstance.operator == operator_type).all()
    return [a.__dict__ for a in res]