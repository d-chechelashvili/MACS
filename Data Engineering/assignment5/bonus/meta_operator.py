from airflow.models import BaseOperator
from airflow.sensors.base import apply_defaults
from airflow.utils.trigger_rule import TriggerRule
from meta_hook import MetaHook

class MetaOperator(BaseOperator):
    def __init__(self,
                 dag_id=None,
                 execution_date=None,
                 operator_type=None,
                 theme_color="#cc0000",
                 *args,
                 **kwargs
                 ):
        self.dag_id_=dag_id
        self.execution_date = execution_date
        self.operator_type = operator_type
        self.theme_color = theme_color
        super(MetaOperator, self).__init__(trigger_rule=TriggerRule.ALL_DONE, *args,**kwargs)


    def execute(self, context):
        MetaHook(
            dag_id=context['dag_run'].dag_id,
            execution_date=context['dag_run'].execution_date,
            operator_type=self.operator_type ).execute(context)