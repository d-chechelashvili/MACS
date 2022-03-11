from datetime import datetime
from datetime import timedelta

from airflow.models import DAG
from airflow.providers.apache.spark.operators.spark_jdbc import SparkJDBCOperator
from airflow.providers.apache.spark.operators.spark_sql import SparkSqlOperator
from airflow.providers.apache.spark.operators.spark_submit import SparkSubmitOperator
from airflow.operators.trigger_dagrun import TriggerDagRunOperator
from airflow.operators.python_operator import PythonOperator
#from airflow.sensors.external_task_sensor import ExternalTaskSensor

from external_task import ExternalTaskSensor



import os

with DAG(
    dag_id="filter_data",
    schedule_interval="*/2 * * * *",
    start_date=datetime(2021, 1, 1),
    catchup=False,
    dagrun_timeout=timedelta(seconds=180),
    tags=['FreeUni'],
    max_active_runs=2,
) as dag:
    external_task_sensor = ExternalTaskSensor(
        task_id='external_task_sensor',
        poke_interval=10,
        external_task_id='submit_job',
        external_dag_id=['hdfs_sensor_1', 'hdfs_sensor_2'],
        dag=dag
    )

    submit_job = SparkSubmitOperator(
        application="/airflow/jobs/filterData.py",
        task_id="submit_job",
    )

    external_task_sensor >> submit_job
