from datetime import datetime

from airflow.models import DAG
from airflow.providers.apache.spark.operators.spark_jdbc import SparkJDBCOperator
from airflow.providers.apache.spark.operators.spark_sql import SparkSqlOperator
from airflow.providers.apache.spark.operators.spark_submit import SparkSubmitOperator
from airflow.operators.trigger_dagrun import TriggerDagRunOperator


import os

with DAG(
    dag_id="fakestream_2",
    schedule_interval="* * * * *",
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=['FreeUni'],
    max_active_runs=1,
) as dag:
    submit_job_2 = SparkSubmitOperator(
        application="/airflow/jobs/fakestream_2_job.py", task_id="submit_job_2"
    )

    
    trigger = TriggerDagRunOperator(
        task_id="trigger_dagrun",
        trigger_dag_id="hdfs_sensor_2",
    )

    submit_job_2 >> trigger