from datetime import datetime

from airflow.models import DAG
from airflow.providers.apache.spark.operators.spark_jdbc import SparkJDBCOperator
from airflow.providers.apache.spark.operators.spark_sql import SparkSqlOperator
from airflow.providers.apache.spark.operators.spark_submit import SparkSubmitOperator
from airflow.operators.trigger_dagrun import TriggerDagRunOperator


import os

with DAG(
    dag_id="pokemon_report",
    schedule_interval="*/15 * * * *",
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=['FreeUni'],
    max_active_runs=1,
) as dag:
    submit_job = SparkSubmitOperator(
        application="/airflow/jobs/pokemon_report_job.py",
        task_id="submit_job",
        application_args=['{{ dag_run.conf["message"] if (dag_run and "message" in dag_run.conf) else "type1" }}'],
    )