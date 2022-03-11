from datetime import datetime

from airflow.models import DAG
from airflow.providers.apache.spark.operators.spark_submit import SparkSubmitOperator


import os

with DAG(
    dag_id="staging",
    schedule_interval=None,
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=['FreeUni'],
    max_active_runs=1,
) as dag:
    submit_job = SparkSubmitOperator(
        application="/airflow/jobs/staging_job.py",
        task_id="submit_job",
    )