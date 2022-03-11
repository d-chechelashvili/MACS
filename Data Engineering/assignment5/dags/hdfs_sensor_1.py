from datetime import datetime

from airflow.models import DAG

from airflow.providers.apache.hdfs.sensors.hdfs import HdfsRegexSensor
from airflow.providers.apache.spark.operators.spark_jdbc import SparkJDBCOperator
from airflow.providers.apache.spark.operators.spark_sql import SparkSqlOperator
from airflow.providers.apache.spark.operators.spark_submit import SparkSubmitOperator

import os
import re

with DAG(
    dag_id="hdfs_sensor_1",
    schedule_interval=None,
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=['FreeUni'],
    max_active_runs=1,
) as dag:
    sensor = HdfsRegexSensor(dag=dag,
        task_id="sense",
        filepath="/DataLake/fakestream_1",
        regex=re.compile("^(part-).*"),
        hdfs_conn_id="hdfs_conn")

    submit_job = SparkSubmitOperator(
        application="/airflow/jobs/DLToStaging_1.py",
        task_id="submit_job",
    )

    sensor >> submit_job