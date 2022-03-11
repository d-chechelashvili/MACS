from datetime import datetime

from airflow.models import DAG
from airflow.providers.apache.spark.operators.spark_submit import SparkSubmitOperator
from airflow.models import Variable
import math



import os

DEFAULT_IMAGE_LIMIT = "500"

with DAG(
    dag_id="image_downloader",
    schedule_interval=None,
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=['FreeUni'],
    max_active_runs=1,
) as dag:

    num = int(Variable.get("thread_num", default_var="1"))
    limit = int(Variable.get("image_limit", default_var=DEFAULT_IMAGE_LIMIT))
    range_start = 1
    range_end = 0
    jobs = []
    step = int(math.ceil(limit/num))
    for i in range(num):
        range_end = min((range_start + step - 1), limit)
        download_job = SparkSubmitOperator(
            application="/airflow/jobs/download_job.py",
            task_id="download_job_" + str(i),
            application_args = [str(range_start), str(range_end), str(i)]
        )
        range_start = range_end + 1
        jobs.append(download_job)
    
    write_to_hive_job = SparkSubmitOperator(
            application="/airflow/jobs/write_to_hive.py",
            task_id="write_job",
            application_args = [str(num)],
            trigger_rule="all_success",
        )

    jobs >> write_to_hive_job
