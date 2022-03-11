
from pyspark import SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.types import *
import sys


if __name__ == "__main__":
    app_name = "DataEngineering"

    conf = SparkConf()

    hdfs_host = "hdfs://namenode:8020"

    conf.set("hive.metastore.uris", "http://hive-metastore:9083")
    conf.set("spark.kerberos.access.hadoopFileSystem", hdfs_host)
    conf.set("spark.sql.warehouse.dir", f"{hdfs_host}/user/hive/warehouse")
    conf.set("hive.metastore.warehouse.dir",
             f"{hdfs_host}/user/hive/warehouse")
    conf.setMaster("local[*]")

    spark = (
        SparkSession.builder.appName(app_name)
        .config(conf=conf)
        .enableHiveSupport()
        .getOrCreate()
    )
    # ----------------------------------------------------------------------------------------
    num = int(sys.argv[1])

    path = "/airflow/data/poster_paths/path_"
    for i in range(num):
        df = spark.read.option("header", True).csv(path+str(i))
        df.write.saveAsTable('movies.posters', format='parquet', mode='append')

    