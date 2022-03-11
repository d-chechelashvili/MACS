from pyspark import SparkContext, SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.types import *
from pyspark.sql.functions import countDistinct
from pyspark.sql.window import Window
from pyspark.sql.functions import col, row_number
from pyspark.sql.functions import concat_ws
from pyspark.sql.functions import lit

import sys

if __name__ == "__main__":

    app_name = 'DataEngineering'

    conf = SparkConf()

    hdfs_host = 'hdfs://namenode:8020'

    conf.set("hive.metastore.uris", "http://hive-metastore:9083")
    conf.set("spark.kerberos.access.hadoopFileSystem", hdfs_host)
    conf.set("spark.sql.warehouse.dir", f"{hdfs_host}/user/hive/warehouse")
    conf.set("hive.metastore.warehouse.dir",
             f"{hdfs_host}/user/hive/warehouse")
    conf.setMaster("local[*]")

    spark = SparkSession \
        .builder \
        .appName(app_name) \
        .config(conf=conf) \
        .enableHiveSupport()\
        .getOrCreate()

    sc = spark.sparkContext

    path = hdfs_host + "/user/hive/warehouse/staging/filteredData"
    df = spark.read.option("header", True).parquet(path)
    group_col = sys.argv[1]

    df = df.withColumn("is_legendary", df.is_legendary.cast(IntegerType()))
    df_res = df.groupBy(group_col).sum('is_legendary')
    df_res = df.groupBy(group_col).agg(
        countDistinct("abilities")).join(df_res, group_col)

    against_cols = [
        col for col in df.schema.names if col.startswith('against_')]
    mean_min_max_cols = ['attack', 'base_egg_steps', 'base_happiness', 'base_total',
                         'capture_rate', 'defense', 'height_m', 'hp', 'sp_attack', 'sp_defense', 'speed', 'weight_kg']

    for col_name in against_cols:
        df = df.withColumn(col_name, col(col_name).cast(IntegerType()))
        df_tmp = df.groupBy(group_col).sum(col_name)
        df_res = df_tmp.withColumnRenamed(
            df_tmp.schema.names[1], col_name).join(df_res, group_col)

    for col_name in mean_min_max_cols:
        df = df.withColumn(col_name, col(col_name).cast(IntegerType()))

        windowDept = Window.partitionBy(
            group_col).orderBy(col(col_name).desc())
        df_tmp = df.withColumn("row", row_number().over(windowDept)) \
            .filter(col("row") == 1).drop("row") \
            .withColumn("max_" + col_name, concat_ws('_', col(col_name), col("name"))).select(group_col, "max_" + col_name)
        df_res = df_res.join(df_tmp, group_col)

        windowDept = Window.partitionBy(group_col).orderBy(col(col_name).asc())
        df_tmp = df.withColumn("row", row_number().over(windowDept)) \
            .filter(col("row") == 1).drop("row") \
            .withColumn("min_" + col_name, concat_ws('_', col(col_name), col("name"))).select(group_col, "min_" + col_name)
        df_res = df_res.join(df_tmp, group_col)
        df_tmp = df.groupBy(group_col).mean(col_name)
        df_res = df_tmp.withColumnRenamed(
            df_tmp.schema.names[1], "mean_" + col_name).join(df_res, group_col)

    func_list = {}
    for col_name in against_cols:
        func_list[col_name] = "sum"

    for col_name in mean_min_max_cols:
        func_list["min_" + col_name] = "min"
        func_list["max_" + col_name] = "max"
        func_list["mean_" + col_name] = "mean"

    df_tmp = df_res.agg(func_list)
    for name in df_tmp.schema.names:
        df_tmp = df_tmp.withColumnRenamed(
            name, name[name.find('(')+1:name.find(')')])

    df_tmp_1 = df_res.agg(countDistinct("count(DISTINCT abilities)")).withColumnRenamed(
        "count(DISTINCT count(DISTINCT abilities))", "count(DISTINCT abilities)")
    df_tmp_2 = df_res.agg({"sum(is_legendary)": 'sum'}).withColumnRenamed(
        "sum(sum(is_legendary))", "sum(is_legendary)")

    df_tmp = df_tmp.withColumn(group_col, lit('total'))
    df_tmp_1 = df_tmp_1.withColumn(group_col, lit('total'))
    df_tmp_2 = df_tmp_2.withColumn(group_col, lit('total'))
    df_tmp = df_tmp.join(df_tmp_1, group_col)
    df_tmp = df_tmp.join(df_tmp_2, group_col)

    df_res = df_res.union(df_tmp.select(df_res.columns))

    df_res.coalesce(1).write.csv("/airflow/data/reports", mode="append", header=True)
