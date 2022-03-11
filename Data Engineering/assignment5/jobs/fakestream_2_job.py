from pyspark import SparkContext, SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.functions import rand
from random import randint
from pyspark.sql.functions import rand 
from pyspark.sql.functions import row_number,lit
from pyspark.sql.window import Window
from pyspark.sql.types import StringType
from pyspark.sql.functions import concat, col, lit


def random_row(df, num):
    w = Window().orderBy(lit('A'))
    df = df.withColumn("id", row_number().over(w))
    df = df.withColumn("id", df.id.cast(StringType()))
    str = "_" + num
    df = df.withColumn("id", concat(col("id"), lit(str)))
    return df.orderBy(rand()).limit(1)


if __name__ == "__main__":
	app_name = 'DataEngineering'

	conf = SparkConf()

	hdfs_host = 'hdfs://namenode:8020'

	conf.set("hive.metastore.uris", "http://hive-metastore:9083")
	conf.set("spark.kerberos.access.hadoopFileSystem", hdfs_host)
	conf.set("spark.sql.warehouse.dir", f"{hdfs_host}/user/hive/warehouse")
	conf.set("hive.metastore.warehouse.dir", f"{hdfs_host}/user/hive/warehouse")
	conf.setMaster("local[*]")

	spark = SparkSession \
		.builder \
		.appName(app_name) \
		.config(conf=conf) \
		.enableHiveSupport() \
		.getOrCreate()

	sc = spark.sparkContext

	num = "2"#str(randint(1,2))
	path = "/airflow/data/pokemon" + num + ".csv"
	df = spark.read.option("header", True).csv(path)

	random_row_df = random_row(df, num)

	data_lake_path = hdfs_host + "/DataLake/fakestream_2"

	random_row_df.write.parquet(data_lake_path, mode="append")
