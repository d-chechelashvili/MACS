from pyspark.sql import SparkSession
from pyspark import SparkConf

from pyspark.sql.types import *
from pyspark.sql.functions import udf
import shutil # to save it locally
import os


if __name__ == "__main__":
	app_name = 'DataEngineering'

	conf = SparkConf()

	hdfs_host = 'hdfs://namenode:8020'

	conf.set("hive.metastore.uris", "http://hive-metastore:9083")
	conf.set("spark.kerberos.access.hadoopFileSystem", hdfs_host)
	conf.set("spark.sql.warehouse.dir", f"{hdfs_host}/user/hive/warehouse")
	conf.set("hive.metastore.warehouse.dir", f"{hdfs_host}/user/hive/warehouse")
	conf.setMaster("local[*]")
	conf.set("spark.sql.execution.arrow.enabled", "true")
	conf.set("spark.sql.execution.arrow.fallback.enabled.", "true")



	spark = SparkSession \
		.builder \
		.appName(app_name) \
		.config(conf=conf) \
		.enableHiveSupport() \
		.getOrCreate()
	# ----------------------------------------------------------------------------------------

	path = hdfs_host + "/user/hive/warehouse/movies/"
	posters = spark.read.option("header", True).parquet(path + "posters")
	genres = spark.read.option("header", True).parquet(path + "movie_genres")

	# ----------------------------------------------------------------------------------------

	df = genres.join(posters, "movie_id", "left")
	df = df.filter(df.poster_path.isNotNull())

	# ----------------------------------------------------------------------------------------
	def distribute_train_to_folders(genre_name, poster_path):
		img_path = "/airflow/data/image_classifier/train/" + genre_name
		if not os.path.exists(img_path):
			os.makedirs(img_path)
		shutil.copy(poster_path, img_path)
		return img_path

	def distribute_test_to_folders(genre_name, poster_path):
		img_path = "/airflow/data/image_classifier/test/" + genre_name
		if not os.path.exists(img_path):
			os.makedirs(img_path)
		shutil.copy(poster_path, img_path)
		return img_path

	# ----------------------------------------------------------------------------------------

	dfs = df.randomSplit([0.8, 0.2])
	train_distributer = udf(distribute_train_to_folders, StringType())
	test_distributer = udf(distribute_test_to_folders, StringType())
	dfs[0] = dfs[0].withColumn("TMP", train_distributer(dfs[0].genre_name, dfs[0].poster_path))
	dfs[1] = dfs[1].withColumn("TMP", test_distributer(dfs[1].genre_name, dfs[1].poster_path))
	dfs[0].select("TMP").show(1)
	dfs[1].select("TMP").show(1)
	# ----------------------------------------------------------------------------------------

	