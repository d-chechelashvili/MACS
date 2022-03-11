import os
import shutil  # to save it locally

import pyspark.sql.functions as f
import requests  # to get image from the web
from imdb import IMDb, IMDbError
from pyspark import SparkConf
from pyspark.sql import SparkSession, Window
from pyspark.sql.functions import udf
from pyspark.sql.types import *
import sys


def download_poster(imdb_id):
    poster_path = None
    try:
        ia = IMDb()
        movie = ia.get_movie(imdb_id)
        # Set up the image URL and filename
        if "cover url" not in movie:
            return poster_path
        image_url = movie["cover url"]
        filename = "/airflow/data/images/" + image_url.split("/")[-1]
        # Open the url image, set stream to True, this will return the stream content.
        r = requests.get(image_url, stream=True)

        # Check if the image was retrieved successfully
        if r.status_code == 200:
            # Set decode_content value to True, otherwise the downloaded image file's size will be zero.
            r.raw.decode_content = True

            # Open a local file with wb ( write binary ) permission.
            if not os.path.exists(filename):
                with open(filename, "wb") as f:
                    shutil.copyfileobj(r.raw, f)
                    poster_path = filename
            else:
                poster_path = filename
        else:
            print("Not 200")
    except IMDbError as e:
        print("FAILED", e)
    return poster_path


if __name__ == "__main__":
    if not os.path.exists("/airflow/data/images"):
        os.makedirs("/airflow/data/images")

    # ----------------------------------------------------------------------------------------

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

    path = hdfs_host + "/user/hive/warehouse/movies/movies"
    df = spark.read.option("header", True).parquet(path).drop("status", "title")
    df = df.withColumn("imdb_id", f.regexp_replace("imdb_id", r"^[t]*", ""))

    # ----------------------------------------------------------------------------------------

    range_start = int(sys.argv[1])
    range_end = int(sys.argv[2])
    window = Window.orderBy(f.col("movie_id"))
    df = df.withColumn("row", f.row_number().over(window))
    df = df.filter((f.col("row") >= range_start) & (f.col("row") <= range_end))
    df.persist()
    count = df.count()

    # ----------------------------------------------------------------------------------------

    downloader = udf(lambda x: download_poster(x), StringType())
    df = df.withColumn("poster_path", downloader(df.imdb_id)).drop("imdb_id", "row")
    df.write.option("header", True).csv(f"/airflow/data/poster_paths/path_{sys.argv[3]}", mode="overwrite")
