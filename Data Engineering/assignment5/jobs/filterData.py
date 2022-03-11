from pyspark import SparkContext, SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.functions import rand
from random import randint
from pyspark.sql.functions import rand 
from pyspark.sql.functions import row_number,lit
from pyspark.sql.window import Window
from pyspark.sql.types import StringType
from pyspark.sql.functions import concat, col, lit


def hdp_conf():
    hdp = spark._jvm.org.apache.hadoop.fs
    fs = hdp.FileSystem.get(spark._jsc.hadoopConfiguration())
    
    return (fs, hdp.Path)


def list_files(path: str):
    fs, path_type = hdp_conf()
    
    path = path_type(path)
    
    files = fs.listFiles(path, False)
    
    result = []
    
    while files.hasNext():
        item = files.next()
        result.append(item.getPath().toString())
        
    return result

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

	spark._jsc.hadoopConfiguration().set('fs.defaultFS', hdfs_host)

	path = "/user/hive/warehouse/staging"
	parquet_names = [name for name in list_files(path) if name.endswith(".parquet")]
	df = spark.read.parquet(parquet_names[0])
	for name in parquet_names[0:]:
		df_1 = spark.read.parquet(name)
		df = df.union(df_1)
	
	df = df.dropDuplicates(['id'])
	df = df.drop('id')
	df.write.option('path', f'{path}/filteredData').saveAsTable('pokemon', format='parquet',  mode='overwrite')

