from pyspark import SparkContext, SparkConf
from pyspark.sql import SparkSession


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
		.enableHiveSupport()\
		.getOrCreate()

	spark._jsc.hadoopConfiguration().set('fs.defaultFS', hdfs_host)

	path = "/DataLake/fakestream_2"
	filename = None
	for name in list_files(path):
		if name.startswith("hdfs://namenode:8020/DataLake/fakestream_2/part-") and name.endswith(".parquet"):
			filename = name
			break
	
	df = spark.read.parquet(filename)
	staging_path = hdfs_host + "/user/hive/warehouse/staging"
	df.write.parquet(staging_path, mode="append")
	
	
	fs, path_type = hdp_conf()
	path = path_type(filename[len("hdfs://namenode:8020"):])
	fs.delete(path, False)
