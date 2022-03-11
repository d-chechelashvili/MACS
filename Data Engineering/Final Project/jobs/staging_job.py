from pyspark.sql import SparkSession
from pyspark import SparkConf
import pyspark.sql.functions as f

from pyspark.sql.types import *



def replace_escape_chars(df, col_names):
    new_df = df
    for col_name in col_names:
        new_df = new_df.withColumn(col_name, f.regexp_replace(f.col(col_name), ': None', ': null'))\
          .withColumn(col_name, f.regexp_replace(f.col(col_name), "\\\\'", ""))\
          .withColumn(col_name, f.regexp_replace(f.col(col_name), "\\\\", ""))
    return new_df

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
	# ----------------------------------------------------------------------------------------

	cast_schema = ArrayType(StructType([StructField('cast_id', IntegerType(), True),
										StructField('character', StringType(), True),
										StructField('credit_id', StringType(), True),
										StructField('gender', IntegerType(), True),
										StructField('id', IntegerType(), True),
										StructField('name', StringType(), True),
										StructField('order', IntegerType(), True),
										StructField('profile_path', StringType(), True)]))

	creds = spark \
		.read \
		.option("quote", "\"") \
		.option("escape", "\"") \
		.csv('/airflow/data/credits.csv', header=True)

	creds = replace_escape_chars(creds, ["cast"])

	creds = creds \
		.withColumn('cast_members', f.from_json(creds.cast, cast_schema))

	creds = creds.select(creds.id, f.explode(creds['cast_members']))

	creds = creds.select(creds.id, creds["col"].getItem("id"), creds['col'].getItem('name'))
	creds = creds.withColumnRenamed("col.name", "name").withColumnRenamed("id", "movie_id")\
			.withColumnRenamed("col.id", "actor_id")

	actors = creds.select(creds.actor_id, creds.name).dropDuplicates()
	cast = creds.drop("name")

	# ----------------------------------------------------------------------------------------

	keywords_schema = ArrayType(StructType([StructField('id', IntegerType(), True),
											StructField('name', StringType(), True)]))

	keywords = spark \
		.read \
		.option("quote", "\"") \
		.option("escape", "\"") \
		.csv('/airflow/data/keywords.csv', header=True)

	keywords = replace_escape_chars(keywords, ["keywords"])

	keywords = keywords \
		.withColumn('kws', f.from_json(keywords.keywords, keywords_schema))

	keywords = keywords.select(keywords.id, f.explode(keywords['kws']))
	keywords = keywords.select(keywords.id, keywords["col"].getItem("name"))\
		.withColumnRenamed("id", "movie_id").withColumnRenamed("col.name", "keyword")

	# ----------------------------------------------------------------------------------------

	ratings = spark.read.csv('/airflow/data/ratings.csv', header=True).select("movieId", "rating").withColumnRenamed("movieId",
																									   "movie_id")

	# ----------------------------------------------------------------------------------------
	movie = spark \
		.read \
		.option("multiLine", "true") \
		.option("quote", '"') \
		.option("header", "true") \
		.option("escape", '"') \
		.option("wholeFile", True) \
		.csv("/airflow/data/movies_metadata.csv", header=True)

	movie = movie.where("imdb_id not in ('0', 'tt0113002', 'tt2423504', 'tt2622826')")
	movie = replace_escape_chars(movie, ['genres', 'production_companies', 'production_countries', 'spoken_languages'])
	movie = movie.drop("adult", "belongs_to_collection", "budget", "original_language", "poster_path", "tagline",
					   "video", "vote_average", "vote_count").withColumnRenamed("id", "movie_id")

	# ----------------------------------------------------------------------------------------

	genres_schema = ArrayType(StructType([StructField('id', IntegerType(), True),
										  StructField('name', StringType(), True)]))

	companies_schema = ArrayType(StructType([StructField('id', IntegerType(), True),
											 StructField('name', StringType(), True)]))

	countries_schema = ArrayType(StructType([StructField('iso_3166_1', StringType(), True),
											 StructField('name', StringType(), True)]))

	languages_schema = ArrayType(StructType([StructField('iso_639_1', StringType(), True),
											 StructField('name', StringType(), True)]))

	movie = movie.withColumn('genres', f.from_json(movie.genres, genres_schema))
	movie = movie.withColumn('production_companies', f.from_json(movie.production_companies, companies_schema))
	movie = movie.withColumn('production_countries', f.from_json(movie.production_countries, countries_schema))
	movie = movie.withColumn('spoken_languages', f.from_json(movie.spoken_languages, languages_schema))

	# ----------------------------------------------------------------------------------------

	companies = movie.select("movie_id", f.explode(movie["production_companies"]))
	companies = companies.select(companies.movie_id, companies["col"].getItem("name")).withColumnRenamed("col.name",
																										 "production_company_name")
	movie = movie.drop("production_companies")

	# ----------------------------------------------------------------------------------------

	countries = movie.select("movie_id", f.explode(movie["production_countries"]))
	countries = countries.select(countries.movie_id, countries["col"].getItem("name")).withColumnRenamed("col.name",
																										 "production_country_name")
	movie = movie.drop("production_countries")

	# ----------------------------------------------------------------------------------------

	genres = movie.select("movie_id", f.explode(movie["genres"]))
	genres = genres.select(genres.movie_id, genres["col"].getItem("name")).withColumnRenamed("col.name", "genre_name")
	movie = movie.drop("genres")

	# ----------------------------------------------------------------------------------------

	languages = movie.select("movie_id", f.explode(movie["spoken_languages"]))
	languages = languages.select(languages.movie_id, languages["col"].getItem("name")).withColumnRenamed("col.name",
																										 "spoken_language_name")
	movie = movie.drop("spoken_languages")

	# ----------------------------------------------------------------------------------------

	additional_info = movie.select("movie_id", "homepage", "original_title", "overview")
	movie = movie.drop("homepage", "original_title", "overview")

	# ----------------------------------------------------------------------------------------

	stats = movie.select("movie_id", "popularity", "release_date", "revenue", "runtime")
	movie = movie.drop("popularity", "release_date", "revenue", "runtime")

	# ----------------------------------------------------------------------------------------

	votes = ratings.groupBy("movie_id").agg(f.mean('rating'), f.count('rating')) \
		.withColumnRenamed("avg(rating)", "vote_average") \
		.withColumnRenamed("count(rating)", "vote_count")
	stats = stats.join(votes, "movie_id", 'left')

	# ----------------------------------------------------------------------------------------

	path = hdfs_host + "/user/hive/warehouse/"
	spark.sql(f"create schema if not exists movies location '{path}/movies'")

	keywords.write.saveAsTable('movies.keywords', format='parquet', mode='overwrite')
	ratings.write.saveAsTable('movies.ratings', format='parquet', mode='overwrite')
	cast.write.saveAsTable('movies.cast', format='parquet', mode='overwrite')
	actors.write.saveAsTable('movies.actors', format='parquet', mode='overwrite')
	movie.write.saveAsTable('movies.movies', format='parquet', mode='overwrite')
	companies.write.saveAsTable('movies.movie_companies', format='parquet', mode='overwrite')
	countries.write.saveAsTable('movies.movie_countries', format='parquet', mode='overwrite')
	genres.write.saveAsTable('movies.movie_genres', format='parquet', mode='overwrite')
	languages.write.saveAsTable('movies.movie_languages', format='parquet', mode='overwrite')
	stats.write.saveAsTable('movies.stats', format='parquet', mode='overwrite')
	additional_info.write.saveAsTable('movies.additional_info', format='parquet', mode='overwrite')

