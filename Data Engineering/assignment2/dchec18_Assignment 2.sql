
---------------------------------------- 1 ----------------------------------------

--ავტივრთე მოცემული csv ფაილები  /data_lake დირექტორიაში
CREATE SCHEMA staging;
CREATE DATABASE movies;

---------------------------------------- 2 ----------------------------------------

-- create table staging.users
CREATE TABLE staging.users
(
    user_id string,
    first_name string,
    last_name string,
    birth_date string, 
    country string
)
ROW FORMAT DELIMITED
FIELDS TERMINATED BY ','
STORED AS TEXTFILE;


-- fill staging.users
LOAD DATA INPATH 'hdfs:///data_lake/users.csv' INTO TABLE staging.users;



-- create table staging.ratings
CREATE TABLE staging.ratings
(
    user_id string,
    movie_id string,
    rating string,
    created_at string
)
ROW FORMAT DELIMITED
FIELDS TERMINATED BY ','
STORED AS TEXTFILE;


-- fill staging.ratings
LOAD DATA INPATH 'hdfs:///data_lake/ratings.csv' INTO TABLE staging.ratings;



-- create table staging.movies
CREATE TABLE staging.movies
(
    movie_id string,
    title string,
    genres string
)
ROW FORMAT DELIMITED
FIELDS TERMINATED BY ','
STORED AS TEXTFILE;


-- fill staging.movies
LOAD DATA INPATH 'hdfs:///data_lake/movies.csv' INTO TABLE staging.movies;



-- create table staging.tags
CREATE TABLE staging.tags
(
    user_id string,
    movie_id string,
    tag string,
    created_at string
)
ROW FORMAT DELIMITED
FIELDS TERMINATED BY ','
STORED AS TEXTFILE;


-- fill staging.tags
LOAD DATA INPATH 'hdfs:///data_lake/tags.csv' INTO TABLE staging.tags;



-- create table movies.users
CREATE EXTERNAL TABLE movies.users
(
    user_id int,
    first_name string,
    last_name string,
    birth_date date,
    registration_date date
)
PARTITIONED BY
    (country string, yearp int)
STORED AS PARQUET
LOCATION 'hdfs:///user/hive/warehouse/movies.db/partitioning/users_p';


set hive.exec.dynamic.partition = true;
set hive.exec.dynamic.partition.mode = nonstrict;


-- fill table movies.users
INSERT INTO movies.users
PARTITION (country, yearp)
SELECT 
    cast(u.user_id as int) user_id,
    u.first_name,
    u.last_name,
    to_date(from_unixtime(unix_timestamp(u.birth_date, 'M/dd/yyyy'))) birth_date,
    to_date(from_unixtime(unix_timestamp(u.registration_date, 'M/dd/yyyy'))) registration_date,
    u.country,
    cast(substr(u.registration_date, -4) as int) yearp
FROM staging.users u;



-- create table movies.ratings
CREATE EXTERNAL TABLE movies.ratings
(
    user_id int,
    movie_id int,
    rating double,
    created_at timestamp
)
STORED AS PARQUET;


-- fill table movies.ratings
INSERT INTO movies.ratings
SELECT 
    cast(r.user_id as int) user_id,
    cast(r.movie_id as int) movie_id,
    cast(r.rating as double) rating,
    from_unixtime(unix_timestamp(r.created_at, 'M/dd/yyyy HH:mm:ss')) created_at
FROM staging.ratings r;



-- create table movies.tags
CREATE EXTERNAL TABLE movies.tags
(
    user_id int,
    movie_id int,
    tags array<string>,
    created_at timestamp
)
STORED AS PARQUET;


-- fill table movies.tags
INSERT INTO movies.tags 
SELECT
    cast(t.user_id as int) user_id,
    cast(t.movie_id as int) movie_id,
    collect_set(tag) as tags,
    from_unixtime(unix_timestamp(t.created_at, 'M/dd/yyyy HH:mm')) as created_at
FROM staging.tags t
GROUP BY t.user_id, t.movie_id, t.created_at;



-- create table movies.genres
CREATE TABLE movies.genres
(
    genre_id int,
    genre_name string
)
STORED AS PARQUET;


-- fill table movies.genres
WITH genre_names as
(
    SELECT DISTINCT
        names 
    FROM staging.movies LATERAL VIEW explode(split(genres, '\\|')) exploded as names
)
INSERT INTO movies.genres
SELECT
    row_number() over(),
    *
FROM genre_names;



-- create table movies.movies
CREATE EXTERNAL TABLE movies.movies
(
    movie_id int,
    title string,
    year int,
    genres array<int>
)
STORED AS PARQUET;


-- fill table movies.movies
WITH movies_tmp as(
SELECT
    cast(m.movie_id as int) movie_id,
    substr(m.title, 1, length(m.title) - 6) as title,
    cast(substr(m.title, -5, 4) as int) as year, 
    split(m.genres, '\\|') as genres
FROM staging.movies m
),

-- tmp2 as(
--     SELECT
--         *,
--         split(mt.title, ':') parts 
--     FROM movies_tmp as mt
-- ),

-- tmp3 as(
--     SELECT
--         *,
--         trim(t2.parts[size(t2.parts)-1]) as article
--     FROM tmp2 as t2
--     WHERE trim(t2.parts[size(t2.parts)-1]) in ('The', 'A', 'An', 'Les', 'El', 'La')
-- )
--თეგები კი ამოაქვს ამას სწორად, მაგრამ მერე ამ თეგების შიგნით პოვნა და წაშლა ცალსახა რომ იყოს საკმაოდ რთული დასაწერია, მგონი
-- SELECT * from tmp3

exploded_movies as(
    SELECT 
        * 
    FROM movies_tmp LATERAL VIEW explode(genres) exploded as genre_name
),

tmp as(
SELECT 
    em.movie_id,
    em.title,
    em.year,
    collect_list(g.genre_id) as genres
FROM exploded_movies em 
LEFT JOIN genres g ON (em.genre_name = g.genre_name)
GROUP BY em.movie_id, em.title, em.year
)

INSERT INTO movies.movies
SELECT
    em.movie_id,
    CASE
        WHEN substring_index(em.title, ': The ', 1) != em.title THEN 'The ' || rtrim(substring_index(em.title, ': The ', 1))
        WHEN substring_index(em.title, ': A ', 1) != em.title THEN 'A ' || rtrim(substring_index(em.title, ': A ', 1))
        WHEN substring_index(em.title, ': An ', 1) != em.title THEN 'An ' || rtrim(substring_index(em.title, ': An ', 1))
        WHEN substring_index(em.title, ': Les ', 1) != em.title THEN 'Les ' || rtrim(substring_index(em.title, ': Les ', 1))
        WHEN substring_index(em.title, ': El ', 1) != em.title THEN 'El ' || rtrim(substring_index(em.title, ': El ', 1))
        WHEN substring_index(em.title, ': La ', 1) != em.title THEN 'La ' || rtrim(substring_index(em.title, ': La ', 1))
        ELSE rtrim(em.title)
    END as title,
    em.year,
    em.genres
FROM tmp as em;


---------------------------------------- 3 ----------------------------------------

-- 3.1)
SELECT 
    r.movie_id,
    year(r.created_at) as year,
    avg(r.rating) as avg_rating
FROM movies.ratings as r
GROUP BY r.movie_id, year(r.created_at);



-- 3.2)
SELECT 
    m.title,
    r1.average_rating,
    r1.rating_count 
FROM movies.movies m RIGHT JOIN 
    (
        SELECT 
            r.movie_id,
            avg(r.rating) as average_rating,
            count(*) as rating_count 
        FROM movies.ratings r
        GROUP BY r.movie_id
        HAVING average_rating > 4 AND rating_count >= 100
    ) r1 
ON m.movie_id = r1.movie_id;



-- 3.3)
-- (count(*) from movies.users) - (count(distinct user_id) from movies.ratings) ვერ გავმართე სინტაქსურად ამიტომ ასე ვწერ
SELECT 
    count(*)
FROM movies.users u
LEFT JOIN movies.ratings r ON (u.user_id = r.user_id)
WHERE r.user_id IS NULL;



-- 3.4)
WITH exploded_tags as(
    SELECT 
    *
    FROM movies.tags t LATERAL VIEW explode(t.tags) exploded_tags as tag
),

duplicates as(
    SELECT 
        et.user_id,
        et.movie_id,
        et.tag
    FROM exploded_tags as et
    GROUP BY et.user_id, et.movie_id, et.tag
    HAVING count(*) > 1
)

SELECT
    d.user_id, 
    concat_ws(' ', u.first_name, u.last_name) as user_name,
    d.movie_id,
    m.title
FROM duplicates d 
LEFT JOIN movies.users u ON (d.user_id = u.user_id)
LEFT JOIN movies.movies m ON (d.movie_id = m.movie_id);



-- 3.5)
WITH tag_sizes as
(
    SELECT 
        t.movie_id,
        size(t.tags) as sz,
        t.tags
    FROM movies.tags t
    ORDER BY sz DESC
),

top_3 as
(
    SELECT 
        t.movie_id,
        max(size(t.tags)) as sz
        
    FROM movies.tags t
    GROUP BY t.movie_id
    ORDER BY sz DESC
    LIMIT 3
),

full_top_3 as 
(
    SELECT 
        t3.movie_id,
        t3.sz,
        ts.tags
    FROM top_3 t3 
    INNER JOIN tag_sizes ts 
    ON (t3.movie_id = ts.movie_id AND t3.sz = ts.sz)
)

SELECT 
    m.title,
    ft3.sz,
    ft3.tags
FROM full_top_3 ft3
LEFT JOIN movies.movies m ON (ft3.movie_id = m.movie_id);



-- 3.6)
WITH top_films as (
SELECT 
    r.movie_id,
    sum(r.rating) as rating_sum,
    avg(r.rating) as rating_avg
FROM movies.ratings r
WHERE year(r.created_at) >= 2005 and year(r.created_at) <= 2015
GROUP BY r.movie_id
HAVING  rating_avg >= 4
ORDER BY rating_sum DESC
LIMIT 10
)

SELECT 
    m.title,
    tf.rating_sum,
    tf.rating_avg
FROM top_films tf
LEFT JOIN movies.movies m ON (tf.movie_id = m.movie_id);

