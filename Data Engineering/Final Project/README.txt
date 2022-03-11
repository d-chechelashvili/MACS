ერთი მთლიანი პროცესი თუ უნდა იყოს ეს დავალება შემეძლო ერთ დაგში გამეერთიანებინა ტასკები უბრალოდ
ან დამეტრიგერებინა ერთი დაგის მორჩენისას მეორე, მაგრამ ასე ვტოვებ, რადგან გასატესტად უფრო მოსახერხებელია წესით.


Airflow Connections-ში ჩავანაცვლე spark_default-ში host - local; Extra - { "deploy_mode": "client", "spark.root.logger": "ALL", "verbose": "true"}
Variables - thread_num : 5, image_limit : 2000

scheduler და webserver კონტეინერებში გავუშვი შემდეგი ბრძანებები:
pip install --upgrade pip
pip install IMDbPY
pip install tensorflow
apt-get update && apt-get install -y python3-opencv
pip install opencv-python

predict.py-ს გაშვების ინსტრუქცია:

1)predict.py უნდა იყოს data ფოლდერში

2) ვუშვებ scheduler-ის CLI-დან შემდეგნაირად
    2.1) cd data/
    2.2) python predict.py path_to_image