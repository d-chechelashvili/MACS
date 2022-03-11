from datetime import datetime

from airflow.models import DAG
from airflow.providers.apache.spark.operators.spark_submit import SparkSubmitOperator
from airflow.operators.python import PythonOperator


import os

with DAG(
    dag_id="train_model",
    schedule_interval=None,
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=['FreeUni'],
    max_active_runs=1,
) as dag:
    prepare_data = SparkSubmitOperator(
        application="/airflow/jobs/prepare_data_job.py",
        task_id="prepare_data",
    )

    def train_model():
        from keras.models import Sequential
        from keras.layers import Dense, Conv2D , MaxPool2D , Flatten , Dropout 
        from keras.preprocessing.image import ImageDataGenerator
        from tensorflow.keras.optimizers import Adam
        import tensorflow as tf
        import cv2
        import os
        import numpy as np
        
        labels = ['Action', 'Adventure', 'Animation', 'Comedy', 'Crime', 'Documentary', 'Drama', 'Family', 'Fantasy', 'Foreign', 
         'History', 'Horror', 'Music', 'Mystery', 'Romance', 'Science Fiction', 'Thriller', 'TV Movie', 'War', 'Western']

        def get_data(data_dir):
            data = [] 
            for label in labels: 
                path = os.path.join(data_dir, label)
                if not os.path.exists(path):
                    print("PATH DOESNT EXIST:", path)
                    continue
                class_num = labels.index(label)
                for img in os.listdir(path):
                    try:
                        img_arr = cv2.imread(os.path.join(path, img))[...,::-1] #convert BGR to RGB format
                        resized_arr = cv2.resize(img_arr, (101, 150)) # Reshaping images to preferred size
                        data.append([resized_arr, class_num])
                    except Exception as e:
                        print(e)
            return np.array(data)

        train = get_data('/airflow/data/image_classifier/train')
        test = get_data('/airflow/data/image_classifier/test')

        x_train = []
        y_train = []
        x_test = []
        y_test = []

        for feature, label in train:
            x_train.append(feature)
            y_train.append(label)

        for feature, label in test:
            x_test.append(feature)
            y_test.append(label)

        # Normalize the data
        x_train = np.array(x_train) / 255
        x_test = np.array(x_test) / 255

        x_train.reshape(-1, 150, 101, 1)
        y_train = np.array(y_train)

        x_test.reshape(-1, 150, 101, 1)
        y_test = np.array(y_test)

        datagen = ImageDataGenerator(
        featurewise_center=False,  # set input mean to 0 over the dataset
        samplewise_center=False,  # set each sample mean to 0
        featurewise_std_normalization=False,  # divide inputs by std of the dataset
        samplewise_std_normalization=False,  # divide each input by its std
        zca_whitening=False,  # apply ZCA whitening
        rotation_range = 30,  # randomly rotate images in the range (degrees, 0 to 180)
        zoom_range = 0.2, # Randomly zoom image 
        width_shift_range=0.1,  # randomly shift images horizontally (fraction of total width)
        height_shift_range=0.1,  # randomly shift images vertically (fraction of total height)
        horizontal_flip = True,  # randomly flip images
        vertical_flip=False)  # randomly flip images


        datagen.fit(x_train)

        model = Sequential()
        model.add(Conv2D(32,3,padding="same", activation="relu", input_shape=(150,101,3)))
        model.add(MaxPool2D())

        model.add(Conv2D(32, 3, padding="same", activation="relu"))
        model.add(MaxPool2D())

        model.add(Conv2D(64, 3, padding="same", activation="relu"))
        model.add(MaxPool2D())
        model.add(Dropout(0.4))

        model.add(Flatten())
        model.add(Dense(128,activation="relu"))
        model.add(Dense(len(labels), activation="softmax"))

        opt = Adam(learning_rate=0.00001)
        model.compile(optimizer = opt , loss = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True) , metrics = ['accuracy'])

        model.fit(x_train,y_train, epochs = 2 , validation_data = (x_test, y_test))

        path = "/airflow/data/image_classifier/model"
        model.save(path)

   
    train_model_task =  PythonOperator(
        task_id='train_model',
        python_callable=train_model,
        trigger_rule="all_success",
    )

    prepare_data >> train_model_task