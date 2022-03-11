import cv2
import sys
import os
import tensorflow as tf
import numpy as np




def main():
    path = sys.argv[1]
    if not os.path.exists(path):
        print("File not found")
        return
    
    labels = ['Action', 'Adventure', 'Animation', 'Comedy', 'Crime', 'Documentary', 'Drama', 'Family', 'Fantasy', 'Foreign', 
         'History', 'Horror', 'Music', 'Mystery', 'Romance', 'Science Fiction', 'Thriller', 'TV Movie', 'War', 'Western']

    img_arr = cv2.imread(path)[...,::-1] #convert BGR to RGB format
    resized_arr = cv2.resize(img_arr, (101, 150)) # Reshaping images to preferred size
    normalized_arr = resized_arr / 255
    normalized_arr = np.array([normalized_arr])
    normalized_arr.reshape(-1, 150, 101, 1)

    model = tf.keras.models.load_model("image_classifier/model")
    predictions = model.predict(normalized_arr)
    image_name  = path.split("/")[-1]
    print(f'\n{image_name} most likely to be a {labels[np.argmax(predictions[0])]} movie \n')
    print('all probabilities: ')
    for i in range(len(labels)):
        formatted_output="{:.4f}".format(predictions[0][i])
        print(f'{labels[i]} : {formatted_output}')

if __name__ == "__main__":
    main()
   
    