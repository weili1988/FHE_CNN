import numpy as np
np.random.seed(123)  # for reproducibility
 
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution2D, MaxPooling2D, AveragePooling2D
from keras.utils import np_utils
from keras.datasets import mnist

# 4. Load pre-shuffled MNIST data into train and test sets
(X_train, y_train), (X_test, y_test) = mnist.load_data()
print('before one-hot: \n')
print('Xtrain Shape: ',X_train.shape)
print('Xtrain Shape: ',y_train.shape)
print('Xtrain Shape: ',X_test.shape)
print('Xtrain Shape: ',y_test.shape)
# 5. Preprocess input data
X_train = X_train.reshape(X_train.shape[0], 28, 28, 1)
X_test = X_test.reshape(X_test.shape[0], 28, 28, 1)
X_train = X_train.astype('float32')
X_test = X_test.astype('float32')
X_train /= 255
X_test /= 255

# 6. Preprocess class labels
Y_train = np_utils.to_categorical(y_train, 10)
Y_test = np_utils.to_categorical(y_test, 10)
print('after one-hot: \n')
print('Xtrain Shape: ',X_train.shape)
print('Ytrain Shape: ',Y_train.shape)
print('Xtest Shape: ',X_test.shape)
print('Ytest Shape: ',Y_test.shape)

# define model
model = Sequential()
model.add(Convolution2D(8, 3, 3, activation='relu', input_shape=(28,28,1)))
model.add(AveragePooling2D(pool_size=(2,2)))
model.add(Convolution2D(16, 3, 3, activation='relu'))
model.add(AveragePooling2D(pool_size=(4,4)))
model.add(Flatten())
model.add(Dense(10, activation='softmax'))

print('model summarry: \n')
print(model.summary())

# 8. Compile model
model.compile(loss='categorical_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])
# 9. Fit model on training data
model.fit(X_train, Y_train, 
          batch_size=64, nb_epoch=5, verbose=1)
# 10. Evaluate model on test data
score = model.evaluate(X_test, Y_test, verbose=0)
print("Loss: %.2f" % (score[0]))
print("Accuracy: %.2f%%" % (score[1]*100))
######################################################
print(model.layers[0].get_weights()[0].shape)
print(model.layers[0].get_weights()[1].shape)
print(model.layers[2].get_weights()[0].shape)
print(model.layers[2].get_weights()[1].shape)
print(model.layers[5].get_weights()[0].shape)
print(model.layers[5].get_weights()[1].shape)

CNN1_W = model.layers[0].get_weights()[0]
CNN1_b = model.layers[0].get_weights()[1]
CNN2_W = model.layers[2].get_weights()[0]
CNN2_b = model.layers[2].get_weights()[1]
FC_W = model.layers[5].get_weights()[0]
FC_b = model.layers[5].get_weights()[1]

CNN1_W = CNN1_W.reshape(CNN1_W.shape[0]*CNN1_W.shape[1]*CNN1_W.shape[2], CNN1_W.shape[3]) # 3*3 * nb_nodes
print(CNN1_W.shape)
print(CNN1_b.shape)
CNN2_W = CNN2_W.reshape(CNN2_W.shape[0]*CNN2_W.shape[1]*CNN2_W.shape[2], CNN2_W.shape[3]);
print(CNN2_W.shape)
print(CNN2_b.shape)
print(FC_W.shape)
print(FC_b.shape)

np.savetxt("./data/MNIST_CNN_weights/CNN1_W.txt", CNN1_W, fmt='%10.5f', delimiter=',')
np.savetxt("./data/MNIST_CNN_weights/CNN1_b.txt", CNN1_b, fmt='%10.5f', delimiter=',')
np.savetxt("./data/MNIST_CNN_weights/CNN2_W.txt", CNN2_W, fmt='%10.5f', delimiter=',')
np.savetxt("./data/MNIST_CNN_weights/CNN2_b.txt", CNN2_b, fmt='%10.5f', delimiter=',')
np.savetxt("./data/MNIST_CNN_weights/FC_W.txt", FC_W, fmt='%10.5f', delimiter=',')
np.savetxt("./data/MNIST_CNN_weights/FC_b.txt", FC_b, fmt='%10.5f', delimiter=',')
