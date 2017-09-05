Deep CNN Implemented on Fully Homomophic Encrypted Data
====================

Overview
--------
A deep convolutional neural network (CNN) has been implement with [Simple Encrypted Arithmetic Library - SEAL](https://sealcrypto.codeplex.com/). A classical image classfication example has been used to demonstrates the capability of CNN-based inference on fully homomophic encrpted (FHE) handwritted digts from [the MNIST database](http://yann.lecun.com/exdb/mnist/). The regular CNN model based on [Keras](https://keras.io/) has been developed to generate the weights. The Python implementation is also included.

Requirements
------------
* [Simple Encrypted Arithmetic Library - SEAL](https://sealcrypto.codeplex.com/)
* Visual Studio Community 2017
* Python 2.6 or up or Python 3

Installation
-----
* Install [SEAL](https://sealcrypto.codeplex.com/)
* INSIDE the SEAL folder, clone this respo. This step will create a project folder named "FHE\_CNN" insdie SEAL 
* Open the SEAL project with VS 2017, and add the folder "FHE\_CNN" as an exsiting project
* Set the project "FHE\_CNN" as the startup project
