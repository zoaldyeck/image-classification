# image-classification
classify 

An image classification system to classify an image as a poster(with text in it) or an original(without text).

The whole process goes in 2 steps:
1. Training phase: Train a SVM model on a training data set of 30k labeled posters and originals.
2. Testing phase: Test new image with the SVM model and k centroids, returning either +1(poster) or -1(original).
