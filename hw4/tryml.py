from uwimg import *

def softmax_model(inputs, outputs):
    l = [make_layer(inputs, outputs, SOFTMAX)]
    return make_model(l)

def neural_net(inputs, outputs):
    print(inputs)
    l = [   make_layer(inputs, 64, LRELU),
            make_layer(64, 32, LRELU),
            make_layer(32, outputs, SOFTMAX)]
    return make_model(l)

print("loading data...")
train = load_classification_data(b"cifar.train", b"mnist.labels", 1)
test  = load_classification_data(b"cifar.test", b"mnist.labels", 1)
print("done")

print("training model...")
batch = 128
iters = 3000
rate = 0.1
momentum = .9
decay = 0.01

#m = softmax_model(train.X.cols, train.y.cols)
m = neural_net(train.X.cols, train.y.cols)
train_model(m, train, batch, iters, rate, momentum, decay)
print("done")

print("evaluating model...")
print("training accuracy: %f", accuracy_model(m, train))
print("test accuracy:     %f", accuracy_model(m, test))


## Questions ##

# 2.2.1 Why might we be interested in both training accuracy and testing accuracy? What do these two numbers tell us about our current model?
# Ans: This is because we want to test on how well the model does on our training data versus our testing data so that we can determine
# if the model overfits. If the training accuracy is high but test accuracy is not - it is an indication of overfitting. The model will not
# generalise well in production use.

# 2.2.2 Try varying the model parameter for learning rate to different powers of 10 (i.e. 10^1, 10^0, 10^-1, 10^-2, 10^-3) and training the model. What patterns do you see and how does the choice of learning rate affect both the loss during training and the final model accuracy?
# rate 0.01 - 0.9091 | 0.1 - 0.9171 | 1.0 - 0.8463
# The model loses accuracy quickly as the learning rate increases. Probably because it is overshooting the local minima in gradient descent.

# 2.2.3 Try varying the parameter for weight decay to different powers of 10: (10^0, 10^-1, 10^-2, 10^-3, 10^-4, 10^-5). How does weight decay affect the final model training and test accuracy?
# Increasing weight decay leads to a decreased test accuracy

# 2.3.1 Currently the model uses a logistic activation for the first layer. Try using a the different activation functions we programmed. How well do they perform? What's best?
# LOGISTIC = 0.9445 | RELU = 0.9439 | LRELU = 0.9512
# Leaky RELU performs best with a test accuracy of 95%

# 2.3.2 Using the same activation, find the best (power of 10) learning rate for your model. What is the training accuracy and testing accuracy?
# 0.1 is the best. The train accuracy is 0.956 and test is 0.9512

# 2.3.3 Right now the regularization parameter `decay` is set to 0. Try adding some decay to your model. What happens, does it help? Why or why not may this be?
# does not seem to help as it slightly decreases the test accuracy. However, it is actually preventing the model from overfitting.

# 2.3.4 Modify your model so it has 3 layers instead of two. The layers should be `inputs -> 64`, `64 -> 32`, and `32 -> outputs`. Also modify your model to train for 3000 iterations instead of 1000. Look at the training and testing error for different values of decay (powers of 10, 10^-4 -> 10^0). Which is best? Why?
# Performs better than previous model with 96% accuracy and with a decay of 0.01.

# 3.2.1 How well does your network perform on the CIFAR dataset?
# about the same

