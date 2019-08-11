#include <math.h>
#include <stdlib.h>
#include "image.h"
#include "matrix.h"

// Run an activation function on each element in a matrix,
// modifies the matrix in place
// matrix m: Input to activation function
// ACTIVATION a: function to run
void activate_matrix(matrix m, ACTIVATION a)
{
    int i, j;
    for(i = 0; i < m.rows; ++i)
    {
        double sum = 0;

        for(j = 0; j < m.cols; ++j)
        {
            double x = m.data[i][j];

            if(a == LOGISTIC) {
                // TODO
                m.data[i][j] = 1.0 / (1.0 + exp(-x));
            } 
            else if (a == RELU){
                // TODO
                m.data[i][j] = MAX(0, x);
            } 
            else if (a == LRELU){
                // TODO
                m.data[i][j] = (x > 0 ? x : 0.01*x);
            } 
            else if (a == SOFTMAX){
                // TODO
                m.data[i][j] = exp(x);
            }
            sum += m.data[i][j];
        }
        if (a == SOFTMAX) 
        {
            // TODO: have to normalize by sum if we are using SOFTMAX
            for (j = 0; j < m.cols; ++j){
                m.data[i][j] /= sum;
            }
        }
    }
}

// Calculates the gradient of an activation function and multiplies it into
// the delta for a layer
// matrix m: an activated layer output
// ACTIVATION a: activation function for a layer
// matrix d: delta before activation gradient
void gradient_matrix(matrix m, ACTIVATION a, matrix d)
{
    int i, j;
    for(i = 0; i < m.rows; ++i)
    {
        for(j = 0; j < m.cols; ++j)
        {
            double x = m.data[i][j];

            // TODO: multiply the correct element of d by the gradient
            if (a == LOGISTIC) {
                d.data[i][j] *= x * (1 - x);
            }
            else if (a == RELU) {
                d.data[i][j] *= x > 0 ? 1 : 0;
            }
            else if (a == LRELU) {
                d.data[i][j] *= x > 0 ? 1 : 0.1;
            }
        }
    }
}

// Forward propagate information through a layer
// layer *l: pointer to the layer
// matrix in: input to layer
// returns: matrix that is output of the layer
matrix forward_layer(layer *l, matrix in)
{
    l->in = in;  // Save the input for backpropagation

    // TODO: fix this! multiply input by weights and apply activation function.
    matrix out = matrix_mult_matrix(in, l->w);
    activate_matrix(out, l->activation);

    free_matrix(l->out);// free the old output
    l->out = out;       // Save the current output for gradient calculation
    return out;
}

// Backward propagate derivatives through a layer
// layer *l: pointer to the layer
// matrix delta: partial derivative of loss w.r.t. output of layer
// returns: matrix, partial derivative of loss w.r.t. input to layer
matrix backward_layer(layer *l, matrix delta)
{
    // 1.4.1
    // delta is dL/dy
    // TODO: modify it in place to be dL/d(xw)
    gradient_matrix(l->out, l->activation, delta);


    // 1.4.2
    // TODO: then calculate dL/dw and save it in l->dw
    free_matrix(l->dw);
    matrix in_transposed = transpose_matrix(l->in);
    matrix dw = matrix_mult_matrix(in_transposed, delta);
    l->dw = dw;
    free_matrix(in_transposed);

    
    // 1.4.3
    // TODO: finally, calculate dL/dx and return it.
    matrix wt = transpose_matrix(l->w);
    matrix dx = matrix_mult_matrix(delta, wt);
    free_matrix(wt);

    return dx;
}

// Update the weights at layer l
// layer *l: pointer to the layer
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_layer(layer *l, double rate, double momentum, double decay)
{
    // TODO:
    // Calculate Δw_t = dL/dw_t - λw_t + mΔw_{t-1}
    // save it to l->v
    matrix axpy_wt = axpy_matrix(-decay, l->w, l->dw);
    matrix delta_wt = axpy_matrix(momentum, l->v, axpy_wt);
    free_matrix(l->v);
    l->v = delta_wt;

    // Update l->w
    matrix wt = axpy_matrix(rate, delta_wt, l->w);
    free_matrix(l->w);
    l->w = wt;

    // Remember to free any intermediate results to avoid memory leaks
    free_matrix(axpy_wt);
}

// Make a new layer for our model
// int input: number of inputs to the layer
// int output: number of outputs from the layer
// ACTIVATION activation: the activation function to use
layer make_layer(int input, int output, ACTIVATION activation)
{
    layer l;
    l.in  = make_matrix(1,1);
    l.out = make_matrix(1,1);
    l.w   = random_matrix(input, output, sqrt(2./input));
    l.v   = make_matrix(input, output);
    l.dw  = make_matrix(input, output);
    l.activation = activation;
    return l;
}

// Run a model on input X
// model m: model to run
// matrix X: input to model
// returns: result matrix
matrix forward_model(model m, matrix X)
{
    int i;
    for(i = 0; i < m.n; ++i){
        X = forward_layer(m.layers + i, X);
    }
    return X;
}

// Run a model backward given gradient dL
// model m: model to run
// matrix dL: partial derivative of loss w.r.t. model output dL/dy
void backward_model(model m, matrix dL)
{
    matrix d = copy_matrix(dL);
    int i;
    for(i = m.n-1; i >= 0; --i){
        matrix prev = backward_layer(m.layers + i, d);
        free_matrix(d);
        d = prev;
    }
    free_matrix(d);
}

// Update the model weights
// model m: model to update
// double rate: learning rate
// double momentum: amount of momentum to use
// double decay: value for weight decay
void update_model(model m, double rate, double momentum, double decay)
{
    int i;
    for(i = 0; i < m.n; ++i){
        update_layer(m.layers + i, rate, momentum, decay);
    }
}

// Find the index of the maximum element in an array
// double *a: array
// int n: size of a, |a|
// returns: index of maximum element
int max_index(double *a, int n)
{
    if(n <= 0) return -1;
    int i;
    int max_i = 0;
    double max = a[0];
    for (i = 1; i < n; ++i) {
        if (a[i] > max){
            max = a[i];
            max_i = i;
        }
    }
    return max_i;
}

// Calculate the accuracy of a model on some data d
// model m: model to run
// data d: data to run on
// returns: accuracy, number correct / total
double accuracy_model(model m, data d)
{
    matrix p = forward_model(m, d.X);
    int i;
    int correct = 0;
    for(i = 0; i < d.y.rows; ++i){
        if(max_index(d.y.data[i], d.y.cols) == max_index(p.data[i], p.cols)) ++correct;
    }
    return (double)correct / d.y.rows;
}

// Calculate the cross-entropy loss for a set of predictions
// matrix y: the correct values
// matrix p: the predictions
// returns: average cross-entropy loss over data points, 1/n Σ(-ylog(p))
double cross_entropy_loss(matrix y, matrix p)
{
    int i, j;
    double sum = 0;
    for(i = 0; i < y.rows; ++i){
        for(j = 0; j < y.cols; ++j){
            sum += -y.data[i][j]*log(p.data[i][j]);
        }
    }
    return sum/y.rows;
}


// Train a model on a dataset using SGD
// model m: model to train
// data d: dataset to train on
// int batch: batch size for SGD
// int iters: number of iterations of SGD to run (i.e. how many batches)
// double rate: learning rate
// double momentum: momentum
// double decay: weight decay
void train_model(model m, data d, int batch, int iters, double rate, double momentum, double decay)
{
    int e;
    for(e = 0; e < iters; ++e){
        data b = random_batch(d, batch);
        matrix p = forward_model(m, b.X);
        fprintf(stderr, "%06d: Loss: %f\n", e, cross_entropy_loss(b.y, p));
        matrix dL = axpy_matrix(-1, p, b.y); // partial derivative of loss dL/dy
        backward_model(m, dL);
        update_model(m, rate/batch, momentum, decay);
        free_matrix(dL);
        free_data(b);
    }
}
