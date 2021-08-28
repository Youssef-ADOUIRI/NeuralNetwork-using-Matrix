#ifndef DEF_NN
#define DEF_NN

#include <iostream>
#include <vector>
#include "matrix.h"

#define default_learningRate 0.5


class NeuralNetwork
{
public:
    NeuralNetwork(uint input_nodes_num, uint hidden_nodes_num, uint output_nodes_num); //after : unint hidden layer num
    ~NeuralNetwork();
    std::vector<double> forward_propagation(std::vector<double> input_data);
    void train(std::vector<double> listOfInputs , std::vector <double> targets);
    bool test(std::vector<double> Inputs , std::vector <double> targets);

private:
    uint m_input_nodes_num;
    uint m_hidden_nodes_num;
    uint m_output_nodes_num;
    float learning_rate;
    Matrix W_I_h;
    Matrix W_h_o;
};

void sigmoid(double &x);

#endif