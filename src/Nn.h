#ifndef DEF_NN
#define DEF_NN

#include <iostream>

typedef struct
{

} node;

class NeuralNetwork
{
public:
    NeuralNetwork(uint input_nodes_num, uint hidden_nodes_num, uint output_nodes_num); //after : unint hidden layer num
    void train();

private:
    uint m_input_nodes_num;
    uint m_hidden_nodes_num;
    uint m_output_nodes_num;
};

#endif