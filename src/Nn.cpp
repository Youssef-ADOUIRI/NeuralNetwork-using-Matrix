
#include <iostream>
#include <math.h>
#include "Nn.h"

using namespace std;

NeuralNetwork::NeuralNetwork(uint input_nodes_num, uint hidden_nodes_num, uint output_nodes_num) : m_input_nodes_num(input_nodes_num), m_hidden_nodes_num(hidden_nodes_num), m_output_nodes_num(output_nodes_num), learning_rate(default_learningRate)
{

    W_I_h.Reset(m_hidden_nodes_num, m_input_nodes_num);
    W_h_o.Reset(m_output_nodes_num, m_hidden_nodes_num);

    W_I_h.nRand(0.0f, pow(m_hidden_nodes_num, -0.5));
    W_h_o.nRand(0.0f, pow(m_output_nodes_num, -0.5));
}

void sigmoid(double &x)
{
    x = 1 / (1 + exp(-x));
}

double max_vect = 0;
double sum = 0;

void softmax(double &x){
    x = exp(x-max_vect)/sum;
}

std::vector<double> NeuralNetwork::forward_propagation(vector<double> input_data)
{
    Matrix input(5, 4);
    input.fromVector(input_data);
    cout << "input : " << endl;
    input.print();

    Matrix X1 = W_I_h * input;
    cout << "X : " << endl;
    X1.print();

    Matrix H1 = X1;
    H1.apply(sigmoid);

    H1.print();

    /*second layer */

    Matrix X2 = W_h_o * H1;

    Matrix H2 = X2;

    H2.apply(sigmoid);

    Matrix output = H2;

    //softmax

    max_vect = output.get(0, 0);

    for (int i = 0; i < m_output_nodes_num; i++)
    {
        if (output.get(0, i) > max_vect)
        {
            max_vect = output.get(0, i);
        }
    }
    for (int i = 0; i < m_output_nodes_num; i++)
        sum += exp(output.get(0, i) - max_vect);

    output.apply(softmax);
    

    cout << "Final outputs : " << endl;
    output.print();

    vector<double> vecOutput = output.toVector();

    return vecOutput;
}

void NeuralNetwork:: train(vector<double> Inputs ,vector <double> targets){
    //prep
    // vector<double> f_result = forward_propagation(Inputs);
    //feed forward
    
    Matrix input(5, 4);
    input.fromVector(Inputs);
    Matrix X1 = W_I_h * input;
    Matrix H1 = X1;
    H1.apply(sigmoid);
    cout<<" \n\nokf\n\n";

    /*second layer */
    Matrix X2 = W_h_o * H1;
    

    Matrix H2 = X2;
    H2.apply(sigmoid);
    Matrix output = H2;
    //softmax
    max_vect = output.get(0, 0);
    for (int i = 0; i < m_output_nodes_num; i++)
    {
        if (output.get(0, i) > max_vect)
        {
            max_vect = output.get(0, i);
        }
    }
    for (int i = 0; i < m_output_nodes_num; i++)
        sum += exp(output.get(0, i) - max_vect);

    output.apply(softmax);

    cout<<" \n\nok\n\n";

    /* backprobagation: */
    Matrix m_target;
    m_target.fromVector(targets);

    
    Matrix m_result = output;
    m_result.multiply_by(-1.0);
    Matrix out_err = m_target + m_result;
    W_h_o.T();
    Matrix hidden_error = W_h_o * out_err;
    W_h_o.T();

    //gradient decent
    //first layer 
    Matrix term1 = Matrix::normal_multiplication(out_err , output);
    m_result.add_by(1.0); //1 - outputs_final
    Matrix term2 = Matrix::normal_multiplication(term1 , m_result);
    Matrix CopyH2 = H2; //for the others
    H2.T();
    Matrix delta_hidden_output = term2 * H2;
    delta_hidden_output.multiply_by(learning_rate);
    W_h_o = W_h_o + delta_hidden_output;

    //second layer
    term1 =  Matrix::normal_multiplication(hidden_error , CopyH2);
    CopyH2.multiply_by(-1.0);
    CopyH2.add_by(1.0); // 1 - hidden_output
    term2 = Matrix::normal_multiplication(term1 , CopyH2);
    Matrix inputs_cpy = input;
    inputs_cpy.T();
    delta_hidden_output = term2 * inputs_cpy;
    delta_hidden_output.multiply_by(learning_rate);
    W_I_h = W_I_h + delta_hidden_output;
    
}


bool  NeuralNetwork::test(vector<double> Inputs ,vector <double> targets){

    vector<double> results = forward_propagation(Inputs);
    if(results.size() != targets.size()){
        cerr<<"error in sizes"<<endl;
        return false;
    }
    uint max1 = 0; 
    uint max2 = 0;

    for(int i = 0 ; i < results.size() ; i++ ){
        if(results[i]> results[max1]){
            max1 = i;
        }
        if(targets[i]> targets[max2]){
            max2 = i;
        }
    }

    if(max1 == max2){
        return true;
    }
    else{
        return false;
    }

}



NeuralNetwork::~NeuralNetwork(){}