#include <iostream>
#include "Nn.h"
#include <math.h>
#include <vector>
#include <fstream>
#include <list>
#include <sstream>
//#include "CImg.h"

using namespace std;

int main()
{
  uint intp_num = 784;
  uint out_num = 10;
  NeuralNetwork nn(intp_num, 100, out_num);

  /*prepare data*/
  string const File_Name = "/home/youssef/Documents/DataSet/MMsit/mnist_test_10.csv";
  ifstream Data_file;
  string line;
  list<string> lines_list;
  Data_file.open(File_Name.c_str());
  cout<<"start"<<endl;

  if (Data_file.is_open())
  {
    while (getline(Data_file, line))
    { //get all lines in a list
      lines_list.push_back(line);
    }

    for (string lin : lines_list)
    {
      //extract a line to numbers
      string num_string;
      const char dilm = ',';
      stringstream string_stream(lin.c_str());
      vector<double> elements;
      int temp, count = 0;
      bool isTarget = true;
      double target_index;
      while (getline(string_stream, num_string, dilm))
      {
        if (stringstream(num_string) >> temp)
        {
          if (isTarget)
          {
            target_index = (double)temp;
            isTarget = false;
          }
          else
          {
            elements.push_back((double)temp);
            count++;
          }
        }
      }

     

      cout<<"target : "<<target_index<<endl;


      vector<double> targets(10 , 0.01);
      targets[target_index] = 1.0;

       if(count != 784 || targets.size() != out_num || elements.size() != intp_num ){
        cerr<<"error"<<endl;
        exit(0);
      }

      nn.train(elements , targets);

    }

    cout<<"finish"<<endl;



    //nn.test();





    /*
    //show the result
    cimg_library::CImg<unsigned char> theImage(28, 28, 1, 3, 1);

    for (int i = 0; i < 28; i++)
    {
      for (int j = 0; j < 28; j++)
      {
        theImage(j, i) = element[(i * 28 + j) + 1];
      }
    }

    theImage.save_bmp("output.bmp" ); // write it
    */
    Data_file.close();
  }
  else
  {
    cerr << "Unable to open the file!" << endl;
  }

  return 0;
}
