//
//  Parameters.hpp
//  CS-776_Assignment_3
//
//  Created by Scott S Forer on 9/27/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#ifndef Parameters_hpp
#define Parameters_hpp

#include <stdio.h>

using namespace std;

class Parameters
{
    friend class GA;
    friend class Policy;
public:
    //Parameters
    int pop_size = 100;
    int gen_max = 300;
    double xover_prob;
    double mutation_prob;
    double xover_prob_low = 0.2;
    double mutation_prob_low = 0.0001;
    double xover_prob_nom = 0.67;
    double mutation_prob_nom = 0.001;
    double xover_prob_high = 0.99;
    double mutation_prob_high = 0.01;
    int prob_case = 2;          //use 0 for low, 1 for nominal, 2 for high
    int function_num = 1;       //use 1 for DeJong 1, 2 for DeJong 2, 2 for DeJong 3, 4 for DeJong 4
    int num_x_vals;
    int bit_length;
    int string_size;
    int bits_to_flip = 1;
};

#endif /* Parameters_hpp */
