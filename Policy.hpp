//
//  Policy.hpp
//  CS-776_Assignment_3
//
//  Created by Scott S Forer on 9/27/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#ifndef Policy_hpp
#define Policy_hpp

#include <stdio.h>

using namespace std;

class Policy
{
    friend class GA;
    friend class Parameters;
public:
    
    vector<int> binary;
    double prob_selection;
    double fitness;
    vector<double> x_vals;

};

#endif /* Policy_hpp */
