//
//  main.cpp
//  CS-776_Assignment_3
//
//  Created by Scott S Forer on 9/24/16.
//  Copyright © 2016 Scott S Forer. All rights reserved.
//


#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <ctime>
#include <random>



using namespace std;

#include "Parameters.hpp"
#include "Policy.hpp"
#include "GA.hpp"






//-----------------------------------------------------------
//-----------------------------------------------------------
int main()
{
    srand(time(NULL));
    for (int sr=0; sr<30; sr++)
    {
        cout << "SR" << sr << endl;
        Parameters P;
        GA G;
        G.pP = &P;
        G.Run_GA(sr);
    }
    return 0;
}
