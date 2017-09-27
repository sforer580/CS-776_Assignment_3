//
//  GA.hpp
//  CS-776_Assignment_3
//
//  Created by Scott S Forer on 9/27/17.
//  Copyright © 2017 Scott S Forer. All rights reserved.
//

#ifndef GA_hpp
#define GA_hpp

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

class GA
{
    friend class Policy;
    friend class Parameters;
public:
    Parameters* pP;
    vector<Policy> parent;
    vector<Policy> child;
    
    vector<double> max_fitness;
    vector<double> ave_fitness;
    vector<double> min_fitness;
    
    void Select_Case();
    void Build_Pop();
    void Decode();
    void DeJong_1();
    void DeJong_2();
    void DeJong_3();
    void DeJong_4();
    void Run_DeJong_Function();
    void Get_Probability_Selection();
    void Crossover();
    void Mutation();
    struct Greater_than_par_fitness;
    void Sort_Policies_By_Fitness();
    void Save_Data();
    void Write_Data_To_File();
    void Delete_Text_Files();
    
    void Run_GA(int sr);
};

//-----------------------------------------------------------
//selects the experiment
void GA::Select_Case()
{
    if (pP->prob_case == 0)
    {
        pP->bit_length = 30;
        pP->xover_prob = pP->xover_prob_low;
        pP->mutation_prob = pP->mutation_prob_low;
    }
    if (pP->prob_case == 1)
    {
        pP->xover_prob = pP->xover_prob_nom;
        pP->mutation_prob = pP->mutation_prob_nom;
    }
    if (pP->prob_case == 2)
    {
        pP->xover_prob = pP->xover_prob_high;
        pP->mutation_prob = pP->mutation_prob_high;
    }
    if (pP->function_num==1)
    {
        pP->bit_length = 30;
        pP->num_x_vals = 3;
        pP->string_size = pP->bit_length/pP->num_x_vals;
    }
    if (pP->function_num==2)
    {
        pP->bit_length = 24;
        pP->num_x_vals = 2;
        pP->string_size = pP->bit_length/pP->num_x_vals;
    }
    if (pP->function_num==3)
    {
        pP->bit_length = 50;
        pP->num_x_vals = 5;
        pP->string_size = pP->bit_length/pP->num_x_vals;
    }
    if (pP->function_num==4)
    {
        pP->bit_length = 30;
        pP->num_x_vals = 3;
        pP->string_size = pP->bit_length/pP->num_x_vals;
    }
}


//-----------------------------------------------------------
//builds pop_size number of solutions
void GA::Build_Pop()
{
    for (int p=0; p<pP->pop_size; p++)
    {
        Policy Pa;
        parent.push_back(Pa);
        Policy C;
        child.push_back(C);
        
        //randomly assigns an bit value (0 or 1)
        for (int num=0; num<pP->bit_length; num++)
        {
            parent.at(p).binary.push_back(rand() % 2);
            child.at(p).binary.push_back(0);
        }
        for (int x=0; x<pP->num_x_vals; x++)
        {
            parent.at(p).x_vals.push_back(0);
            child.at(p).x_vals.push_back(0);
        }
    }
}

//-----------------------------------------------------------
//Decodes the binary vector into x values
void GA::Decode()
{
    for (int p=0; p<pP->pop_size; p++)
    {
        for(int x=0; x<pP->num_x_vals; x++)
        {
            parent.at(p).x_vals.at(x) = 0;
            double num = 0;
            for(int i=0; i<pP->string_size; i++)
            {
                num += parent.at(p).binary.at(i+x*pP->string_size)*pow(2,i);
            }
            parent.at(p).x_vals.at(x) = num;
        }
    }
}




//-----------------------------------------------------------
//DeJong function 1
//sum of x_i*x_i
void GA::DeJong_1()
{
    for (int p=0; p<pP->pop_size; p++)
    {
        for (int x=0; x<pP->num_x_vals; x++)
        {
            parent.at(p).x_vals.at(x) = (512-parent.at(p).x_vals.at(x))/100;
            assert(parent.at(p).x_vals.at(x)<=5.12 && parent.at(p).x_vals.at(x)>=-5.12);
        }
    }

    
    for(int p=0; p<pP->pop_size; p++)
    {
        parent.at(p).fitness = 0;
        double sum = 0;
        for (int x=0; x<pP->num_x_vals; x++)
        {
            //cout << "cp" << endl;
            //cout << par.at(pop).val.at(v) << endl;
            sum += (parent.at(p).x_vals.at(x)*parent.at(p).x_vals.at(x));
        }
        parent.at(p).fitness = sum;
    }
}


//-----------------------------------------------------------
//DeJong function 2
//100*((x1*x1-x2)*(x1*x1-x2))+((1-x1)*(1-x1))
void GA::DeJong_2()
{
    for (int p=0; p<pP->pop_size; p++)
    {
        for (int x=0; x<pP->num_x_vals; x++)
        {
            parent.at(p).x_vals.at(x) = (2048-parent.at(p).x_vals.at(x))/1000;
            assert(parent.at(p).x_vals.at(x)<=2.048 && parent.at(p).x_vals.at(x)>=-2.048);
        }
    }
    
    for(int p=0; p<pP->pop_size; p++)
    {
        parent.at(p).fitness = 0;
        double sum = 0;
        double x1 = 0;
        double x2 = 0;
        x1 = parent.at(p).x_vals.at(0);
        x2 = parent.at(p).x_vals.at(1);
        sum = 100*(((x1*x1)-x2)*((x1*x1)-x2))+((1-x1)*(1-x1));
        parent.at(p).fitness = sum;
    }
}

//-----------------------------------------------------------
//DeJong function 3
//sum of
void GA::DeJong_3()
{
    
}


//-----------------------------------------------------------
//DeJong function 4
//sum of
void GA::DeJong_4()
{
    
}



//-----------------------------------------------------------
//Runs DeJong fucntion
void GA::Run_DeJong_Function()
{
    if (pP->function_num==1)
    {
        DeJong_1();
    }
    if (pP->function_num==2)
    {
        DeJong_2();
    }
    if (pP->function_num==3)
    {
        DeJong_3();
    }
    if (pP->function_num==4)
    {
        DeJong_4();
    }
}



//-----------------------------------------------------------
//calculates the prob of selection for each parent
void GA::Get_Probability_Selection()
{
    double total_fitness = 0;
    for (int p=0; p<pP->pop_size; p++)
    {
        total_fitness += parent.at(p).fitness;
    }

    double check_prob;
    for (int p=0; p<pP->pop_size; p++)
    {
        parent.at(p).prob_selection = 0;
        parent.at(p).prob_selection = parent.at(p).fitness/total_fitness;
        check_prob += parent.at(p).prob_selection;
    }
    //cout << check_prob << endl;
    //assert(check_prob == 1);
}


//-----------------------------------------------------------
//selects two parents based on their prob of selection then runs xover
void GA::Crossover()
{
    //sort(par.begin(), par.end(), less_than_prob_selection());
    for (int p=0; p<pP->pop_size/2; p++)
    {
        //selects parent based on its probaility of selection
        int index_1 = 0;
        double stat_1 = 0;
        double r = ((double) rand() / (RAND_MAX));
        for (int i=0; i<pP->pop_size; i++)
        {
            stat_1 = stat_1 + parent.at(i).prob_selection;
            if(stat_1 > r)
            {
                index_1 = i;
                break;
            }
        }
        //cout << "index 1" << "\t" << "is" << "\t" << index_1 << endl;
        int index_2 = 0;
        double stat_2 = 0;
        double rr = ((double) rand() / (RAND_MAX));
        for (int i=0; i<pP->pop_size; i++)
        {
            stat_2 = stat_2 + parent.at(i).prob_selection;
            if(stat_2 > rr)
            {
                index_2 = i;
                break;
            }
        }
        
        child.at(2*p) = parent.at(index_1);
        child.at(2*p+1) = parent.at(index_2);
        
        
        //gene cross over
        double rrr = ((double) rand() / (RAND_MAX));
        if (rrr <= pP->xover_prob)
        {
            int r_splice = (int)rand() % pP->bit_length;
            //cout << "splice location" << "\t" << r_splice << endl;
            vector<int> place_holder_1;
            vector<int> place_holder_2;
            for (int i=0; i<pP->bit_length-r_splice; i++)
            {
                place_holder_1.push_back(child.at(2*p).binary.at(r_splice+i));
                place_holder_2.push_back(child.at(2*p+1).binary.at(r_splice+i));
            }
            //cout << "cp" << endl;
            for (int i=0; i<pP->bit_length-r_splice; i++)
            {
                child.at(p).binary.at(r_splice+i) = place_holder_2.at(i);
                child.at(p+1).binary.at(r_splice+i) = place_holder_1.at(i);
            }
            //cout << "cp" << endl;
        }
    }
}


//-----------------------------------------------------------
//mutates the copies if the probability of mutation statement is met
void GA::Mutation()
{
    for (int p=0; p<pP->pop_size; p++)
    {
        double r_mutate = ((double) rand() / (RAND_MAX));
        if (r_mutate < pP->mutation_prob)
        {
            //cout << "mutation" << endl;
            Policy M;
            M = child.at(p);
            //Mutates the child if the probability of mutation statement is met
            
            for (int i=0; i<pP->bits_to_flip; i++)
            {
                int bitdex = (int)rand() % pP->bit_length;
                if (M.binary.at(bitdex) == 0)
                {
                    M.binary.at(bitdex) = 1;
                }
                else
                {
                    M.binary.at(bitdex) = 0;
                }
            }
            child.at(p) = M;
        }
    }
}


//-----------------------------------------------------------
//sorts the population based on their fitness from highest to lowest
struct GA::Greater_than_par_fitness
{
    inline bool operator() (const Policy& struct1, const Policy& struct2)
    {
        return (struct1.fitness > struct2.fitness);
    }
};


//-------------------------------------------------------------------------
//Sorts population
void GA::Sort_Policies_By_Fitness()
{
    for (int p=0; p<pP->pop_size; p++)
    {
        sort(parent.begin(), parent.end(), Greater_than_par_fitness());
    }
}


//-----------------------------------------------------------
//Saves statistical data
void GA::Save_Data()
{
    max_fitness.push_back(parent.at(0).fitness);
    
    double sum = 0;
    for (int p=0; p<pP->pop_size; p++)
    {
        sum+= parent.at(p).fitness;
    }
    ave_fitness.push_back(sum/double(pP->pop_size));
    
    min_fitness.push_back(parent.back().fitness);
    
}


//-----------------------------------------------------------
//Writes the statistical data to a txt file
void GA::Write_Data_To_File()
{
    //
    ofstream File1;
    File1.open("Max_Fitness.txt", ios_base::app);
    for (int gen=0; gen<pP->gen_max; gen++)
    {
        File1 << max_fitness.at(gen) << "\t";
    }
    File1 << endl;
    File1.close();
    
    
    //
    ofstream File2;
    File2.open("Ave_Fitness.txt", ios_base::app);
    for (int gen=0; gen<pP->gen_max; gen++)
    {
        File2 << ave_fitness.at(gen) << "\t";
    }
    File2 << endl;
    File2.close();
    
    
    //
    ofstream File3;
    File3.open("Min_Fitness.txt", ios_base::app);
    for (int gen=0; gen<pP->gen_max; gen++)
    {
        File3 << min_fitness.at(gen) << "\t";
    }
    File3 << endl;
    File3.close();

    
    //
    ofstream File4;
    File4.open("Parameters.txt");
    File4 << "GA PARAMETERS" << endl;
    File4 << "Population Size" << "\t" << "\t" << pP->pop_size << endl;
    File4 << "Number of Generations" << "\t" << pP->gen_max << endl;
    File4 << "Crossover Probability" << "\t" << pP->xover_prob << endl;
    File4 << "Mutation Probability" << "\t" << pP->mutation_prob << endl;
    File4 << "Bits to Flip" << "\t" << pP->bits_to_flip << endl;
    File4 << "DeJong Function" << "\t" << "\t" << pP->function_num << endl;
    File4 << "Number of x values" << "\t" << pP->num_x_vals << endl;
    File4 << "Bit Length" << "\t" << pP->bit_length << endl;
    File4 << "String Size" << "\t" << pP->string_size << endl;
    File4.close();
}


//-----------------------------------------------------------
//Deletes the txt files
void GA::Delete_Text_Files()
{
    //
    if( remove( "Max_Fitness.txt" ) != 0 )
        perror( "ERROR DELETING FILE Max_Fitness" );
    else
        puts( "Max_Fitness FILE SUCCESSFULLY DELETED" );
    cout << endl;
    
    
    //
    if( remove( "Ave_Fitness.txt" ) != 0 )
        perror( "ERROR DELETING FILE Ave_Fitness" );
    else
        puts( "Ave_Fitness FILE SUCCESSFULLY DELETED" );
    cout << endl;
    
    
    //
    if( remove( "Min_Fitness.txt" ) != 0 )
        perror( "ERROR DELETING FILE Min_Fitness" );
    else
        puts( "Min_Fitness FILE SUCCESSFULLY DELETED" );
    cout << endl;
}



//-----------------------------------------------------------
//Runs entire GA
void GA::Run_GA(int sr)
{
    if (sr==0)
    {
        Delete_Text_Files();
    }
    Select_Case();
    Build_Pop();
    for (int gen=0; gen<pP->gen_max; gen++)
    {
        if (gen%10==0)
        {
            //cout << "SR::" << sr << "\t" << "GEN::" << gen << endl;
        }
        if (gen<=pP->gen_max-1)
        {
            Decode();
            Run_DeJong_Function();
            Get_Probability_Selection();
            Crossover();
            Mutation();
            Sort_Policies_By_Fitness();
            //cout << "BEST FITNESS" << "\t" << parent.at(0).fitness << endl;
            Save_Data();
            parent = child;
        }
        if (gen==pP->gen_max)
        {
            Decode();
            Run_DeJong_Function();
            Sort_Policies_By_Fitness();
            //cout << "BEST FITNESS" << "\t" << parent.at(0).fitness << endl;
            Save_Data();
        }
    }
    Write_Data_To_File();
}

#endif /* GA_hpp */
