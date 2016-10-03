//
//  main.cpp
//  CS-776_Assignment_3
//
//  Created by Scott S Forer on 9/24/16.
//  Copyright © 2016 Scott S Forer. All rights reserved.
//


#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <fstream>

#define bit_length 50           //use 30 for case 0, 24 for case 1, 50 for case 2

using namespace std;


class Parent
{
public:
    int combo[bit_length];
    double prob_selection;
    double fitness;
    vector<double> val;
    double val_1;
    double val_2;
    double val_3;
    double val_4;
    double val_5;
};

class Child
{
public:
    int combo[bit_length];
};


vector<Parent> par;
vector<Child> chi;


//Parameters
int pop_size = 1000;
int gen_max = 400;
double xover_prob;
double mutation_prob;
double xover_prob_low = 0.2;
double mutation_prob_low = 0.0001;
double xover_prob_nom = 0.67;
double mutation_prob_nom = 0.001;
double xover_prob_high = 0.99;
double mutation_prob_high = 0.01;
int prob_case = 2;          //use 0 for low, 1 for nominal, 2 for high
int function_num = 2;       //use 0 for DeJong 1, 1 for DeJong 2, 2 for DeJong 3, 3 for DeJong 4
int num_val;


//Statistics
int num_trials = 30;
int array_split;
vector<double> max_fitness;
vector<double> ave_fitness;
vector<double> min_fitness;
double ave_max_fitness;
double ave_ave_fitness;
double ave_min_fitness;
vector<double> ave_of_ave_max;
vector<double> ave_of_ave_ave;
vector<double> ave_of_ave_min;


//For testing purposes
int target[bit_length];




//-----------------------------------------------------------
//selects the probability of xover and mutation
void select_prob_case()
{
    if (prob_case == 0)
    {
        xover_prob = xover_prob_low;
        mutation_prob = mutation_prob_low;
        array_split = 10;
        num_val = 3;
    }
    if (prob_case == 1)
    {
        xover_prob = xover_prob_nom;
        mutation_prob = mutation_prob_nom;
        array_split = 12;
        num_val = 2;
    }
    if (prob_case == 2)
    {
        xover_prob = xover_prob_high;
        mutation_prob = mutation_prob_high;
        array_split = 10;
        num_val = 5;
    }
}



//-----------------------------------------------------------
//builds pop_size number of solutions
void build_pop()
{
    for (int pop=0; pop<pop_size; pop++)
    {
        Parent P;
        par.push_back(P);
        Child C;
        chi.push_back(C);
        
        //randomly assigns an bit value (0 or 1)
        for (int num=0; num<bit_length; num++)
        {
            par.at(pop).combo[num] = rand() % 2;
            //par.at(pop).combo[num] = 0;
            //if (num == 9)
            //{
                //par.at(pop).combo[num] = 1;
            //}
            //if (num == 19)
            //{
                //par.at(pop).combo[num] = 1;
            //}
            //if (num == 29)
            //{
                //par.at(pop).combo[num] = 1;
            //}
            chi.at(pop).combo[num] = 0;
            //cout << "cp" << endl;
            //cout << par.at(pop).combo[num] << "\t";
        }
        //cout << endl;
        par.at(pop).val.resize(num_val);
        for (int v=0; v<num_val; v++)
        {
            par.at(pop).val.at(v) = 0;
        }
        //cout << endl;
    }
}


//-----------------------------------------------------------
//builds target array --> testing purposes only
void create_target()
{
    for (int num=0; num<bit_length; num++)
    {
        //target[num] = rand() % 2;
        target[num] = 1;
    }
    //cout << "target" << endl;
    //for (int num=0; num<bit_length; num++)
    //{
    //cout << target[num] << "\t";
    //}
    //cout << endl;
    //cout << endl;
}



//-----------------------------------------------------------
//compares the parent bit by bit to the target --> testing purposes only
void easy_eval()
{
    for(int pop=0; pop<pop_size; pop++)
    {
        par.at(pop).fitness = 0;
        int fit = 0;
        for (int num=0; num<bit_length; num++)
        {
            if(par.at(pop).combo[num] == target[num])
            {
                fit ++;
            }
            else
            {
                continue;
            }
        }
        par.at(pop).fitness = fit;
    }
}


//-----------------------------------------------------------
//DeJong function 1
//sum of x_i*x_i
void DeJong_1()
{
    for (int pop=0; pop<pop_size; pop++)
    {
        par.at(pop).val_1 = 0;
        par.at(pop).val_2 = 0;
        par.at(pop).val_3 = 0;
        for (int num=0; num<10; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_1 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_1 += pow(2,(num));
            }
            //cout << par.at(pop).val_1 << endl;
        }
        for (int num=10; num<20; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_2 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_2 += pow(2,(num-10));
            }
        }
        for (int num=20; num<30; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_3 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_3 += pow(2,(num-20));
            }
        }
        //cout << par.at(pop).val_1 << endl;
        //cout << par.at(pop).val_2 << endl;
        //cout << par.at(pop).val_3 << endl;
        par.at(pop).val.clear();
        par.at(pop).val.push_back((512-par.at(pop).val_1)/100);
        par.at(pop).val.push_back((512-par.at(pop).val_2)/100);
        par.at(pop).val.push_back((512-par.at(pop).val_3)/100);
        
        //cout << "x1" << "\t" << par.at(pop).val.at(0) << endl;
        //cout << "x2" << "\t" << par.at(pop).val.at(1) << endl;
        //cout << "x3" << "\t" << par.at(pop).val.at(2) << endl;
    }
    
    for(int pop=0; pop<pop_size; pop++)
    {
        par.at(pop).fitness = 0;
        double sum = 0;
        for (int v=0; v<3; v++)
        {
            //cout << "cp" << endl;
            //cout << par.at(pop).val.at(v) << endl;
            sum += (par.at(pop).val.at(v)*par.at(pop).val.at(v));
        }
        //cout << endl;
        //cout << "sum" << "\t" << sum  << endl;
        par.at(pop).fitness = 1.0/(sum+1);
        //cout << "fitness" << "\t" << par.at(pop).fitness << endl;
    }
}


//-----------------------------------------------------------
//DeJong function 2
//100*((x1*x1-x2)*(x1*x1-x2))+((1-x1)*(1-x1))
void DeJong_2()
{
    for (int pop=0; pop<pop_size; pop++)
    {
        par.at(pop).val_1 = 0;
        par.at(pop).val_2 = 0;
        for (int num=0; num<12; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_1 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_1 += pow(2,(num));
            }
            //cout << par.at(pop).val_1 << endl;
        }
        for (int num=12; num<24; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_2 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_2 += pow(2,(num-12));
            }
        }
        //cout << par.at(pop).val_1 << endl;
        //cout << par.at(pop).val_2 << endl;
        par.at(pop).val.clear();
        par.at(pop).val.push_back((2048-par.at(pop).val_1)/1000);
        par.at(pop).val.push_back((2048-par.at(pop).val_2)/1000);
        
        //cout << "x1" << "\t" << par.at(pop).val.at(0) << endl;
        //cout << "x2" << "\t" << par.at(pop).val.at(1) << endl;
    }
    
    for(int pop=0; pop<pop_size; pop++)
    {
        par.at(pop).fitness = 0;
        double sum = 0;
        double x1 = 0;
        double x2 = 0;
        x1 = par.at(pop).val.at(0);
        x2 = par.at(pop).val.at(1);
        sum = 100*(((x1*x1)-x2)*((x1*x1)-x2))+((1-x1)*(1-x1));
        par.at(pop).fitness = 1/(sum+1);
    }
}


//-----------------------------------------------------------
//DeJong function 3
//sum of interger*(x_i)
void DeJong_3()
{
    for (int pop=0; pop<pop_size; pop++)
    {
        par.at(pop).val_1 = 0;
        par.at(pop).val_2 = 0;
        par.at(pop).val_3 = 0;
        par.at(pop).val_4 = 0;
        par.at(pop).val_5 = 0;
        for (int num=0; num<10; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_1 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_1 += pow(2,(num));
            }
            //cout << par.at(pop).val_1 << endl;
        }
        for (int num=10; num<20; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_2 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_2 += pow(2,(num-10));
            }
        }
        for (int num=20; num<30; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_3 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_3 += pow(2,(num-20));
            }
        }
        for (int num=30; num<40; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_4 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_4 += pow(2,(num-30));
            }
        }
        for (int num=40; num<50; num++)
        {
            if (par.at(pop).combo[num] == 0)
            {
                par.at(pop).val_5 += 1;
            }
            if (par.at(pop).combo[num] == 1)
            {
                par.at(pop).val_5 += pow(2,(num-40));
            }
        }
        //cout << par.at(pop).val_1 << endl;
        //cout << par.at(pop).val_2 << endl;
        //cout << par.at(pop).val_3 << endl;
        //cout << par.at(pop).val_4 << endl;
        //cout << par.at(pop).val_5 << endl;
        par.at(pop).val.clear();
        par.at(pop).val.push_back((512-par.at(pop).val_1)/100);
        par.at(pop).val.push_back((512-par.at(pop).val_2)/100);
        par.at(pop).val.push_back((512-par.at(pop).val_3)/100);
        par.at(pop).val.push_back((512-par.at(pop).val_4)/100);
        par.at(pop).val.push_back((512-par.at(pop).val_5)/100);
        
        //cout << "x1" << "\t" << par.at(pop).val.at(0) << endl;
        //cout << "x2" << "\t" << par.at(pop).val.at(1) << endl;
        //cout << "x3" << "\t" << par.at(pop).val.at(2) << endl;
        //cout << "x4" << "\t" << par.at(pop).val.at(3) << endl;
        //cout << "x5" << "\t" << par.at(pop).val.at(4) << endl;
    }
    
    
    for(int pop=0; pop<pop_size; pop++)
    {
        par.at(pop).fitness = 0;
        double sum = 0;
        for (int v=0; v<num_val; v++)
        {
            //cout << "cp" << endl;
            //cout << par.at(pop).val.at(v) << endl;
            sum += 1*(par.at(pop).val.at(v));
        }
            //cout << "cp" << endl;
        //cout << sum << endl;
        par.at(pop).fitness = 1/(abs(sum)+1);
        //cout << "fitness" << "\t" << par.at(pop).fitness << endl;
    }
}


//-----------------------------------------------------------
//gets the fitness for each parent
void get_fitness()
{
    //for testing only
    //easy_eval();
    
    if (function_num==0)
    {
     DeJong_1();
    }
    if (function_num==1)
    {
        DeJong_2();
    }
    if (function_num==2)
    {
        DeJong_3();
    }
    if (function_num==3)
    {
        //DeJong_4();
    }
    //cout << endl;
}


//-----------------------------------------------------------
//calculates the prob of selection for each parent
void get_probability_selection()
{
    double total_fitness = 0;
    for (int pop=0; pop<pop_size; pop++)
    {
        //cout << par.at(pop).fitness << "\t";
        total_fitness += par.at(pop).fitness;
    }
    //cout << endl;
    //cout << "total fitness" << "\t" << total_fitness << endl;
    double check_prob = 0;
    for (int pop=0; pop<pop_size; pop++)
    {
        par.at(pop).prob_selection = 0;
        par.at(pop).prob_selection = par.at(pop).fitness/total_fitness;
        //cout << "prob_selcection" << "\t" << par.at(pop).prob_selection << endl;
        check_prob += par.at(pop).prob_selection;
        //cout << "parent" << "\t" << pop << "\t" << "fitness" << "\t" << par.at(pop).fitness << "\t" << "prob" << "\t" << par.at(pop).prob_selection << endl;
    }
    //cout << "check probability" << "\t" << check_prob << endl;
}


//-----------------------------------------------------------
//selects two parents based on their prob of selection then runs xover
void xover()
{
    //sort(par.begin(), par.end(), less_than_prob_selection());
    for (int s=0; s<pop_size/2; s++)
    {
        //selects parent based on its probaility of selection
        int index_1 = 0;
        double stat_1 = 0;
        double r = ((double) rand() / (RAND_MAX));
        for (int ii=0; ii<pop_size; ii++)
        {
            stat_1 = stat_1 + par.at(ii).prob_selection;
            if(stat_1 >= r)
            {
                index_1 = ii;
                break;
            }
        }
        //cout << "index 1" << "\t" << "is" << "\t" << index_1 << endl;
        int index_2 = 0;
        double stat_2 = 0;
        double rr = ((double) rand() / (RAND_MAX));
        for (int ii=0; ii<pop_size; ii++)
        {
            stat_2 = stat_2 + par.at(ii).prob_selection;
            if(stat_2 >= rr)
            {
                index_2 = ii;
                break;
            }
        }
        //will not allow index_1 and index_2 to be the same
        //while (index_1 == index_2)
        //{
            //double stat_2 = 0;
            //double rr = ((double) rand() / (RAND_MAX));
            //for (int ii=0; ii<pop_size; ii++)
            //{
                //stat_2 = stat_2 + par.at(ii).prob_selection;
                //if(stat_2 >= rr)
                //{
                    //index_2 = ii;
                    //break;
                //}
            //}
        //}
        //cout << "index 2" << "\t" << "is" << "\t" << index_2 << endl;
        //cout << endl;
        
        //copies parents genes
        for(int ii=0; ii<bit_length; ii++)
        {
            chi.at(s).combo[ii] = par.at(index_1).combo[ii];
            chi.at(s+1).combo[ii] = par.at(index_2).combo[ii];
        }
        
        //cout << "child" << "\t" << s << endl;
        //for(int ii=0; ii<bit_length; ii++)
        //{
        //cout << chi.at(s).combo[ii] << "\t";
        //}
        //cout << endl;
        //cout << "child" << "\t" << s+1 << endl;
        //for(int ii=0; ii<bit_length; ii++)
        //{
        //cout << chi.at(s+1).combo[ii] << "\t";
        //}
        //cout << endl;
        //cout << endl;
        
        
        
        
        //gene cross over
        double rrr = ((double) rand() / (RAND_MAX));
        if (rrr <= xover_prob)
        {
            int r_splice = (int)rand() % bit_length;
            //cout << "splice location" << "\t" << r_splice << endl;
            int place_holder[bit_length];
            for (int ii=0; ii<bit_length; ii++)
            {
                place_holder[ii] = par.at(index_2).combo[ii];
                //cout << place_holder[ii] << "\t";
            }
            //cout << endl;
            for (int sp=r_splice; sp<bit_length; sp++)
            {
                chi.at(s).combo[sp] = par.at(index_1).combo[sp];
                chi.at(s).combo[sp] = place_holder[sp];
            }
        }
    }
}



//-----------------------------------------------------------
//mutates the copies if the probability of mutation statement is met
void mutation(Child &M)
{
    //int bits_to_flip = (int)rand() % 1;
    int bits_to_flip = 1;
    for (int p=0; p<bits_to_flip; p++)
    {
        int bitdex = (int)rand() % bit_length;
        for (int iii=0; iii< bit_length; iii++)
        {
            if (iii == bitdex)
            {
                if (M.combo[iii] == 0)
                {
                    M.combo[iii] = 1;
                }
                else
                {
                    M.combo[iii] = 0;
                }
            }
        }
    }
}


//-----------------------------------------------------------
//builds the next generation
void build_next_gen()
{
    xover();
    for (int ii=0; ii<pop_size; ii++)
    {
        double r_mutate = ((double) rand() / (RAND_MAX));
        if (r_mutate < mutation_prob)
        {
            //cout << "mutation" << endl;
            Child M;
            M = chi.at(ii);
            //Mutates the child if the probability of mutation statement is met
            mutation(M);
        }
    }
}


//-----------------------------------------------------------
//runs the selection and child generation process
void natural_selection()
{
    get_probability_selection();
    build_next_gen();
}


//-----------------------------------------------------------
//turns children into parents
void Puberty()
{
    for (int pop=0; pop<pop_size; pop++)
    {
        for (int num=0; num<bit_length; num++)
        {
            par.at(pop).combo[num] = chi.at(pop).combo[num];
        }
    }
}


//-----------------------------------------------------------
//sorts the population based on their fitness from highest to lowest
struct greater_than_par_fitness
{
    inline bool operator() (const Parent& struct1, const Parent& struct2)
    {
        return (struct1.fitness > struct2.fitness);
    }
};


//-----------------------------------------------------------
//copies the max, ave, and min fitness
void run_scoreboard()
{
    //cout << "check" << endl;
    //cout << par.at(0).fitness <<endl;
    max_fitness.push_back(par.at(0).fitness);
    double sum = 0;
    for (int p=0; p<pop_size; p++)
    {
        sum += par.at(p).fitness;
    }
    ave_fitness.push_back(sum/pop_size);
    min_fitness.push_back(par.at(pop_size-1).fitness);
    
}


//-----------------------------------------------------------
//gets the ave of the ave
void get_ave()
{
    for (int gen=0; gen<gen_max; gen++)
    {
        double ave_max_sum = 0;
        double ave_ave_sum = 0;
        double ave_min_sum = 0;
        for (int trial=0; trial<num_trials; trial++)
        {
            ave_max_sum += max_fitness.at(gen+trial*gen_max);
        }
        ave_max_fitness = ave_max_sum/num_trials;
        ave_of_ave_max.push_back(ave_max_fitness);
        for (int trial=0; trial<num_trials; trial++)
        {
            ave_ave_sum += ave_fitness.at(gen+trial*gen_max);
        }
        ave_ave_fitness = ave_ave_sum/num_trials;
        ave_of_ave_ave.push_back(ave_ave_fitness);
        for (int trial=0; trial<num_trials; trial++)
        {
            ave_min_sum += min_fitness.at(gen+trial*gen_max);
        }
        ave_min_fitness = ave_min_sum/num_trials;
        ave_of_ave_min.push_back(ave_min_fitness);
    }
}


//-----------------------------------------------------------
//writes the max, ave, and min fitness for each generation to a txt file
void write_pop_info_all_gens_to_text()
{
    ofstream File1;
    File1.open("Max Fitness For All Gens.txt");
    ofstream File2;
    File2.open("Ave Fitness For All Gens.txt");
    ofstream File3;
    File3.open("Min Fitness For All Gens.txt");
    ofstream File4;
    File4.open("Ave Max Fitness For All Trails");
    ofstream File5;
    File5.open("Ave Ave Fitness For All Trails");
    ofstream File6;
    File6.open("Ave Min Fitness For All Trails");
    //cout << "fitness across generations" << endl;
    //cout << "max fitness" << endl;
    for (int i=0; i<gen_max*num_trials; i++)
    {
        File1 << max_fitness.at(i) << endl;
    }
    //cout << endl;
    //cout << "ave fitness" << endl;
    for (int i=0; i<gen_max*num_trials; i++)
    {
        File2 << ave_fitness.at(i) << endl;
    }
    //cout << endl;
    //cout << "min fitness" << endl;
    for (int i=0; i<gen_max*num_trials; i++)
    {
        File3 << min_fitness.at(i) << endl;
    }
    for (int i=0; i<gen_max; i++)
    {
        File4 << ave_of_ave_max.at(i) << endl;
    }
    for (int i=0; i<gen_max; i++)
    {
        File5 << ave_of_ave_ave.at(i) << endl;
    }
    for (int i=0; i<gen_max; i++)
    {
        File6 << ave_of_ave_min.at(i) << endl;
    }
    //cout << endl;
    //cout << endl;
    File1.close();
    File2.close();
    File3.close();
    File4.close();
    File5.close();
    File6.close();
}


//-----------------------------------------------------------
//runs program
void run_GA()
{
    select_prob_case();
    for (int trial=0; trial<num_trials; trial++)
    {
        cout << "RUN GA" << endl;
        cout << "-------------------------------------" << endl;
        cout << endl;
        build_pop();
        //create_target();              //for testing only
        for (int gen=0; gen<gen_max; gen++)
        {
            if (gen % 100 == 0)
            {
                cout << "generation" << "\t" << gen << endl;
                cout << endl;
            }
            get_fitness();
            sort(par.begin(), par.end(), greater_than_par_fitness());
            run_scoreboard();
            //cout << endl;
            natural_selection();
            Puberty();
            //get_fitness();
            //sort(par.begin(), par.end(), greater_than_par_fitness());
            
            if (gen == gen_max-1)
            {
                cout << "generation" << "\t" << gen << endl;
                cout << endl;
                get_fitness();
                sort(par.begin(), par.end(), greater_than_par_fitness());
                run_scoreboard();
                
                cout << "-------------------------------------" << endl;
                cout << "Final Population" << endl;
                //cout << "Population Fitness" << endl;
                //for (int p=0; p<pop_size; p++)
                //{
                //cout << par.at(p).fitness << "\t";
                //}
                //cout << endl;
                //cout << endl;
                cout << "Best Individual" << endl;
                cout << "Fitness" << endl;
                cout << par.at(0).fitness << endl;
                cout << "string" << endl;
                for (int num=0; num<bit_length/array_split; num++)
                {
                    for(int bit=0; bit<array_split; bit++)
                    {
                        cout << par.at(0).combo[(array_split*num)+bit] << "\t";
                    }
                    cout << endl;
                }
                cout << endl;
                cout << endl;
                //cout << "worst Individual" << endl;
                //cout << par.at(pop_size-1).fitness << endl;
                //cout << endl;
                //double sum = 0;
                //for (int p=0; p<pop_size; p++)
                //{
                //sum = sum + par.at(p).fitness;
                //}
                //double ave;
                //ave = sum/pop_size;
                //cout << "Average Fitness" << "\t" << ave << endl;
                //for (int num=0; num<bit_length; num++)
                //{
                //cout << A.combo[l] << "\t";
                //}
                //cout << endl;
            }
        }
    }
    get_ave();
    write_pop_info_all_gens_to_text();
}




//-----------------------------------------------------------
//-----------------------------------------------------------
int main()
{
    srand(time(NULL));
    run_GA();
    
    return 0;
}
