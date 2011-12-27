/*
 * Copyright (C) 2005-2011 Hervé Rouault <rouault@lps.ens.fr>
 *
 * This file is part of Genherite.
 *
 * Genherite is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Genherite is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Genherite.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>

#include <gsl/gsl_rng.h>

using namespace std;


#include "cmdline.h"
#include "cellule.hpp"
#include "system.hpp"
#include "genherite.hpp"
#include "evolution.hpp"
#include "integr.hpp"


int mode = 0;

void printaftergen(Cellule & cell, int netindex);

void printaftergen(System & sys, int netindex, int mode);

gsl_rng * rndm;

double frand2()
{
    return exp(2 * gsl_rng_uniform(rndm) - 1.0);
}

double frand()
{
    return gsl_rng_uniform(rndm);
}

vd operator+(const vd & vec1, const vd & vec2)
{
    vd vec;
    civd iv2 = vec2.begin();
    for (civd iv1 = vec1.begin(); iv1 != vec1.end(); iv1++) {
        vec.push_back(*iv1 + *iv2);
        iv2++;
    }
    return vec;
}

vd operator-(const vd & vec1, const vd & vec2)
{
    vd vec;
    civd iv2 = vec2.begin();
    for (civd iv1 = vec1.begin(); iv1 != vec1.end(); iv1++) {
        vec.push_back(*iv1 - *iv2);
        iv2++;
    }
    return vec;
}

double
fabs(const vd & vec)
{
    double res = 0;
    for (civd iv = vec.begin() + nb_steps + nb_steps / 5; iv != vec.end(); iv++) {
        res += fabs(*iv);
    }
    return res;
}

unsigned int nb_cells = 4;


vd times_conv;
unsigned int nb_conv = 0;

bool print_mode = false;

void
generations(Milieu & milieu)
{
    mode = 1;
    static unsigned int counter = 0;
    nb_promoters_max = 8;
    nb_proteins_max = 12;
    for (int ngen = 0; ngen < args_info.nb_generations_arg; ngen++) {
        milieu.evolution();
        milieu.selection();
        cout << milieu.cellules[0]->score << " " << milieu.cellules[0]->score_auxi << endl;
        if (milieu.cellules[0]->score < 250.0) {
            nb_promoters_max = 50;
            nb_proteins_max = 15;
        }
        if (ngen > 400 && milieu.cellules[0]->score > 9000.0) break;
        if (ngen % 100 == 0) printaftergen(*(milieu.cellules[0]), ngen);
        if (milieu.cellules[0]->score < 2.0) {
            printaftergen(*(milieu.cellules[0]), ngen);
            times_conv.push_back(ngen);
            counter++;
            break;
        }
    }
    cout << 1 << " " << milieu.cellules[0]->promoters.size() << " " << milieu.cellules[0]->proteins.size() << endl;
}

int
generations_system(Milieu_System & milieu)
{
    mode = 1;
    //   static unsigned int counter=0;
    nb_promoters_max = 8;
    nb_proteins_max = 12;
    //   nb_promoters_max=15;
    //   nb_proteins_max=20;
    if (behavior == 2) {
        nb_cells = 4;
    } else if (behavior == 22) {
        nb_cells = 2;
    }
    for (int ngen = 0; ngen < args_info.nb_generations_arg; ngen++) {
        milieu.evolution();
        //      milieu.evolution();
        milieu.selection3c();
        //      nb_cells=4;
        //      if (ngen%50==0 && ngen>0){
        //         printaftergen(*(milieu.systems[0]),ngen,0);
        //      }
        //      nb_cells=4;
        //      cout<<"Generation "<<ngen<<endl;
        //      cout << milieu.systems[0]->score << endl;
        if (milieu.systems[0]->score < 1.0) {
            opti_stage = 1;
            printaftergen(*(milieu.systems[0]), ngen, 1);
            System * sys = milieu.systems[0]->optimisation();
            opti_stage = 0;
            cout << "opti done" << endl;
            printaftergen(*sys, ngen, 2);
            int resrecept = 0;
            recept_test = 1;
            if (nb_cells == 2 && sys->essential_recept()) resrecept++;
            if (nb_cells == 2) printaftergen(*(sys), 0, 3);
            print_mode = 1;
            recept_test = 0;
            double final_test = 0;
            if (nb_cells == 2) {
                final_test = sys->calculscoresystem2c();
            } else {
                final_test = sys->calculscoresystem3c();
            }
            print_mode = 0;
            //         delete sys;
            if (nb_cells == 2) {
                if (final_test > 1.0) {
                    return 1 + resrecept;
                } else {
                    return 0;
                }
            } else {
                if (final_test < 0.1) {
                    return 1;
                } else {
                    return 0;
                }
            }
        }
        //      milieu.systems[0]->printsystem(cout);
        /*      if (milieu.systems[0]->score < 1200.0){
                 stringstream sout;
                 sout << "goodintegr" << counter << ".dat";
                 ofstream fout(sout.str().c_str());
                 Celleff ceff(*(milieu.systems[0]));
                 ceff.prtmultistable(fout);
                 fout.close();

                 stringstream sout2;
                 sout2 << "goodcell" << counter << ".dot";
                 ofstream fout2(sout2.str().c_str());
                 milieu.systems[0]->printgraph(sout2);
                 fout2.close();

                 System * sys=milieu.systems[0]->optimisation();
                 for (int i=0;i<taille_milieu;i++){
                    delete milieu.systems[i];
                 }
                 milieu.systems[0]=sys;
                 for (int i=1;i<taille_milieu;i++){
                    milieu.systems[i]=milieu.systems[0]->copysystem();
                 }
                 mode=1;
                 counter++;
                 break;
              }*/
    }
    return 0;
}

ofstream outintegr;
void
printaftergen(Cellule & cell, int netindex)
{
    stringstream sout;
    sout << "graph" << netindex << ".dot";
    ofstream fout(sout.str().c_str());
    cell.printgraph(fout);
    fout.close();
    stringstream sout2;
    sout2 << "integr" << netindex << ".dat";
    outintegr.open(sout2.str().c_str());
    print_mode = true;
    cell.calculscore();
    print_mode = false;
    outintegr.close();
    stringstream sout3;
    sout3 << "cell" << netindex << ".dat";
    ofstream fout3(sout3.str().c_str());
    cell.printcelluleintegr(fout3);
    fout3.close();
}

void
printaftergen(System & sys, int netindex, int mode)
{
    stringstream sout;
    string suffix;
    if (mode == 1) {
        suffix = "-befconv-";
    } else if (mode == 2) {
        suffix = "-afterconv-";
    } else if (mode == 3) {
        suffix = "-recept-";
    } else {
        suffix = "";
    }
    sout << "graph" << suffix << netindex << ".dot";
    ofstream fout(sout.str().c_str());
    sys.printgraph(fout);
    fout.close();
    stringstream sout2;
    sout2 << "integr" << suffix << netindex << ".dat";
    outintegr.open(sout2.str().c_str());
    print_mode = true;
    if (behavior == 22) {
        sys.calculscoresystem2c();
    } else if (behavior == 2) {
        sys.calculscoresystem3c();
    }
    print_mode = false;
    outintegr.close();
    //  sys.prtint3c(fout2);
    /*  double state[ceff.esps.size()];
      ceff.savestate(state);
      ceff.prtbiooscill(fout2,final_time);
      ceff.recoverstate(state);
      ceff.prtbiooscill(fout2,final_time-2.5);*/
    /*
      stringstream sout3;
      sout3 << "system" << netindex << ".dat";
      ofstream fout3(sout3.str().c_str());
      sys.printsystem(fout3);
      fout3.close();
    */
}

struct gengetopt_args_info args_info;

unsigned int nb_steps;
unsigned int pop_size;
unsigned int nb_genes;
unsigned int behavior;
unsigned int recept_satur;
unsigned int anchorsignal;

unsigned int recept_test;

unsigned int opti_stage = 0;

double t_new_recept = 0;

int main(int argc, char ** argv)
{
    /* Initialisation du générateur de nombres aléatoires */
    if (cmdline_parser(argc, argv, & args_info) != 0)
        exit(1);
    behavior = args_info.behavior_arg;
    recept_test = 0;
    if (args_info.anchor_signal_flag) {
        anchorsignal = 1;
    } else {
        anchorsignal = 0;
    }
    unsigned long int seed;
    FILE * devrandom;
    rndm = gsl_rng_alloc(gsl_rng_default);
    if (args_info.rand_seed_flag) {
        devrandom = fopen("/dev/urandom", "r");
        fread(&seed, sizeof(seed), 1, devrandom);
        fclose(devrandom);
        gsl_rng_set(rndm, seed);
    }
    if (args_info.recept_satur_flag) {
        recept_satur = 1;
    } else {
        recept_satur = 0;
    }
    nb_steps = 200;
    if (args_info.lambda_flag) {
        cout << "lambda + lambda" << endl;
        pop_size = 2 * args_info.pop_size_arg;
    } else {
        cout << "1 + lambda" << endl;
        pop_size = 1 + args_info.pop_size_arg;
    }
    nb_genes = args_info.nb_genes_arg;
    /*
    if (behavior==2){
       ofstream scores("scores.lst");
       for (int i=0;i<args_info.nb_net_arg;i++){
          Milieu_System milieu;
          stringstream comsave;
          comsave << "./save.sh " << i;
          cout << comsave.str() << "\n";
          generations_system(milieu);
          system(comsave.str().c_str());
          scores << milieu.systems[0]->score << "\n";
       }
       scores.close();
    } else if (behavior==22){
    */
    if (behavior == 2 || behavior == 22) {
        ofstream scores("scores.lst");
        vd interact_proba;
        interact_proba.push_back(0); //r1
        //      interact_proba.push_back(0.05); //r2
        //      interact_proba.push_back(0.1); //r3
        //      interact_proba.push_back(0.25); //r4
        //      interact_proba.push_back(0.5); //r5
        //      interact_proba.push_back(1.0); //r6
        //      interact_proba.push_back(2.0); //r7
        //      interact_proba.push_back(5.0); //r8
        //      interact_proba.push_back(10.0); //r9
        //      interact_proba.push_back(20.0); //r10
        for (ivd iv = interact_proba.begin(); iv != interact_proba.end(); iv++) {
            t_new_recept = *iv;
            unsigned int nb_converged = 0;
            unsigned int nb_essrecept = 0;
            for (int i = 0; i < args_info.nb_net_arg; i++) {
                Milieu_System milieu;
                stringstream comsave;
                comsave << "./save.sh -" << *iv << "-" << i;
                cout << comsave.str() << "\n";
                int result = generations_system(milieu);
                cout << milieu.systems[0]->score << endl;
                unsigned int ess_recept = 0;
                unsigned int converged = 0;
                if (result) {
                    converged = 1;
                    nb_converged++;
                    if (result == 2) {
                        ess_recept++;
                        nb_essrecept++;
                    }
                }
                comsave << "-" << ess_recept;
                if (converged) {
                    system(comsave.str().c_str());
                } else {
                    system("rm *.dat *.dot");
                }
                scores << milieu.systems[0]->score << endl;
                //         if (converged);
            }
            cout << t_new_recept << " " << (double)nb_essrecept / (double)nb_converged << " " << nb_converged << endl;
        }
        scores.close();
    } else {
        ofstream scores("scores.lst");
        for (int i = 0; i < args_info.nb_net_arg; i++) {
            Milieu milieu;
            generations(milieu);
            if (milieu.cellules[0]->score < 50.0) {
                stringstream comsave;
                comsave << "./save.sh " << i;
                cout << comsave.str() << "\n";
                system(comsave.str().c_str());
            }
            scores << milieu.cellules[0]->score << "\n";
        }
        scores.close();
        //         printaftergen(*(milieu.cellules[0]),i);
    }
    gsl_rng_free(rndm);
    return 0;
}
