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
#include "genherite.hpp"
#include "evolution.hpp"
#include "integr.hpp"

void printaftergen(Cellule & cell, int netindex);

gsl_rng *rndm;

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

double fabs(const vd & vec)
{
    double res = 0;
    for (civd iv = vec.begin() + nb_steps + nb_steps / 5; iv != vec.end();
         iv++) {
        res += fabs(*iv);
    }
    return res;
}

unsigned int nb_cells = 4;

vd times_conv;
unsigned int nb_conv = 0;

bool print_mode = false;

void generations(Milieu & milieu)
{
    static unsigned int counter = 0;

    nb_promoters_max = 8;
    nb_proteins_max = 12;

    for (int ngen = 0; ngen < args_info.nb_generations_arg; ngen++) {
        milieu.evolution();
        milieu.selection();
        cout << milieu.cellules[0]->score << " " << milieu.cellules[0]->
            score_auxi << endl;
        if (milieu.cellules[0]->score < 250.0) {
            nb_promoters_max = 50;
            nb_proteins_max = 15;
        }
        if (ngen > 400 && milieu.cellules[0]->score > 9000.0)
            break;
        if (ngen % 100 == 0)
            printaftergen(*(milieu.cellules[0]), ngen);
        if (milieu.cellules[0]->score < 2.0) {
            printaftergen(*(milieu.cellules[0]), ngen);
            times_conv.push_back(ngen);
            counter++;
            break;
        }
    }
    cout << 1 << " " << milieu.cellules[0]->promoters.size() << " " << milieu.
        cellules[0]->proteins.size() << endl;
}

ofstream outintegr;
void printaftergen(Cellule & cell, int netindex)
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

int main(int argc, char **argv)
{
    /*  Parse the command line */
    if (cmdline_parser(argc, argv, &args_info) != 0)
        exit(1);
    behavior = args_info.behavior_arg;
    recept_test = 0;
    if (args_info.anchor_signal_flag) {
        anchorsignal = 1;
    } else {
        anchorsignal = 0;
    }

    /* Initialize the random number generator */
    unsigned long int seed;
    FILE *devrandom;
    rndm = gsl_rng_alloc(gsl_rng_default);
    if (args_info.rand_seed_flag) {
        devrandom = fopen("/dev/urandom", "r");
        fread(&seed, sizeof(seed), 1, devrandom);
        fclose(devrandom);
        gsl_rng_set(rndm, seed);
    }

    /* Initialise the evolutionary algorithm mode */
    nb_steps = 200;
    if (args_info.lambda_flag) {
        cout << "lambda + lambda" << endl;
        pop_size = 2 * args_info.pop_size_arg;
    } else {
        cout << "1 + lambda" << endl;
        pop_size = 1 + args_info.pop_size_arg;
    }
    nb_genes = args_info.nb_genes_arg;

    /* Start the evolution */
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
    //
    gsl_rng_free(rndm);
    return 0;
}
