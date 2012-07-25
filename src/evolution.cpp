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
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <pthread.h>

#include "cmdline.h"
#include "genherite.hpp"
#include "evolution.hpp"
#include "cellule.hpp"
#include "integr.hpp"

Milieu::Milieu()
{
    Cellule *cell;
    for (unsigned int i = 0; i < pop_size; i++) {
        cell = new Cellule();
        cellules.push_back(cell);
        if (use_net) {
            char lab = 'a';
            for (int j = 0; j < 2; j++) {
                cell->addgene();
                cell->genes[j]->label = lab;
                cell->genes[j]->arn->label = lab;
                cell->genes[j]->arn->protein->label = toupper(lab);
                lab++;
            }
            cell->adddimere(cell->proteins[0], cell->proteins[1]);
            cell->addpromotion(cell->genes[0], cell->proteins[1]);
        } else {
            char lab = 'a';
            for (int j = 0; j < args_info.nb_genes_arg; j++) {
                cell->addgene();
                cell->genes[j]->label = lab;
                cell->genes[j]->arn->label = lab;
                cell->genes[j]->arn->protein->label = toupper(lab);
                lab++;
            }
            if (args_info.behavior_arg == 3) {
                Protein *prot = new Protein(0, 1);
                prot->phosphosite = new Phosphosite();
                prot->composants.push_back(new Composant(prot, 1));
                prot->label = "X";
                cell->proteins.push_back(prot);
                prot = new Protein(0, 1);
                prot->phosphosite = new Phosphosite();
                prot->composants.push_back(new Composant(prot, 1));
                prot->label = "Y";
                cell->proteins.push_back(prot);
            }
            if (behavior == 5 || behavior == 7) {
                Protein *prot = new Protein(0, 1);
                prot->phosphosite = new Phosphosite();
                prot->composants.push_back(new Composant(prot, 1));
                prot->label = "X";
                cell->proteins.push_back(prot);
            }
        }
    }
}

Milieu::~Milieu()
{
    vector < Cellule * >::iterator icell;
    for (icell = cellules.begin(); icell != cellules.end(); icell++) {
        delete *icell;
    }
}

void Milieu::optimisation()
{
    vector < Cellule * >::iterator icell;
    for (icell = cellules.begin(); icell != cellules.end(); icell++) {
        (*icell)->optievolution();
    }
}

void Milieu::evolution()
{
    vector < Cellule * >::iterator icell;
    if (args_info.lambda_flag) {
        for (icell = cellules.begin() + pop_size / 2; icell != cellules.end();
             icell++) {
            (*icell)->evolution();
        }
    } else {
        for (icell = cellules.begin() + 1; icell != cellules.end(); icell++) {
            (*icell)->evolution();
        }
    }
}

void Milieu::selection()
{
    vector < Cellule * >::iterator icell;
    int count = 0;
    vector < pthread_t > idths;
    for (icell = cellules.begin(); icell != cellules.end(); icell++) {
        if ((*icell)->score < 1e-8) {
            pthread_t idth;
            Cellule *pcell = *icell;
            pthread_create(&idth, NULL, calcsc_thr, pcell);
            idths.push_back(idth);
            if (!args_info.multi_threading_flag) {
                pthread_join(idth, NULL);
            }
        }
        if ((*icell)->score > 1e8)
            count++;
    }
    if (args_info.multi_threading_flag) {
        for (vector < pthread_t >::iterator ipth = idths.begin();
             ipth != idths.end(); ipth++) {
            pthread_join(*ipth, NULL);
        }
    }
    if (args_info.simulated_annealing_flag) {
        double beta = 1.0;
        unsigned int nbcells = 0;
        vcell copy1(cellules);
        vcell copy2;
        while (nbcells < pop_size / 2) {
            vd probs;
            double z = 0;
            ivcell ic;
            for (ic = copy1.begin(); ic != copy1.end(); ic++) {
                double prob = exp(-(*ic)->score * beta + 50);
                probs.push_back(prob);
                z += prob;
            }
            int cell = 0;
            double rand_cell = z * frand();
            double cursor = 0;
            ic = copy1.begin();
            while (cursor < rand_cell) {
                cursor += probs[cell];
                cell++;
                ic++;
            }
            ic--;
            copy2.push_back(*ic);
            copy1.erase(ic);
            nbcells++;
        }
        for (ivcell ic = copy1.begin(); ic != copy1.end(); ic++) {
            delete *ic;
        }
        cellules.clear();
        for (ivcell ic = copy2.begin(); ic != copy2.end(); ic++) {
            cellules.push_back(*ic);
            cellules.push_back((*ic)->copycellule());
        }
    } else {
        std::sort(cellules.begin(), cellules.end(), compcellules);
        if (args_info.lambda_flag) {
            for (icell = cellules.begin();
                 icell != (cellules.begin() + pop_size / 2); icell++) {
                delete *(icell + pop_size / 2);
                *(icell + pop_size / 2) = (*icell)->copycellule();
            }
        } else {
            for (icell = cellules.begin() + 1; icell != cellules.end();
                 icell++) {
                delete *icell;
                *icell = (*(cellules.begin()))->copycellule();
            }
        }
    }
    //(*(cellules.begin()))->prtbioscill();
}

void Milieu::selection_temp(double temperature)
{
    vector < Cellule * >::iterator icell;
    Milieu milieu;
    for (icell = milieu.cellules.begin(); icell != milieu.cellules.end();
         icell++) {
        delete *icell;
    }
    milieu.cellules.clear();
    for (icell = cellules.begin(); icell != cellules.end(); icell++) {
        if ((*icell)->score == 0)
            (*icell)->calculscore();
    }
    double scoremin = cellules[0]->score;
    for (icell = cellules.begin(); icell != cellules.end(); icell++) {
        if ((*icell)->score < scoremin)
            scoremin = (*icell)->score;
    }
    for (unsigned int i = 0; i < pop_size / 2; i++) {
        vector < double >scores;
        for (icell = cellules.begin(); icell != cellules.end(); icell++) {
            double scorenorm = (*icell)->score - scoremin;
            scores.push_back(exp(-scorenorm / temperature));
        }
        double tot_score = 0;
        for (vector < double >::iterator iscore = scores.begin();
             iscore != scores.end(); iscore++) {
            tot_score += *iscore;
        }
        double ncell = tot_score * frand();
        int j = 0;
        double inf = 0;
        do {
            inf += scores[j];
            j++;
        } while (ncell > inf);
        j--;
        milieu.cellules.push_back(cellules[j]);
        cellules.erase(cellules.begin() + j);
    }
    for (icell = cellules.begin(); icell != cellules.end(); icell++)
        delete *icell;
    cellules.clear();
    for (icell = milieu.cellules.begin(); icell != milieu.cellules.end();
         icell++) {
        cellules.push_back(*icell);
    }
    std::sort(cellules.begin(), cellules.end(), compcellules);
    for (icell = milieu.cellules.begin(); icell != milieu.cellules.end();
         icell++) {
        cellules.push_back((*icell)->copycellule());
    }
    milieu.cellules.clear();
}

void Milieu::optiselection()
{
    vector < Cellule * >::iterator icell;
    for (icell = cellules.begin(); icell != cellules.end(); icell++) {
        (*icell)->opticalculscore();
    }
    std::sort(cellules.begin(), cellules.end(), compcellules);
    for (icell = cellules.begin(); icell != (cellules.begin() + pop_size / 2);
         icell++) {
        delete *(icell + pop_size / 2);
        *(icell + pop_size / 2) = (*icell)->copycellule();
    }
}

bool compcellules(const Cellule * cell1, const Cellule * cell2)
{
    return (cell1->score < cell2->score);
};
