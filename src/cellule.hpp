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

#ifndef Cellule_H
#define Cellule_H

#include <iostream>

#include "composants.hpp"

extern int mode;

typedef vector<double> vd;
typedef vd::iterator ivd;

class Cellule
{
    public:
        double score;
        double score_auxi;
        vector<Protein *> proteins;
        vector<Gene *> genes;
        vector<Promoter *> promoters;
        vector<Arn *> arns;
        vector<Reaction *> reactions;

        Cellule();
        ~Cellule();
        Protein * addprotein(double c = frand2(), double q = frand2());
        Protein * fusion(Protein * prot1, Protein * prot2);
        Reaction * adddimere(Protein * prot1, Protein * prot2, double c = 0.2 * frand2(),
                             double q = frand2(), double cdimerisation = frand2());
        Reaction * addreaction(Espece * react1, Espece * react2, Espece * prod1,
                               Espece * prod2, double c = frand2());
        Reaction * addrandreact();
        void phosphorylation(Protein * prot);
        void addphospho();
        Reaction * addrandactiv();
        Gene * addgene();
        Promoter * addpromotion(Gene * gene, Protein * prot, double e = frand2() * 0.2,
                                double t = frand2());
        void rmreaction(Reaction * r);
        bool rmrandprot();
        bool rmrandpromo();
        bool compprot(Protein & prot1, Protein & prot2);
        void evolution();
        bool evolution_modifq();
        bool evolution_ajout();
        void optievolution();
        void calculscore();
        double scorefunction(vector<double> fonction, vector<double> resultatintegr);
        double calculscorebistable();
        double calculscoremultistable();
        double calculscoremultistable2();
        double calculscoreporte();
        double calculscoreoscill();
        double scorebioscill();
        Cellule * copycellule();
        void printgraph(ostream & out);
        void printinternpart(ostream & out);
        void printproperties(ostream & out);
        void printcellule(ostream & out);
        void printcelluleshort(ostream & out);
        void printcelluleintegr(ostream & out);
        void printintegration(ostream & out);
        //void printintegrationmultistable(ostream &out);
        void prtbioscill();
        bool addclivage();
        Cellule * optimisation();
        Protein * existe(Protein * prot);
        Protein * searchcopyprot(Protein * prot);
        void opticalculscore();
        /*double calculscoreporteinter();*/
        int especenum(Espece * pesp);
        string especestring(Espece * pesp);
        int reactionnum(Reaction * react);
        string reactionstring(Reaction * react);
        void printintegr(ostream & stre, double fint);
        void copygene(Gene * gene);
};

typedef vector<Cellule *> vcell;

typedef vcell::iterator ivcell;

void * calcsc_thr(void * cell);

#endif /* Cellule_H */
