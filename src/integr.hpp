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

#ifndef Integr_H
#define Integr_H

#include <vector>

#include <gsl/gsl_odeiv.h>

using namespace std;

#include "cellule.hpp"
#include "system.hpp"

class Reacteff;
class Especeeff;
class Recepteff;
class Resultatsys;
class Geneeff;
class Promoeff;
class Integrsys;

class Celleff {
  public:
    vector < Especeeff > esps;
    vector < Reacteff > reactions;
    vector < Geneeff > genes;
    vector < Recepteff > recepts;
    vd concinit;

     Celleff(Cellule & cell);
     Celleff(System & system);
    void derivqtiteeuler();
    void derivqtiterk4();
    void derivqtite();
     vector < double > integrbistable(int nprot);
     vector < double > integrbioscill(Cellule & cell, int state);
    //      vector<double> integr(int nprot);
    Resultatsys integrsys3c(System & sys);
     vector < vector < double > > integrmultistable(int prot0, int prot1,
                                                  int prot2);
     vector < double >integrporte(int nprotA, int nprotX, int nprotY,
                                  bool protX1, bool protY1, bool protX2,
                                  bool protY2);
    int integrsys3c(System & sys, Resultatsys & res, vd conci);
    int integrfunc(vd & integrA, vd & integrB, double fintime);
    int integrunifunc(double result1[], double fintime);
    double scorebiooscill(double finosc, double valstate);
    double scoremultistable();
    double scoretristagrad();
    double scorebistable();
    double scorebistable2();
    double scoremonostable();
    void prtbiooscill(ostream & out, double finosc);
    void prtmultistable(ostream & out);
    void prttristagrad(ostream & out);
    void prtbistable(ostream & out);
    void prtoscill(ostream & out);
    void prtmonostable(ostream & out);
    vd scoreoscill();
    void recoverstate(double y[]);
    void savestate(double y[]);
    void initconc(int state);
    void printlong();
    int integrgate(vd & res);
    double scoregate(bool i, bool j);
    double scorelogicgate();
    int integrgrad(double conc, vd & integrA, vd & integrB, double &score);
    void down_recept();
    void initintegr(Integrsys & intsys);
    void initstarty(double *&y);
};

class Especeeff {
  public:
    double cste;
    double qtite;
    double meanc;
    int dirinit;

     Especeeff(double c, double q);
};
typedef vector < Especeeff > vespeff;
typedef vespeff::iterator ivespeff;

class Reacteff {
  public:
    double cste;
    unsigned int r1;
    unsigned int r2;
    unsigned int p1;
    unsigned int p2;

    Reacteff(double c, unsigned int re1, unsigned int re2, unsigned int pr1,
             unsigned int pr2);
};

typedef vector < Reacteff > vreaceff;
typedef vreaceff::iterator ivreaceff;

class Geneeff {
  public:
    unsigned int arn;
    double transrate;
    vector < Promoeff > promos;

    Geneeff(unsigned int narn, double dtransrate, Gene * gn, Cellule & cell);
};
typedef vector < Geneeff > vgeff;
typedef vgeff::iterator ivgeff;

class Promoeff {
  public:
    double eqcte;
    double transrate;
    unsigned int prot;

    Promoeff(unsigned int p, double t, double e);
};

class Recepteff {
  public:
    double eqcte;
    double kphospho;
    unsigned int prots;
    unsigned int protc;
    unsigned int protp;

    Recepteff(unsigned int ps, unsigned int pc, unsigned int pp, double k,
              double e);
};
typedef vector < Recepteff > vrecepte;
typedef vrecepte::iterator ivrecepte;

class Resultatsys {
  public:
    vector < double >result1a;
    vector < double >result1b;
    vector < double >result2a;
    vector < double >result2b;
    vector < double >result3a;
    vector < double >result3b;
    vector < double >result4a;
    vector < double >result4b;
};

class Initconc {
  public:
    double cell1;
    double cell2;
    double cell3;
    double cell4;
};

class Integrsys {
  public:
    gsl_odeiv_step * step;
    gsl_odeiv_control *control;
    gsl_odeiv_evolve *evolve;
    gsl_odeiv_system system;
};

double mini(double d1, double d2);

void prtrescells(Resultatsys res, ostream & fc1a, ostream & fc1b,
                 ostream & fc1c, ostream & fc2a, ostream & fc2b,
                 ostream & fc2c, ostream & fc3a, ostream & fc3b,
                 ostream & fc3c);

void prtresult(vector < double >dbvect, ostream & file);

double scoreres(vector < double >vectdb, double value);

void prtfunctomatchbioosc(ostream & out);

double distcust(double pos1, double pos2);
double distcust(const vd & pos1, const vd & pos2);

#endif                          /* Integr_H */
