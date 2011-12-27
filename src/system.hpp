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

#ifndef System_H
#define System_H

#include "composants.hpp"
#include "cellule.hpp"

class System {
  public:
    double score;
    Cellule *cellule;
     vector < Recepteur * >recepteurs;

    vd concinit;

     System();
    ~System();
    void evolution();
    void addrandrecept();
    System *copysystem();
    void calculscoresystem();
    void printintegration(string filename);
    void prtint3c(ostream & fres);
    void printgraph(ostream & out);
    void printinternpart(ostream & out);
    void printproperties(ostream & out);
    System *optimisation();
    void optievolution();
    void opticalculscore();
    bool rmrandprot();
    void printsystem(ostream & out);
    double calculscoresystem3c();
    double calculscoresystem2c();
    int essential_recept();
};

void *calcsys3c_thr(void *system);

#endif                          /* System_H */
