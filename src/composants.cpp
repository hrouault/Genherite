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

#include <cmath>

#include "composants.hpp"
#include "cellule.hpp"

Composant::Composant(Protein * prot, int n)
:protein(prot), nb(n)
{
}

Protein::Protein(double c, double q)
:cste(c), phosphosite(NULL)
{
    qtite = q;
    meanc = frand2();
    dirinit = (int) (2 * frand());
    if (dirinit == 0)
        dirinit = -1;
}

Protein::~Protein()
{
    if (phosphosite)
        delete phosphosite;
    vector < Composant * >::iterator itercomp;
    for (itercomp = composants.begin(); itercomp != composants.end();
         itercomp++) {
        delete *itercomp;
    }
}

Gene::Gene()
{
    qtite = 1;
    arn = new Arn(this);
    transrate = frand2();
}

Arn::Arn(Gene * gen, double c, double q)
:gene(gen)
{
    cste = c;
    qtite = q;
    gene = gen;
    protein = new Protein();
    translation = new Reaction(this, NULL, this, protein);
}

Promoter::Promoter(Gene * gen, Protein * prot)
:  gene(gen), protein(prot)
{
    transrate = frand2();
    eqcte = frand2();
}

Reaction::Reaction(Espece * react1, Espece * react2, Espece * prod1,
                   Espece * prod2, double constante)
:cste(constante), reactif1(react1), reactif2(react2), produit1(prod1),
produit2(prod2)
{
}

Recepteur::Recepteur(Protein * prots, Protein * protc, Protein * protp,
                     double c, double eq)
:cste(c), a0(eq), protsource(prots), protcible(protc), protphospho(protp)
{
}

bool Protein::eqprot(Protein * prot)
{
    if (this->composants.size() != prot->composants.size())
        return false;
    vector < Composant * >::iterator icomp1, icomp2;
    for (icomp1 = this->composants.begin(); icomp1 != this->composants.end();
         icomp1++) {
        for (icomp2 = prot->composants.begin();
             icomp2 != prot->composants.end(); icomp2++) {
            if ((*icomp1)->protein == (*icomp2)->protein
                && (*icomp1)->nb == (*icomp2)->nb)
                break;
        }
        if (icomp2 == prot->composants.end())
            return false;
    }
    return true;
}

Protein *Protein::copyprot()
{
    vector < Composant * >::iterator icomp;
    Protein *prot = new Protein(cste, qtite);
    prot->dirinit = dirinit;
    copie = prot;
    prot->copie = this;
    prot->label = label;
    for (icomp = composants.begin(); icomp != composants.end(); icomp++) {
        prot->composants.push_back(new
                                   Composant(static_cast <
                                             Protein *
                                             >((*icomp)->protein->copie),
                                             (*icomp)->nb));
    }
    if (phosphosite)
        prot->phosphosite = new Phosphosite();
    return prot;
}

bool Protein::contain(Protein * prot)
{
    vector < Composant * >::iterator icomp;
    if (this == prot)
        return true;
    for (icomp = composants.begin(); icomp != composants.end(); icomp++) {
        if ((*icomp)->protein == prot) {
            return true;
        }
    }
    return false;
}

bool Espece::contain(Protein * prot, Cellule * cell)
{
    for (vector < Protein * >::iterator iprot = cell->proteins.begin();
         iprot != cell->proteins.end(); iprot++) {
        if (this == *iprot) {
            if (static_cast < Protein * >(this)->contain(prot))
                return true;
        }
    }
    return false;
}

void Protein::addcomposant(Protein * prot)
{
    vector < Composant * >::iterator icompi;
    vector < Composant * >::iterator icompf;
    if (prot->composants.empty()) {
        for (icompf = composants.begin(); icompf != composants.end(); icompf++) {
            if ((*icompf)->protein == prot) {
                (*icompf)->nb++;
                break;
            }
        }
        if (icompf == composants.end())
            composants.push_back(new Composant(prot));
    } else {
        for (icompi = prot->composants.begin();
             icompi != prot->composants.end(); icompi++) {
            for (icompf = composants.begin(); icompf != composants.end();
                 icompf++) {
                if ((*icompf)->protein == (*icompi)->protein) {
                    (*icompf)->nb += (*icompi)->nb;
                    break;
                }
            }
            if (icompf == composants.end())
                composants.push_back(new
                                     Composant((*icompi)->protein,
                                               (*icompi)->nb));
        }
    }
}

void Protein::modifdegrad(double ampli)
{
    //Modifie le taux de degradation d'une proteine quelconque
    double c = pow(ampli, 2 * frand() - 1);
    cste *= c;
    if (cste < seuil_degrad_inf)
        cste = seuil_degrad_inf;
    if (cste > seuil_degrad_sup)
        cste = seuil_degrad_sup;
}

void Promoter::modifeq(double ampli)
{
    eqcte *= pow(ampli, 2 * frand() - 1);
    if (eqcte < seuil_degrad_inf)
        eqcte = seuil_degrad_inf;
    if (eqcte > seuil_degrad_sup)
        eqcte = seuil_degrad_sup;
}

void Promoter::modiftrans(double ampli)
{
    transrate *= pow(ampli, 2 * frand() - 1);
    if (transrate < seuil_degrad_inf)
        transrate = seuil_degrad_inf;
    if (transrate > seuil_degrad_sup)
        transrate = seuil_degrad_sup;
}

void Gene::modiftrans(double ampli)
{
    transrate *= pow(ampli, 2 * frand() - 1);
    if (transrate < seuil_degrad_inf)
        transrate = seuil_degrad_inf;
    if (transrate > seuil_degrad_sup)
        transrate = seuil_degrad_sup;
}

void Protein::mutinit()
{
    // inverse l'initialisation des quantités
    dirinit = -dirinit;
}

void Arn::modifdegrad(double ampli)
{
    //Modifie le taux de degradation d'un ARN quelconque
    double c = pow(ampli, 2 * frand() - 1);
    cste *= c;
    if (cste < seuil_react_inf)
        cste = seuil_react_inf;
    if (cste > seuil_react_sup)
        cste = seuil_react_sup;
}

void Reaction::modifcinetique(double ampli)
{
    double c = pow(ampli, 2 * frand() - 1);
    cste *= c;
    if (cste < seuil_react_inf)
        cste = seuil_react_inf;
    if (cste > seuil_react_sup)
        cste = seuil_react_sup;
}

void Recepteur::modifcinetique(double ampli)
{
    double modif = 2 * frand();
    double c = pow(ampli, 2 * frand() - 1);
    if (recept_satur) {
        if (modif < 1.5) {
            cste *= c;
            if (cste < seuil_degrad_inf)
                cste = seuil_degrad_inf;
            if (cste > seuil_degrad_sup)
                cste = seuil_degrad_sup;
        } else {
            a0 *= c;
            if (a0 < seuil_degrad_inf)
                a0 = seuil_degrad_inf;
            if (a0 > seuil_degrad_sup)
                a0 = seuil_degrad_sup;
        }
    } else {
        cste *= c;
        if (cste < seuil_degrad_inf)
            cste = seuil_degrad_inf;
        if (cste > seuil_degrad_sup)
            cste = seuil_degrad_sup;
    }
}

void Protein::modifqtite(double ampli)
{
    double c = pow(ampli, 2 * frand() - 1);
    qtite *= c;
    if (qtite < seuil_qtite_inf)
        qtite = seuil_qtite_inf;
    if (qtite > seuil_qtite_sup)
        qtite = seuil_qtite_sup;
}

void Arn::modifqtite(double ampli)
{
    double c = pow(ampli, 2 * frand() - 1);
    qtite *= c;
    if (qtite < seuil_qtite_inf)
        qtite = seuil_qtite_inf;
    if (qtite > seuil_qtite_sup)
        qtite = seuil_qtite_sup;
}

void Protein::clivage(Protein & prot1, Protein & prot2)
{
    int i, size;
    vector < Composant * >::iterator icomp;
    for (icomp = composants.begin(); icomp != composants.end(); icomp++) {
        size = (*icomp)->nb;
        i = (int) ((size + 1) * frand());
        if (i) {
            prot1.composants.push_back(new Composant((*icomp)->protein, i));
        }
        if (size - i) {
            prot2.
                composants.push_back(new
                                     Composant((*icomp)->protein, size - i));
        }
    }
}

Protein *Protein::addphospho()
{
    vector < Composant * >::iterator icomp;
    Protein *prot = new Protein();
    prot->phosphosite = new Phosphosite();
    for (icomp = composants.begin(); icomp != composants.end(); icomp++) {
        prot->
            composants.push_back(new
                                 Composant((*icomp)->protein, (*icomp)->nb));
    }
    if (composants.empty())
        prot->composants.push_back(new Composant(this, 1));
    prot->label = label;
    prot->label += "*";
    return prot;
}
