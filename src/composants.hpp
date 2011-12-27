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

#ifndef Composants_H
#define Composants_H

#include <vector>
#include <string>

using namespace std;

#include "const.hpp"
#include "genherite.hpp"

class Composant;
class Arn;
class Reaction;
class Protein;
class Cellule;

class Espece
{
        //type Espece chimique: protéïne, gène ou arn
    public:
        double qtite;
        Espece * copie;
        string label;

        bool contain(Protein * prot, Cellule * cell);
};

class Phosphosite
{
};

class Protein : public Espece
{
        //type Protéïne : constante de dégradation et composants (null si un seul)
    public:
        double cste;
        vector<Composant *> composants;
        Phosphosite * phosphosite;
        double meanc;
        int dirinit;

        Protein(double c = frand2(), double q = frand2());
        ~Protein();
        void modifdegrad(double ampli = defampli);
        void addcomposant(Protein * prot);
        void modifqtite(double ampli = defampli);
        void clivage(Protein & prot1, Protein & prot2);
        bool eqprot(Protein * prot);
        Protein * addphospho();
        Protein * copyprot();
        bool contain(Protein * prot);
        void mutinit();
};

class Composant
{
    public:
        Protein * protein;
        int nb;

        Composant(Protein * prot, int n = 1);
};

class Gene : public Espece
{
    public:
        Arn * arn;
        double transrate;

        Gene();
        void modiftrans(double ampli = defampli);
};

class Arn : public Espece
{
    public:
        double cste;
        Gene * gene;
        Protein * protein;
        Reaction * translation;

        Arn(Gene * gen, double c = frand2(), double q = frand2());
        void modifdegrad(double ampli = defampli);
        void modifqtite(double ampli = defampli);
};

class Reaction
{
        //Reaction chimique:Une constante, deux reactifs,deux produits; null s'ils manquent
    public:
        double cste;
        Espece * reactif1;
        Espece * reactif2;
        Espece * produit1;
        Espece * produit2;
        Reaction * copie;

        Reaction(Espece * react1, Espece * react2, Espece * prod1, Espece * prod2, double constante = frand2());
        void modifcinetique(double ampli = defampli);
};

class Promoter : public Espece
{
    public:
        Gene * gene;
        Protein * protein;
        double eqcte;
        double transrate;

        Promoter(Gene * gen, Protein * prot);
        void modifeq(double ampli = defampli);
        void modiftrans(double ampli = defampli);
};

class Recepteur
{
    public:
        double cste;
        double a0;
        Protein * protsource;
        Protein * protcible;
        Protein * protphospho;

        Recepteur(Protein * prots, Protein * protc, Protein * protp, double c = frand2(), double eqcte = frand2());
        void modifcinetique(double ampli = frand2());
};

typedef vector<Gene *> vgene;
typedef vector<Arn *> varn;
typedef vector<Protein *> vprotein;
typedef vector<Reaction *> vreaction;
typedef vector<Promoter *> vpromoter;
typedef vector<Recepteur *> vrecept;

typedef vgene::iterator ivgene;
typedef varn::iterator ivarn;
typedef vprotein::iterator ivprotein;
typedef vreaction::iterator ivreaction;
typedef vpromoter::iterator ivpromoter;
typedef vrecept::iterator ivrecept;

#endif /* Composants_H */
