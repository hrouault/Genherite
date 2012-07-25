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
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "cmdline.h"
#include "evolution.hpp"
#include "composants.hpp"
#include "cellule.hpp"
#include "integr.hpp"

unsigned int nb_promoters_max = 5;
unsigned int nb_proteins_max = 15;

Cellule::Cellule()
:  score(0)
{
}

Cellule::~Cellule()
{
    vector < Gene * >::iterator igene;
    vector < Arn * >::iterator iarn;
    vector < Protein * >::iterator iprot;
    vector < Promoter * >::iterator ipromo;
    vector < Reaction * >::iterator ireact;
    for (igene = genes.begin(); igene != genes.end(); igene++) {
        delete *igene;
    }
    for (iarn = arns.begin(); iarn != arns.end(); iarn++) {
        delete *iarn;
    }
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        delete *iprot;
    }
    for (ipromo = promoters.begin(); ipromo != promoters.end(); ipromo++) {
        delete *ipromo;
    }
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        delete *ireact;
    }
}

Protein *Cellule::existe(Protein * prot)
{
    vector < Protein * >::iterator iprot;
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        if (prot->eqprot(*iprot)) {
            return *iprot;
        }
    }
    return NULL;
}

Protein *Cellule::searchcopyprot(Protein * prot)
{
    vector < Protein * >::iterator iprot;
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        if ((static_cast < Protein * >(*iprot)->copie) == prot)
            return *iprot;
    }
    return NULL;
}

string Cellule::especestring(Espece * pesp)
{
    if (pesp == NULL)
        return "NULL";
    vector < Gene * >::iterator igene;
    vector < Arn * >::iterator iarn;
    vector < Protein * >::iterator iprot;
    vector < Promoter * >::iterator ipromo;
    int i = 0;
    for (igene = genes.begin(); igene != genes.end(); igene++) {
        if (pesp == *igene) {
            stringstream sgene;
            sgene << "cell->genes[" << i << "]";
            return sgene.str();
        }
        i++;
    }
    i = 0;
    for (iarn = arns.begin(); iarn != arns.end(); iarn++) {
        if (pesp == *iarn) {
            stringstream sarn;
            sarn << "cell->arns[" << i << "]";
            return sarn.str();
        }
        i++;
    }
    i = 0;
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        if (pesp == *iprot) {
            stringstream sprot;
            sprot << "cell->proteins[" << i << "]";
            return sprot.str();
        }
        i++;
    }
    i = 0;
    for (ipromo = promoters.begin(); ipromo != promoters.end(); ipromo++) {
        if (pesp == *ipromo) {
            stringstream spromo;
            spromo << "cell->promoters[" << i << "]";
            return spromo.str();
        }
        i++;
    }
    cout << "Erreur : espèce introuvable!\n";
    return "-1";
}

void Cellule::copygene(Gene * gene)
{
    Gene *cgene;
    Protein *cprot;
    Espece *react1, *react2, *prod1, *prod2;
    Reaction *react;
    vector < Gene * >::iterator igene;
    vector < Protein * >::iterator iprot;
    vector < Promoter * >::iterator ipromo;
    vector < Composant * >::iterator icomp;
    vector < Reaction * >::iterator ireact;
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        (*ireact)->copie = 0;
    }
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        (*iprot)->copie = (*iprot);
    }
    char lastgene = 'a';
    for (igene = genes.begin(); igene != genes.end(); igene++) {
        lastgene++;
    }
    cgene = addgene();
    cgene->copie = gene;
    gene->copie = cgene;
    cgene->label = lastgene;
    cgene->transrate = gene->transrate;
    gene->arn->copie = cgene->arn;
    cgene->arn->copie = gene->arn;
    cgene->arn->qtite = gene->arn->qtite;
    cgene->arn->cste = gene->arn->cste;
    cgene->arn->label = lastgene;
    gene->arn->protein->copie = cgene->arn->protein;
    cgene->arn->protein->copie = gene->arn->protein;
    cgene->arn->protein->cste = gene->arn->protein->cste;
    cgene->arn->protein->qtite = gene->arn->protein->qtite;
    cgene->arn->protein->label = toupper(lastgene);
    cgene->arn->translation->cste = gene->arn->translation->cste;
    gene->arn->translation->copie = cgene->arn->translation;
    cgene->transrate = gene->transrate;
    cgene->arn->translation->copie = cgene->arn->translation;
    cprot = cgene->arn->protein;
    vector < Protein * >tempprotvec;
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        for (icomp = (*iprot)->composants.begin();
             icomp != (*iprot)->composants.end(); icomp++) {
            if ((*icomp)->protein == gene->arn->protein) {
                Protein *cpprot = (*iprot)->copyprot();
                cpprot->copie = *iprot;
                (*iprot)->copie = cpprot;
                tempprotvec.push_back(cpprot);
                break;
            }
        }
    }
    for (iprot = tempprotvec.begin(); iprot != tempprotvec.end(); iprot++) {
        Protein *prot = *iprot;
        for (string::iterator istr = prot->label.begin();
             istr != prot->label.end(); istr++) {
            if (*istr == gene->arn->protein->label[0]) {
                *istr = cgene->arn->protein->label[0];
            }
        }
        proteins.push_back(*iprot);
    }
    vector < Promoter * >temppromovec;
    for (ipromo = promoters.begin(); ipromo != promoters.end(); ipromo++) {
        Protein *prot = (*ipromo)->protein;
        if (prot->contain(gene->arn->protein)) {
            temppromovec.push_back(*ipromo);
        }
    }
    for (ipromo = temppromovec.begin(); ipromo != temppromovec.end(); ipromo++) {
        addpromotion((*ipromo)->gene,
                     static_cast < Protein * >((*ipromo)->protein->copie),
                     (*ipromo)->eqcte, (*ipromo)->transrate);
    }
    vector < Promoter * >temppromovec2;
    for (ipromo = promoters.begin(); ipromo != promoters.end(); ipromo++) {
        if ((*ipromo)->gene == gene) {
            temppromovec2.push_back(*ipromo);
        }
    }
    for (ipromo = temppromovec2.begin(); ipromo != temppromovec2.end();
         ipromo++) {
        addpromotion(cgene, static_cast < Protein * >((*ipromo)->protein),
                     (*ipromo)->eqcte, (*ipromo)->transrate);
    }
    vector < Reaction * >tempreactvec;
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        if ((*ireact)->copie == 0) {
            react = *ireact;
            bool test = false;
            if (react->reactif1->contain(gene->arn->protein, this)) {
                react1 = react->reactif1->copie;
                test = true;
            } else
                react1 = react->reactif1;
            if (react->reactif2->contain(gene->arn->protein, this)) {
                react2 = react->reactif2->copie;
                test = true;
            } else
                react2 = react->reactif2;
            if (react->produit1->contain(gene->arn->protein, this)) {
                prod1 = react->produit1->copie;
                test = true;
            } else
                prod1 = react->produit1;
            if (react->produit2->contain(gene->arn->protein, this)) {
                prod2 = react->produit2->copie;
                test = true;
            } else
                prod2 = react->produit2;
            if (test) {
                Reaction *creact =
                    new Reaction(react1, react2, prod1, prod2,
                                 (*ireact)->cste);
                tempreactvec.push_back(creact);
                creact->copie = react;
                react->copie = creact;
            }
        }
    }
    for (ireact = tempreactvec.begin(); ireact != tempreactvec.end(); ireact++) {
        reactions.push_back(*ireact);
    }
}

int Cellule::reactionnum(Reaction * react)
{
    if (react == NULL)
        return 0;
    vector < Reaction * >::iterator ireact;
    int i = 1;
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        if (react == *ireact)
            return i;
        i++;
    }
    cout << "Erreur : réaction introuvable!\n";
    return -1;
}

string Cellule::reactionstring(Reaction * react)
{
    if (react == NULL)
        return "NULL";
    vector < Reaction * >::iterator ireact;
    int i = 0;
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        if (react == *ireact) {
            stringstream sreact;
            sreact << "cell->reactions[" << i << "]";
            return sreact.str();
        }
        i++;
    }
    cout << "Erreur : réaction introuvable!\n";
    return "-1";
}

Cellule *Cellule::copycellule()
{
    Cellule *cell = new Cellule();
    Gene *gene;
    Protein *prot;
    Promoter *promo;
    Espece *react1, *react2, *prod1, *prod2;
    Reaction *react;
    vector < Gene * >::iterator igene;
    vector < Protein * >::iterator iprot;
    vector < Promoter * >::iterator ipromo;
    vector < Composant * >::iterator icomp;
    vector < Reaction * >::iterator ireact;
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        (*ireact)->copie = 0;
    }
    for (igene = genes.begin(); igene != genes.end(); igene++) {
        gene = cell->addgene();
        (*igene)->copie = gene;
        gene->copie = *igene;
        gene->label = (*igene)->label;
        (*igene)->arn->copie = gene->arn;
        gene->arn->copie = (*igene)->arn;
        gene->arn->qtite = (*igene)->arn->qtite;
        gene->arn->label = (*igene)->arn->label;
        gene->arn->cste = (*igene)->arn->cste;
        (*igene)->arn->protein->copie = gene->arn->protein;
        gene->arn->protein->copie = (*igene)->arn->protein;
        gene->arn->protein->cste = (*igene)->arn->protein->cste;
        gene->arn->protein->qtite = (*igene)->arn->protein->qtite;
        gene->arn->protein->label = (*igene)->arn->protein->label;
        gene->arn->translation->cste = (*igene)->arn->translation->cste;
        (*igene)->arn->translation->copie = gene->arn->translation;
        gene->transrate = (*igene)->transrate;
    }
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        if (!((*iprot)->composants.empty())) {
            prot = (*iprot)->copyprot();
            cell->proteins.push_back(prot);
        }
    }
    for (ipromo = promoters.begin(); ipromo != promoters.end(); ipromo++) {
        promo = new Promoter(static_cast < Gene * >((*ipromo)->gene->copie),
                             static_cast <
                             Protein * >((*ipromo)->protein->copie));
        cell->promoters.push_back(promo);
        promo->label = (*ipromo)->label;
        promo->eqcte = (*ipromo)->eqcte;
        promo->transrate = (*ipromo)->transrate;
        (*ipromo)->copie = promo;
        promo->copie = *ipromo;
    }
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        if ((*ireact)->copie == 0) {
            react = *ireact;
            (react->reactif1) ? react1 = react->reactif1->copie : react1 =
                NULL;
            (react->reactif2) ? react2 = react->reactif2->copie : react2 =
                NULL;
            (react->produit1) ? prod1 = react->produit1->copie : prod1 = NULL;
            (react->produit2) ? prod2 = react->produit2->copie : prod2 = NULL;
            Reaction *creact =
                cell->addreaction(react1, react2, prod1, prod2,
                                  (*ireact)->cste);
            creact->copie = react;
            react->copie = creact;
        }
    }
    cell->score = score;
    return cell;
}

int Cellule::especenum(Espece * pesp)
{
    if (pesp == NULL)
        return 0;
    vector < Arn * >::iterator iarn;
    vector < Protein * >::iterator iprot;
    int narns = arns.size();
    int i = 1;
    for (iarn = arns.begin(); iarn != arns.end(); iarn++) {
        if (pesp == *iarn)
            return i;
        i++;
    }
    i = 1 + narns;
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        if (pesp == *iprot)
            return i;
        i++;
    }
    cout << "Erreur : espèce introuvable (especenum)!" << pesp->label << endl;
    return -1;
}

Protein *Cellule::addprotein(double c, double q)
{
    //Rajoute une proteine exterieure dans le milieu
    Protein *prot = new Protein(c, q);
    proteins.push_back(prot);
    return prot;
}

Protein *Cellule::fusion(Protein * prot1, Protein * prot2)
{
    // rassemble des composants de i et j
    Protein *prot = new Protein();
    prot->addcomposant(prot1);
    prot->addcomposant(prot2);
    return prot;
}

Reaction *Cellule::adddimere(Protein * prot1, Protein * prot2, double c,
                             double q, double cdimerisation)
{
    vector < Composant * >::iterator icomp;
    Protein *prot = fusion(prot1, prot2);
    prot->qtite = q;
    prot->cste = c;
    proteins.push_back(prot);
    for (icomp = prot->composants.begin(); icomp != prot->composants.end();
         icomp++) {
        for (int i = 0; i < (*icomp)->nb; i++)
            prot->label += (*icomp)->protein->label;
    }
    //Ajout de la réaction de dimérisation
    return addreaction(prot1, prot2, prot, NULL, cdimerisation);
}

Reaction *Cellule::addrandactiv()
{
    Protein *anchorprot = proteins[nb_genes];
    Protein *prots;
    Protein *protc;
    if (anchorsignal) {
        prots = proteins[(int) (proteins.size() * frand())];
        protc = proteins[(int) (proteins.size() * frand())];
    } else {
        prots = anchorprot;
        protc = anchorprot;
        while (prots == anchorprot) {
            prots = proteins[(int) (proteins.size() * frand())];
        }
        while (protc == anchorprot) {
            protc = proteins[(int) (proteins.size() * frand())];
        }
    }
    Protein *protp = protc->addphospho();
    proteins.push_back(protp);
    return addreaction(prots, protc, prots, protp);
}

Reaction *Cellule::addreaction(Espece * react1, Espece * react2,
                               Espece * prod1, Espece * prod2, double q)
{
    Reaction *react = new Reaction(react1, react2, prod1, prod2, q);
    reactions.push_back(react);
    return react;
}

Reaction *Cellule::addrandreact()
{
    Protein *randprot1, *randprot2, *anchorprot;
    Reaction *react;
    anchorprot = proteins[nb_genes];
    randprot1 = randprot2 = anchorprot;
    if (anchorsignal) {
        randprot1 = proteins[(int) (frand() * proteins.size())];
        randprot2 = proteins[(int) (frand() * proteins.size())];
    } else {
        randprot1 = randprot2 = anchorprot;
        while (randprot1 == anchorprot) {
            randprot1 = proteins[(int) (frand() * proteins.size())];
        }
        while (randprot2 == anchorprot) {
            randprot2 = proteins[(int) (frand() * proteins.size())];
        }
    }
    react = adddimere(randprot1, randprot2);
    return react;
}

void Cellule::phosphorylation(Protein * prot)
{
    Protein *protcible = prot->addphospho();
    proteins.push_back(protcible);
    addreaction(prot, NULL, protcible, NULL, 0.2 * frand2());
}

void Cellule::addphospho()
{
    Protein *anchorprot = proteins[nb_genes];
    Protein *prot = anchorprot;
    if (anchorsignal) {
        prot = proteins[(int) (frand() * proteins.size())];
    } else {
        while (prot == anchorprot) {
            prot = proteins[(int) (frand() * proteins.size())];
        }
    }
    Protein *protcible = prot->addphospho();
    proteins.push_back(protcible);
    addreaction(prot, NULL, protcible, NULL);
}

Gene *Cellule::addgene()
{
    Gene *gene = new Gene();
    genes.push_back(gene);
    Arn *arn = gene->arn;
    arns.push_back(arn);
    proteins.push_back(arn->protein);
    reactions.push_back(arn->translation);
    return gene;
}

Promoter *Cellule::addpromotion(Gene * gene, Protein * prot, double e,
                                double t)
{
    Promoter *promo = new Promoter(gene, prot);
    promo->label = gene->label;
    promo->transrate = t;
    promo->eqcte = e;
    promoters.push_back(promo);
    return promo;
}

bool Cellule::rmrandpromo()
{
    if (promoters.size() == 0)
        return false;
    Promoter *promo;
    int i = (int) (frand() * promoters.size());
    promo = promoters[i];
    promoters.erase(promoters.begin() + i);
    delete promo;
    return true;
}

bool Cellule::rmrandprot()
{
    int i = (int) (proteins.size() * frand());
    Protein *prot = proteins[i];
    vector < Reaction * >reacts;
    if (prot->composants.empty() && prot->phosphosite == NULL)
        return false;
    if (prot->label == "X")
        return false;
    vector < Reaction * >::iterator ireact;
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        if ((*ireact)->reactif1 == prot || (*ireact)->reactif2 == prot)
            return false;
    }
    for (vector < Promoter * >::iterator ipromo = promoters.begin();
         ipromo != promoters.end(); ipromo++) {
        if ((*ipromo)->protein == prot)
            return false;
    }
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        if ((*ireact)->produit1 == prot || (*ireact)->produit2 == prot) {
            reacts.push_back(*ireact);
        }
    }
    for (unsigned int j = 0; j < reacts.size(); j++) {
        rmreaction(reacts.back());
        reacts.pop_back();
    }
    proteins.erase(proteins.begin() + i);
    delete prot;
    return true;
}

void Cellule::rmreaction(Reaction * r)
{
    reactions.erase(find(reactions.begin(), reactions.end(), r));
    delete r;
}

bool Cellule::evolution_modifq()
{
    const double sommetaux =
        t_mod_r + t_mod_pro + t_mod_qtite + t_rmreact + t_mod_qpromo;
    const double taux[] =
        { t_mod_r, t_mod_pro, t_mod_qtite, t_rmreact, t_mod_qpromo };
    /* Sélection de la modification */
    double r = frand() * sommetaux;
    double inf = 0;
    int i = 0;
    do {
        inf += taux[i];
        i++;
    } while (r > inf);
    i--;
    switch (i) {
    case 0:
        //Modif réaction
        reactions[(int) (frand() * reactions.size())]->modifcinetique();
        return true;
        break;
    case 1:
        //Modif dégradation
        if (frand() < 0.8)
            proteins[(int) (frand() * proteins.size())]->modifdegrad();
        else
            arns[(int) (frand() * arns.size())]->modifdegrad();
        return true;
        break;
    case 2:
        //Modif quantités
        if (frand() < 0.2)
            arns[(int) (frand() * arns.size())]->modifqtite();
        else
            proteins[(int) (frand() * proteins.size())]->modifqtite();
        return true;
        break;
    case 3:
        //Suppression de réactions
        if (frand() < 0.5) {
            if (rmrandpromo())
                return true;
        } else {
            if (rmrandprot())
                return true;
        }
        break;
    case 4:
        //Modification des quantités de promoteurs
        if (promoters.size()) {
            promoters[(int) (frand() * promoters.size())]->modifeq();
            return true;
        }
        break;
    }
    return false;
}

bool Cellule::evolution_ajout()
{
    //procedure d'evolution
    const double sommetaux = t_new_r + t_new_promo + t_phospho;
    const double taux[] = { t_new_r, t_new_promo, t_phospho };
    Gene *randgene;
    Protein *randprot;
    /* Sélection de la modification */
    double r = frand() * sommetaux;
    double inf = 0;
    int i = 0;
    do {
        inf += taux[i];
        i++;
    } while (r > inf);
    i--;
    switch (i) {
    case 0:
        //Dimérisation
        if (reactions.size() < nb_reactions_max
            && proteins.size() < nb_proteins_max) {
            addrandreact();
            return true;
        }
        break;
    case 1:
        //Ajout d'une promotion
        if (reactions.size() < nb_reactions_max - 2
            && promoters.size() < nb_promoters_max) {
            randgene = genes[(int) (frand() * genes.size())];
            randprot = proteins[(int) (frand() * proteins.size())];
            addpromotion(randgene, randprot);
            return true;
        }
        break;
    case 2:
        //Ajout d'une phosphorylation
        if (reactions.size() < nb_reactions_max
            && proteins.size() < nb_proteins_max) {
            if (frand() < 0.5)
                addphospho();
            else
                addrandactiv();
            return true;
        }
        break;
    }
    return false;
}

void Cellule::evolution()
{
    // evolutionary process
    const double sommetaux = t_new_r + t_mod_r + t_mod_pro + t_mod_qtite
        + t_new_promo + t_rmreact + t_phospho + t_mod_qpromo + t_mutdir;
    const double taux[] = { t_new_r, t_mod_r, t_mod_pro, t_mod_qtite,
        t_new_promo, t_rmreact, t_phospho, t_mod_qpromo, t_mutdir
    };
    Gene *randgene;
    Protein *randprot;
    /* Modification selection */
    int j = 1;
    while (j) {
        double r = frand() * sommetaux;
        double inf = 0;
        int i = 0;
        do {
            inf += taux[i];
            i++;
        } while (r > inf);
        i--;
        switch (i) {
        case 0:
            // Dimer formation
            if (reactions.size() < nb_reactions_max
                && proteins.size() < nb_proteins_max) {
                addrandreact();
                j = 0;
            }
            break;
        case 1:
            // Reaction parameter modification
            reactions[(int) (frand() * reactions.size())]->modifcinetique();
            j = 0;
            break;
        case 2:
            // Protein degradation modification
            if (frand() < 0.8)
                proteins[(int) (frand() * proteins.size())]->modifdegrad();
            else
                arns[(int) (frand() * arns.size())]->modifdegrad();
            j = 0;
            break;
        case 3:
            // Initial quantity modification
            if (frand() < 0.2)
                arns[(int) (frand() * arns.size())]->modifqtite();
            else
                proteins[(int) (frand() * proteins.size())]->modifqtite();
            j = 0;
            break;
        case 4:
            // New promotion
            if (reactions.size() < nb_reactions_max - 2
                && promoters.size() < nb_promoters_max) {
                randgene = genes[(int) (frand() * genes.size())];
                randprot = proteins[(int) (frand() * proteins.size())];
                addpromotion(randgene, randprot);
                j = 0;
            }
            break;
        case 5:
            // Reaction removal
            if (frand() < 0.5) {
                if (rmrandpromo())
                    j = 0;
            } else {
                if (rmrandprot())
                    j = 0;
            }
            break;
        case 6:
            // Phosphorylation
            if (reactions.size() < nb_reactions_max
                && proteins.size() < nb_proteins_max) {
                if (frand() < 0.5)
                    addphospho();
                else
                    addrandactiv();
                j = 0;
            }
            break;
        case 7:
            // Promoter parameter modification
            if (promoters.size()) {
                if (frand() < 0.5) {
                    promoters[(int) (frand() * promoters.size())]->modifeq();
                } else {
                    if (frand() < 0.5) {
                        promoters[(int) (frand() * promoters.size())]->
                            modiftrans();
                    } else
                        genes[(int) (frand() * genes.size())]->modiftrans();
                }
                j = 0;
            }
            break;
        case 8:
            //modification du sens des quantités initiales
            proteins[(int) (frand() * proteins.size())]->mutinit();
            j = 0;
            break;
        }
    }
    score = 0;
}

void Cellule::optievolution()
{
    const double sommetaux =
        t_mod_r_opti + t_mod_pro_opti + t_mod_qtite_opti + t_rmreact_opti;
    const double taux[] =
        { t_mod_r_opti, t_mod_pro_opti, t_mod_qtite_opti, t_rmreact_opti };
    /* Sélection de la modification */
    int j = 1;
    while (j) {
        double r = frand() * sommetaux;
        double inf = 0;
        int i = 0;
        do {
            inf += taux[i];
            i++;
        } while (r > inf);
        i--;
        switch (i) {
        case 0:
            //Modif réaction
            reactions[(int) (frand() * reactions.size())]->
                modifcinetique(ampli_opti);
            j = 0;
            break;
        case 1:
            //Modif dégradation
            if (frand() < 0.5) {
                proteins[(int) (frand() * proteins.size())]->
                    modifdegrad(ampli_opti);
            } else {
                arns[(int) (frand() * arns.size())]->modifdegrad(ampli_opti);
            }
            j = 0;
            break;
        case 2:
            //Modif quantités
            if (frand() < 0.3) {
                arns[(int) (frand() * arns.size())]->modifqtite(ampli_opti);
            } else {
                proteins[(int) (frand() * proteins.size())]->
                    modifqtite(ampli_opti);
            }
            j = 0;
            break;
        case 3:
            //Suppression de réactions
            if (frand() < 0.3) {
                if (rmrandpromo())
                    j = 0;
            } else {
                if (rmrandprot())
                    j = 0;
            }
            break;
        }
    }
    score = 0;
}

void Cellule::printgraph(ostream & out)
{
    printinternpart(out);
    printproperties(out);
}

void Cellule::printinternpart(ostream & out)
{
    vector < Gene * >::iterator igene;
    vector < Arn * >::iterator iarn;
    vector < Protein * >::iterator iprot;
    vector < Promoter * >::iterator ipromo;
    vector < Reaction * >::iterator ireact;
    out << "digraph network {" << endl;
    out << "fontname=\"Helvetica\"" << endl;
    out << "node [fontname=\"Helvetica\"];" << endl;
    out << "edge [fontname=\"Helvetica\"]" << endl;
    out << "{ rank=same; node [shape=box,color=red];" << endl;
    for (iarn = arns.begin(); iarn != arns.end(); iarn++) {
        out << "\"" << *iarn << "\" [label=\"" << (*iarn)->
            label << "\\nρ=" << setprecision(2) << (*iarn)->gene->
            transrate << ", δ=" << (*iarn)->cste << ", α=" << (*iarn)->
            translation->cste << "\"]; ";
    }
    out << "}" << endl;
    out << "{" << endl;
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        out << "\"" << *iprot << "\" [label=\"" << (*iprot)->
            label << ", δ=" << setprecision(2) << (*iprot)->cste << "\"]; ";
    }
    out << "}" << endl;
    out << "node [shape=plaintext];" << endl;
    for (iarn = arns.begin(); iarn != arns.end(); iarn++) {
        out << "\"" << *iarn << "\" -> \"" << (*iarn)->protein << "\"; ";
    }
    for (ipromo = promoters.begin(); ipromo != promoters.end(); ipromo++) {
        out << "\"" << (*ipromo)->protein << "\" -> \"" << (*ipromo)->gene->
            arn << "\" [color=red, decorate=true, label=\"K=" <<
            setprecision(2) << (*ipromo)->eqcte << "\\nρ=" << (*ipromo)->
            transrate << "\"";
        if ((*ipromo)->transrate > (*ipromo)->gene->transrate) {
            out << "];";
        } else {
            out << ", arrowhead=\"tee\"];";
        }
    }
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++)
        (*ireact)->copie = NULL;
    for (igene = genes.begin(); igene != genes.end(); igene++) {
        (*igene)->arn->translation->copie = (*igene)->arn->translation;
    }
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        if ((*ireact)->reactif1 && (*ireact)->reactif2 && (*ireact)->produit1
            && (*ireact)->produit2
            && (*ireact)->reactif1 == (*ireact)->produit1) {
            out << "\"" << *ireact << "\" [label=\"" << setprecision(2) <<
                (*ireact)->cste << "\"];";
            out << "\"" << (*ireact)->
                reactif2 << "\" -> \"" << *ireact << "\"; ";
            out << "\"" << *ireact << "\" -> \"" << (*ireact)->
                produit2 << "\"; ";
            out << "\"" << (*ireact)->
                reactif1 << "\" -> \"" << *ireact << "\" [style=dashed]; ";
            (*ireact)->copie = *ireact;
        }
    }
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        if ((*ireact)->copie == NULL) {
            out << "\"" << *ireact << "\" [label=\"" << setprecision(2) <<
                (*ireact)->cste << "\"];";
            if ((*ireact)->reactif1) {
                out << "\"" << (*ireact)->
                    reactif1 << "\" -> \"" << *ireact << "\"; ";
            }
            if ((*ireact)->reactif2) {
                out << "\"" << (*ireact)->
                    reactif2 << "\" -> \"" << *ireact << "\"; ";
            }
            if ((*ireact)->produit1) {
                out << "\"" << *ireact << "\" -> \"" << (*ireact)->
                    produit1 << "\"; ";
            }
            if ((*ireact)->produit2) {
                out << "\"" << *ireact << "\" -> \"" << (*ireact)->
                    produit2 << "\"; ";
            }
        }
    }
    out << endl;
}

void Cellule::printproperties(ostream & out)
{
    out << "label=\"Gene\ttranscription rate\\n";
    for (ivgene ivg = genes.begin(); ivg != genes.end(); ivg++) {
        out << (*ivg)->label << "\t" << (*ivg)->transrate << "\\n";
    }
    out << "\\n";
    out << "ARN\tinitial concentration\tdegradation" << "\\n";
    for (ivarn iva = arns.begin(); iva != arns.end(); iva++) {
        out << (*iva)->label << "\t" << (*iva)->qtite << "\t" << (*iva)->
            cste << "\\n";
    }
    out << "\\n";
    out << "Protein/Complex name\tinitial quantity\tdecay rate" << "\\n";
    for (ivprotein ivprot = proteins.begin(); ivprot != proteins.end();
         ivprot++) {
        out << (*ivprot)->label << "\t" << (*ivprot)->
            qtite << "\t" << (*ivprot)->cste << "\\n";
    }
    out << "\\n";
    out << "Regulated gene Regulator\tEq const.\tTranscription rate" << "\\n";
    for (ivpromoter ivpromo = promoters.begin(); ivpromo != promoters.end();
         ivpromo++) {
        out << (*ivpromo)->protein->label << "\t" << (*ivpromo)->gene->
            label << "\t" << (*ivpromo)->eqcte << "\t" << (*ivpromo)->
            transrate << "\\n";
    }
    out << "\\n";
    out << "Reaction\trate" << "\\n";
    for (ivreaction ivr = reactions.begin(); ivr != reactions.end(); ivr++) {
        string labr1;
        string labr2;
        string labp1;
        string labp2;
        labr1 = (*ivr)->reactif1->label;
        if ((*ivr)->reactif2)
            labr2 = (*ivr)->reactif2->label;
        else
            labr2 = "";
        labp1 = (*ivr)->produit1->label;
        if ((*ivr)->produit2)
            labp2 = (*ivr)->produit2->label;
        else
            labp2 = "";
        out << labr1 << "\t" << labr2 << "\t->\t" << labp1 << "\t" << labp2 <<
            "\t" << (*ivr)->cste << "\\n";
    }
    out << "\\n";
    if (args_info.behavior_arg != 2 && args_info.behavior_arg != 22) {
        out << "\"\\n}\\n";
    }
}

void Cellule::printintegration(ostream & out)
{
    vector < double >resultat;
    Celleff ceff = Celleff(*this);
    switch (args_info.behavior_arg) {
    case 0:
        resultat = ceff.integrbistable(especenum(proteins[0]));
        break;
    case 1:
        //         resultat=ceff.integr(especenum(proteins[0]));
        break;
    case 3:
        resultat =
            ceff.integrporte(especenum(proteins[0]),
                             especenum(proteins[args_info.nb_genes_arg]),
                             especenum(proteins[args_info.nb_genes_arg + 1]),
                             false, false, false, false);
        break;
    }
    vector < double >::iterator idouble;
    for (idouble = resultat.begin(); idouble != resultat.end(); idouble++) {
        out << *idouble << endl;
    }
}

void Cellule::printcellule(ostream & out)
{
    vector < Gene * >::iterator igene;
    vector < Arn * >::iterator iarn;
    vector < Protein * >::iterator iprot;
    vector < Promoter * >::iterator ipromo;
    vector < Composant * >::iterator icomp;
    vector < Reaction * >::iterator ireact;
    //   for (igene=genes.begin();igene!=genes.end();igene++){
    //      out << "gene " << (*igene)->label << endl ;
    //   }
    //   for (iarn=arns.begin();iarn!=arns.end();iarn++){
    //      out << "arn " << (*iarn)->label << " " << (*iarn)->qtite << " " << (*iarn)->cste << " ";
    //      out << especenum((*iarn)->gene) << " " << especenum((*iarn)->protein) << " ";
    //      out << reactionnum((*iarn)->translation) << endl;
    //   }
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        out << "prot " << (*iprot)->label << " " << (*iprot)->
            qtite << " " << (*iprot)->cste << " ";
        if ((*iprot)->phosphosite)
            out << " 1 ";
        out << endl;
    }
    //   for (ipromo=promoters.begin();ipromo!=promoters.end();ipromo++){
    //      out << "promo " << (*ipromo)->label << " " << (*ipromo)->qtite << " ";
    //      out << especenum((*ipromo)->gene) << " " << especenum((*ipromo)->protein) << " ";
    //   }
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        out << "react " << (*ireact)->cste << " ";
        out << especenum((*ireact)->reactif1) << " " << especenum((*ireact)->
                                                                  reactif2) <<
            " ";
        out << especenum((*ireact)->produit1) << " " << especenum((*ireact)->
                                                                  produit2) <<
            endl;
    }
}

void Cellule::printcelluleshort(ostream & out)
{
    vector < Gene * >::iterator igene;
    vector < Arn * >::iterator iarn;
    vector < Protein * >::iterator iprot;
    vector < Promoter * >::iterator ipromo;
    vector < Composant * >::iterator icomp;
    vector < Reaction * >::iterator ireact;
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        out << (*iprot)->label << endl;
    }
    for (ipromo = promoters.begin(); ipromo != promoters.end(); ipromo++) {
        out << (*ipromo)->label << endl;
    }
}

void Cellule::printcelluleintegr(ostream & out)
{
    vector < Gene * >::iterator igene;
    vector < Arn * >::iterator iarn;
    vector < Protein * >::iterator iprot;
    vector < Promoter * >::iterator ipromo;
    vector < Composant * >::iterator icomp;
    vector < Reaction * >::iterator ireact;
    out << "Gene *gene;\n";
    for (igene = genes.begin(); igene != genes.end(); igene++) {
        out << "gene=new Gene();\n";
        out << "cell->genes.push_back(gene);\n";
        out << "gene->label=\"" << (*igene)->label << "\";\n";
        out << "gene->qtite=" << (*igene)->qtite << ";\n";
    }
    out << "Arn *arn;\n";
    for (iarn = arns.begin(); iarn != arns.end(); iarn++) {
        out << "arn=new Arn(" << especestring((*iarn)->gene) << ");\n";
        out << "cell->arns.push_back(arn);\n";
        out << "arn->label=\"" << (*iarn)->label << "\";\n";
        out << "arn->qtite=" << (*iarn)->qtite << ";\n";
        out << "arn->cste=" << (*iarn)->cste << ";\n";
    }
    out << "Protein *prot;\n";
    for (iprot = proteins.begin(); iprot != proteins.end(); iprot++) {
        out << "prot=new Protein();\n";
        out << "cell->proteins.push_back(prot);\n";
        out << "prot->label=\"" << (*iprot)->label << "\";\n";
        out << "prot->qtite=" << (*iprot)->qtite << ";\n";
        out << "prot->cste=" << (*iprot)->cste << ";\n";
        if ((*iprot)->phosphosite)
            out << "prot->phosphosite=new Phosphosite();\n";
        for (icomp = (*iprot)->composants.begin();
             icomp != (*iprot)->composants.end(); icomp++) {
            out << "prot->composants.push_back(new Composant(" <<
                especestring((*icomp)->protein);
            out << "," << (*icomp)->nb << "));\n";
        }
    }
    out << "Promoter *promo;\n";
    for (ipromo = promoters.begin(); ipromo != promoters.end(); ipromo++) {
        out << "promo=new Promoter(" << especestring((*ipromo)->gene);
        out << "," << especestring((*ipromo)->protein) << ");\n";
        out << "cell->promoters.push_back(promo);\n";
        out << "promo->label=\"" << (*ipromo)->label << "\";\n";
        out << "promo->qtite=" << (*ipromo)->qtite << ";\n";
    }
    out << "Reaction *react;\n";
    for (ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        out << "react=new Reaction(";
        out << especestring((*ireact)->
                            reactif1) << "," << especestring((*ireact)->
                                                             reactif2) << ",";
        out << especestring((*ireact)->
                            produit1) << "," << especestring((*ireact)->
                                                             produit2);
        out << "," << (*ireact)->cste << ");\n";
        out << "cell->reactions.push_back(react);\n";
    }
    int i = 0;
    for (igene = genes.begin(); igene != genes.end(); igene++) {
        out << "cell->genes[" << i << "]->arn=" << especestring((*igene)->
                                                                arn) << ";\n";
        i++;
    }
    i = 0;
    for (iarn = arns.begin(); iarn != arns.end(); iarn++) {
        out << "cell->arns[" << i << "]->protein=" << especestring((*iarn)->
                                                                   protein) <<
            ";\n";
        out << "cell->arns[" << i << "]->translation=" <<
            reactionstring((*iarn)->translation) << ";\n";
        i++;
    }
    i = 0;
    for (ipromo = promoters.begin(); ipromo != promoters.end(); ipromo++) {
        i++;
    }
}

void *calcsc_thr(void *cell)
{
    ((Cellule *)cell) -> calculscore();
    return 0;
}

void Cellule::calculscore()
{
    double s = 0;               //score
    double sca = 0;
    Celleff ceff = Celleff(*this);
    vd scores;
    switch (args_info.behavior_arg) {
    case 0:
        s = ceff.scorebistable2();
        //         s+=1.0*proteins.size();
        //         s+=1.0*promoters.size();
        break;
    case 1:
        scores = ceff.scoreoscill();
        s = scores[0];
        sca = scores[1];
        //         s+=1.0*proteins.size();
        //         s+=1.0*promoters.size();
        break;
    }
    if (s != s)
        score = 1.0e10;
    else {
        score = s;
        score_auxi = sca;
    }
}

void Cellule::printintegr(ostream & stre, double fint)
{
    Celleff ceff = Celleff(*this);
    double *result = new double[nb_steps];
    //ceff.integrfunc(result,fint);
    for (unsigned int i = 0; i < nb_steps; i++) {
        stre << result[i] << endl;
    }
    delete[]result;
}

void Cellule::opticalculscore()
{
    double s = 0;               //score
    calculscore();
    s += 5 * promoters.size();
    s += 5 * proteins.size();
    score += s;
}

double Cellule::scorefunction(vector < double >fonction,
                              vector < double >resultatintegr)
{
    double s = 0;
    vector < double >::iterator idbres, idbf;
    for (idbres = resultatintegr.begin(); idbres != resultatintegr.end();
         idbres++) {
        s += gsl_pow_2(*idbres - *idbf);
    }
    return s;
}

Cellule *Cellule::optimisation()
{
    //procedure d'optimisation finale (cstes de réaction et suppression des réactions inutiles
    Milieu milieu;
    vector < Cellule * >::iterator icell;
    for (icell = milieu.cellules.begin(); icell != milieu.cellules.end();
         icell++) {
        delete *icell;
    }
    milieu.cellules.clear();
    for (int i = 0; i < args_info.pop_size_arg; i++) {
        milieu.cellules.push_back(copycellule());
    }
    for (int j = 0; j < nb_generations_opti; j++) {
        milieu.optimisation();
        milieu.optiselection();
        cout << "Generation opti " << j << endl;
        cout << milieu.cellules[0]->score << endl;
        milieu.cellules[0]->printcelluleshort(cout);
        //milieu.cellules[0]->printcelluleintegr(cout);
    }
    return milieu.cellules[0]->copycellule();
}
