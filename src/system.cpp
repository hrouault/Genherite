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
#include <fstream>
#include <string>
#include <gsl/gsl_math.h>

using namespace std;

#include "cmdline.h"
#include "evolution.hpp"
#include "system.hpp"
#include "integr.hpp"

System::System()
    : score(0)
{
    Protein * prot;
    //   if (anchorsignal){
    prot = new Protein(0, 1);
    prot->phosphosite = new Phosphosite();
    prot->composants.push_back(new Composant(prot, 1));
    prot->label = "X";
    //   }
    cellule = new Cellule();
    char lab = 'a';
    for (int i = 0; i < args_info.nb_genes_arg; i++) {
        cellule->addgene();
        cellule->genes[i]->label = lab;
        cellule->genes[i]->arn->label = lab;
        cellule->genes[i]->arn->protein->label = toupper(lab);
        lab++;
    }
    //   if (anchorsignal){
    cellule->proteins.push_back(prot);
    //   }
    //cellule->addpromotion(cellule->genes[0],cellule->proteins[0],1.6,10.9,11.0);
    //cellule->adddimere(cellule->proteins[0],cellule->proteins[1],1,1,20);
}

System::~System()
{
    vector<Recepteur *>::iterator irecept;
    for (irecept = recepteurs.begin(); irecept != recepteurs.end(); irecept++) {
        delete *irecept;
    }
    delete cellule;
}

void System::addrandrecept()
{
    Protein * prots;
    Protein * protc;
    Protein * anchorprot = cellule->proteins[nb_genes];
    if (anchorsignal) {
        prots = cellule->proteins[(int)(cellule->proteins.size() * frand())];
        protc = cellule->proteins[(int)(cellule->proteins.size() * frand())];
    } else {
        prots = anchorprot;
        protc = anchorprot;
        while (prots == anchorprot) {
            prots = cellule->proteins[(int)(cellule->proteins.size() * frand())];
        }
        while (protc == anchorprot) {
            protc = cellule->proteins[(int)(cellule->proteins.size() * frand())];
        }
    }
    //   while (protc==prots || protc==anchorprot){
    //      protc=cellule->proteins[(int)(cellule->proteins.size()*frand())];
    //   }
    Protein * protp = protc->addphospho();
    cellule->proteins.push_back(protp);
    recepteurs.push_back(new Recepteur(prots, protc, protp));
}

void System::evolution()
{
    const double sommetaux = t_new_r + t_mod_r + t_mod_pro + t_mod_qtite
                             + t_new_promo + t_clivage + t_rmreact + t_phospho + t_mod_qpromo + t_dupgene + t_mutdir + t_new_recept + t_mod_recept;
    const double taux[] = {t_new_r, t_mod_r, t_mod_pro, t_mod_qtite,
                           t_new_promo, t_clivage, t_rmreact, t_phospho, t_mod_qpromo, t_dupgene, t_mutdir, t_new_recept, t_mod_recept
                          };
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
                //Dimérisation
                if (cellule->reactions.size() < nb_reactions_max && cellule->proteins.size() < nb_proteins_max) {
                    cellule->addrandreact();
                    j = 0;
                }
                break;
            case 1:
                //Modif réaction
                cellule->reactions[(int)(frand()*cellule->reactions.size())]->modifcinetique();
                j = 0;
                break;
            case 2:
                //Modif dégradation
                if (frand() < 0.8) cellule->proteins[(int)(frand()*cellule->proteins.size())]->modifdegrad();
                else cellule->arns[(int)(frand()*cellule->arns.size())]->modifdegrad();
                j = 0;
                break;
            case 3:
                //Modif quantités
                if (frand() < 0.2) cellule->arns[(int)(frand()*cellule->arns.size())]->modifqtite();
                else cellule->proteins[(int)(frand()*cellule->proteins.size())]->modifqtite();
                j = 0;
                break;
            case 4:
                //Ajout d'une promotion
                if (cellule->reactions.size() < nb_reactions_max - 2 && cellule->promoters.size() < nb_promoters_max) {
                    Gene * randgene = cellule->genes[(int)(frand() * cellule->genes.size())];
                    Protein * anchorprot = cellule->proteins[nb_genes];
                    Protein * randprot;
                    if (anchorsignal) {
                        randprot = cellule->proteins[(int)(frand() * cellule->proteins.size())];
                    } else {
                        randprot = anchorprot;
                        while (randprot == anchorprot) {
                            randprot = cellule->proteins[(int)(frand() * cellule->proteins.size())];
                        }
                    }
                    cellule->addpromotion(randgene, randprot);
                    j = 0;
                }
                break;
            case 5:
                //Clivage d'une protéïne
                if (cellule->reactions.size() < nb_reactions_max && cellule->proteins.size() < nb_proteins_max - 1) {
                    cellule->addclivage();
                    j = 0;
                }
                break;
            case 6:
                //Suppression de réactions
                if (frand() < 0.5) {
                    if (cellule->rmrandpromo()) j = 0;
                } else {
                    if (rmrandprot()) j = 0;
                }
                break;
            case 7:
                //Ajout d'une phosphorylation
                if (cellule->reactions.size() < nb_reactions_max && cellule->proteins.size() < nb_proteins_max) {
                    if (frand() < 0.5) cellule->addphospho();
                    else cellule->addrandactiv();
                    j = 0;
                }
                break;
            case 8:
                //Modification des quantités de promoteurs
                if (cellule->promoters.size()) {
                    if (frand() < 0.5) {
                        cellule->promoters[(int)(frand()*cellule->promoters.size())]->modifeq();
                    } else {
                        if (frand() < 0.5) {
                            cellule->promoters[(int)(frand()*cellule->promoters.size())]->modiftrans();
                        } else cellule->genes[(int)(frand()*cellule->genes.size())]->modiftrans();
                    }
                    j = 0;
                }
                break;
            case 9:
                //Duplication d'un gène
                cellule->copygene(cellule->genes[(int)(frand()*cellule->genes.size())]);
                j = 0;
                break;
            case 10:
                //modification du sens des quantités initiales
                cellule->proteins[(int)(frand()*cellule->proteins.size())]->mutinit();
                j = 0;
                break;
            case 11:
                // Ajout d'un récepteur
                if (recepteurs.size() < nb_recepteurs_max) {
                    addrandrecept();
                    j = 0;
                }
                break;
            case 12:
                // Modification d'un récepteur
                if (recepteurs.size()) {
                    recepteurs[(int)(frand()*recepteurs.size())]->modifcinetique();
                    j = 0;
                }
                break;
        }
    }
    score = 0;
}

void System::printgraph(ostream & out)
{
    printinternpart(out);
    printproperties(out);
}

void System::printinternpart(ostream & out)
{
    vector<Recepteur *>::iterator irecept;
    cellule->printinternpart(out);
    out << "node [shape=plaintext,style=dotted,color=blue];" << endl;
    for (irecept = recepteurs.begin(); irecept != recepteurs.end(); irecept++) {
        out << "\"" << *irecept << "\" [label=\"" << setprecision(2) << (*irecept)->cste << "\"];";
        out << "\"" << (*irecept)->protsource << "\" -> \"" << *irecept << "\" [style=dotted,color=blue]; ";
        out << "\"" << (*irecept)->protcible << "\" -> \"" << *irecept << "\" [style=dashed,color=blue]; ";
        out << "\"" << *irecept << "\" -> \"" << (*irecept)->protphospho << "\" [style=dashed,color=blue]; ";
    }
    out << "\n";
}

void System::printproperties(ostream & out)
{
    cellule->printproperties(out);
    out << "Receptor\trate" << "\\n";
    for (ivrecept ivr = recepteurs.begin(); ivr != recepteurs.end(); ivr++) {
        out << (*ivr)->protsource->label << "\t" << (*ivr)->protcible->label << "\t" << (*ivr)->cste << "\\n";
    }
    out << "\\n\"" << endl;
    out << "}" << endl;
}

void System::prtint3c(ostream & fres)
{
    Resultatsys res;
    vector<double>::iterator idb;
    Celleff sys = Celleff(*this);
    print_mode = true;
    sys.integrsys3c(*this, res, concinit);
    print_mode = false;
    for (unsigned int i = 0; i < 2 * nb_steps; i++) {
        fres << i << " " << res.result1a[i] << " ";
        fres << res.result2a[i] << " ";
        fres << res.result3a[i] << " ";
        if (nb_cells > 3) {
            fres << res.result4a[i] << " ";
        }
        fres << res.result1b[i] << " ";
        fres << res.result2b[i] << " ";
        fres << res.result3b[i] << " ";
        if (nb_cells > 3) {
            fres << res.result4b[i] << " ";
        }
        fres << endl;
    }
}


void System::printsystem(ostream & out)
{
    vector<Recepteur *>::iterator irecept;
    cellule->printcellule(out);
    for (irecept = recepteurs.begin(); irecept != recepteurs.end(); irecept++) {
        out << "recept " << (*irecept)->cste << " ";
        out << cellule->especenum((*irecept)->protsource) << " ";
        out << cellule->especenum((*irecept)->protcible) << " ";
        out << cellule->especenum((*irecept)->protphospho) << endl;
    }
}

void System::opticalculscore()
{
    calcsys3c_thr(this);
    double s = score;
    s += coef_react * cellule->reactions.size();
    s += coef_prot * cellule->proteins.size();
    s += coef_promo * cellule->promoters.size();
    s += coef_recept * recepteurs.size();
    if (s != s) {
        score = 1.0e10;
    } else {
        score = s;
    }
}

void * calcsys3c_thr(void * system)
{
    System * psys = (System *)system;
    System & sys = *psys;
    if (behavior == 2) {
        sys.calculscoresystem3c();
    } else if (behavior == 22) {
        sys.calculscoresystem2c();
    }
    return 0;
}

double System::calculscoresystem3c()
{
    score = 0;
    Resultatsys res;
    Celleff sys = Celleff(*this);
    print_mode = true;
    int flag = sys.integrsys3c(*this, res, concinit);
    print_mode = false;
    if (!flag) {
        //      for (unsigned int i=nb_steps/4;i<2*nb_steps;i++){//nb_steps;i<2*nb_steps;i++){
        //         score+=gsl_pow_2(10.0-res.result1a[i]);
        //         score+=gsl_pow_2(10.0-res.result2b[i]);
        //         score+=gsl_pow_2(0.1-res.result1b[i]);
        //         score+=gsl_pow_2(0.1-res.result2a[i]);
        //         score+=gsl_pow_2(0.1-res.result3a[i]);
        //         score+=gsl_pow_2(0.1-res.result3b[i]);
        //         score+=gsl_pow_2(0.1-res.result4a[i]);
        //         score+=gsl_pow_2(0.1-res.result4b[i]);
        //      }
        //      score+=10*cellule->promoters.size();
        //      score+=10*cellule->proteins.size();
        //      score/=50.0;
        //      score=0;
        for (unsigned int i = nb_steps; i < 2 * nb_steps; i++) { //nb_steps;i<2*nb_steps;i++){
            double dA12 = distcust(res.result1a[i], res.result2a[i]);
            double dB12 = distcust(res.result1b[i], res.result2b[i]);
            double dA13 = distcust(res.result1a[i], res.result3a[i]);
            double dB13 = distcust(res.result1b[i], res.result3b[i]);
            double dA23 = distcust(res.result2a[i], res.result3a[i]);
            double dB23 = distcust(res.result2b[i], res.result3b[i]);
            double dinvA34 = 1 - distcust(res.result3a[i], res.result4a[i]);
            double dinvB34 = 1 - distcust(res.result3b[i], res.result4b[i]);
            score += mini(dA12, dB12) + dinvA34 + dinvB34 + mini(dA13, dB13) + mini(dA23, dB23);
        }
        for (unsigned int i = nb_steps / 2; i < nb_steps; i++) { //nb_steps;i<2*nb_steps;i++){
            score += 0.5 * fabs(res.result1a[nb_steps + i] - res.result1a[3 * nb_steps + i]);
            score += 0.5 * fabs(res.result2a[nb_steps + i] - res.result2a[3 * nb_steps + i]);
            score += 0.5 * fabs(res.result3a[nb_steps + i] - res.result3a[3 * nb_steps + i]);
            score += 0.5 * fabs(res.result1b[nb_steps + i] - res.result1b[3 * nb_steps + i]);
            score += 0.5 * fabs(res.result2b[nb_steps + i] - res.result2b[3 * nb_steps + i]);
            score += 0.5 * fabs(res.result3b[nb_steps + i] - res.result3b[3 * nb_steps + i]);
        }
    } else {
        score = 1e10;
    }
    if (opti_stage) {
        score += coef_react * cellule->reactions.size();
        score += coef_prot * cellule->proteins.size();
        score += coef_promo * cellule->promoters.size();
        score += coef_recept * recepteurs.size();
    }
    if (isnan(score)) cout << "NaA observed!!" << endl;
    if (print_mode) {
        double dA12 = distcust(res.result1a[2 * nb_steps - 1], res.result2a[2 * nb_steps - 1]);
        double dB12 = distcust(res.result1b[2 * nb_steps - 1], res.result2b[2 * nb_steps - 1]);
        double dA13 = distcust(res.result1a[2 * nb_steps - 1], res.result3a[2 * nb_steps - 1]);
        double dB13 = distcust(res.result1b[2 * nb_steps - 1], res.result3b[2 * nb_steps - 1]);
        double dA23 = distcust(res.result2a[2 * nb_steps - 1], res.result3a[2 * nb_steps - 1]);
        double dB23 = distcust(res.result2b[2 * nb_steps - 1], res.result3b[2 * nb_steps - 1]);
        double dinvA34 = 1 - distcust(res.result3a[2 * nb_steps - 1], res.result4a[2 * nb_steps - 1]);
        double dinvB34 = 1 - distcust(res.result3b[2 * nb_steps - 1], res.result4b[2 * nb_steps - 1]);
        score = mini(dA12, dB12) + dinvA34 + dinvB34 + mini(dA13, dB13) + mini(dA23, dB23);
    }
    return score;
}

double System::calculscoresystem2c()
{
    score = 0;
    Resultatsys res;
    Celleff sys = Celleff(*this);
    //   print_mode=true;
    int flag = sys.integrsys3c(*this, res, concinit);
    //   print_mode=false;
    if (!flag) {
        for (unsigned int i = nb_steps / 2; i < 2 * nb_steps; i++) { //nb_steps;i<2*nb_steps;i++){
            double dA12 = distcust(res.result1a[i], res.result2a[i]);
            score += dA12;
        }
    } else {
        score = 1e10;
    }
    if (opti_stage) {
        score += coef_react * cellule->reactions.size();
        score += coef_prot * cellule->proteins.size();
        score += coef_promo * cellule->promoters.size();
        score += coef_recept * recepteurs.size();
    }
    if (print_mode) {
        return fabs(res.result1a[2 * nb_steps - 1] - res.result2a[2 * nb_steps - 1]);
    }
    return score;
}

int System::essential_recept()
{
    score = 0;
    Resultatsys res;
    Celleff sys = Celleff(*this);
    print_mode = true;
    sys.integrsys3c(*this, res, concinit);
    print_mode = false;
    //   print_mode=true;
    if (fabs(res.result1a[2 * nb_steps - 1] - res.result2a[2 * nb_steps - 1]) > 0.2) return 1;
    //   print_mode=false;
    return 0;
}

System * System::copysystem()
{
    System * psys = new System();
    delete psys->cellule;
    vector<Recepteur *>::iterator irecept;
    psys->concinit.push_back(concinit[0]);
    psys->concinit.push_back(concinit[1]);
    psys->cellule = cellule->copycellule();
    for (irecept = recepteurs.begin(); irecept != recepteurs.end(); irecept++) {
        Protein * prots = static_cast<Protein *>((*irecept)->protsource->copie);
        Protein * protc = static_cast<Protein *>((*irecept)->protcible->copie);
        Protein * protp = static_cast<Protein *>((*irecept)->protphospho->copie);
        double c = (*irecept)->cste;
        Recepteur * precept = new Recepteur(prots, protc, protp, c);
        psys->recepteurs.push_back(precept);
    }
    psys->score = score;
    return psys;
}

bool System::rmrandprot()
{
    int i = (int)(cellule->proteins.size() * frand());
    Protein * prot = cellule->proteins[i];
    vector<Reaction *> reacts;
    if (prot->composants.empty() && prot->phosphosite == NULL) return false;
    if (prot->label == "X") return false;
    vector<Reaction *>::iterator ireact;
    for (ireact = cellule->reactions.begin(); ireact != cellule->reactions.end(); ireact++) {
        if ((*ireact)->reactif1 == prot || (*ireact)->reactif2 == prot) return false;
    }
    for (vector<Promoter *>::iterator ipromo = cellule->promoters.begin(); ipromo != cellule->promoters.end(); ipromo++) {
        if ((*ipromo)->protein == prot) return false;
    }
    vector<Recepteur *>::iterator irecept, delrecept;
    for (irecept = recepteurs.begin(); irecept != recepteurs.end(); irecept++) {
        if ((*irecept)->protsource == prot || (*irecept)->protcible == prot) return false;
    }
    for (ireact = cellule->reactions.begin(); ireact != cellule->reactions.end(); ireact++) {
        if ((*ireact)->produit1 == prot || (*ireact)->produit2 == prot) {
            reacts.push_back(*ireact);
        }
    }
    for (unsigned int j = 0; j < reacts.size(); j++) {
        cellule->rmreaction(reacts.back());
        reacts.pop_back();
    }
    delrecept = recepteurs.end();
    for (irecept = recepteurs.begin(); irecept != recepteurs.end(); irecept++) {
        if ((*irecept)->protphospho == prot) {
            delrecept = irecept;
            break;
        }
    }
    if (delrecept != recepteurs.end()) {
        delete *delrecept;
        recepteurs.erase(delrecept);
    }
    cellule->proteins.erase(cellule->proteins.begin() + i);
    delete prot;
    return true;
}

void System::optievolution()
{
    const double sommetaux = t_mod_r_opti + t_mod_pro_opti + t_mod_qtite_opti + t_rmreact_opti + t_mod_recept_opti + t_mod_qpromo;
    const double taux[] = {t_mod_r_opti, t_mod_pro_opti, t_mod_qtite_opti, t_rmreact_opti, t_mod_recept_opti, t_mod_qpromo};
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
                cellule->reactions[(int)(frand()*cellule->reactions.size())]->modifcinetique(ampli_opti);
                j = 0;
                break;
            case 1:
                //Modif dégradation
                if (frand() < 0.8) cellule->proteins[(int)(frand()*cellule->proteins.size())]->modifdegrad(ampli_opti);
                else cellule->arns[(int)(frand()*cellule->arns.size())]->modifdegrad(ampli_opti);
                j = 0;
                break;
            case 2:
                //Modif quantités
                if (frand() < 0.2) cellule->arns[(int)(frand()*cellule->arns.size())]->modifqtite(ampli_opti);
                else cellule->proteins[(int)(frand()*cellule->proteins.size())]->modifqtite(ampli_opti);
                j = 0;
                break;
            case 3:
                //Suppression de réactions
                if (frand() < 0.3) {
                    if (cellule->rmrandpromo()) j = 0;
                } else {
                    if (rmrandprot()) j = 0;
                }
                break;
            case 4:
                if (recepteurs.size()) {
                    recepteurs[(int)(frand()*recepteurs.size())]->modifcinetique();
                    j = 0;
                    break;
                }
            case 5:
                //Modification des quantités de promoteurs
                if (cellule->promoters.size()) {
                    if (frand() < 0.5) {
                        cellule->promoters[(int)(frand()*cellule->promoters.size())]->modifeq();
                    } else {
                        if (frand() < 0.5) {
                            cellule->promoters[(int)(frand()*cellule->promoters.size())]->modiftrans();
                        } else cellule->genes[(int)(frand()*cellule->genes.size())]->modiftrans();
                    }
                    j = 0;
                }
                break;
        }
    }
    score = 0;
}

System * System::optimisation()
{
    //procedure d'optimisation finale (cstes de réaction et suppression des réactions inutiles
    Milieu_System milieu;
    vector<System *>::iterator isys;
    for (isys = milieu.systems.begin(); isys != milieu.systems.end(); isys++) delete *isys;
    milieu.systems.clear();
    for (unsigned int i = 0; i < pop_size; i++) {
        milieu.systems.push_back(copysystem());
    }
    for (int j = 0; j < nb_generations_opti; j++) {
        milieu.optimisation();
        milieu.selection3c();
    }
    return milieu.systems[0]->copysystem();
}
