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

#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "cmdline.h"
#include "integr.hpp"

double
mini(double d1, double d2)
{
    if (d1 < d2) return d1;
    else return d2;
}


Celleff::Celleff(System & system)
{
    Cellule * cell = system.cellule;
    int nbtot = cell->arns.size() + cell->proteins.size();
    vector<Gene *>::iterator igene;
    vector<Arn *>::iterator iarn;
    vector<Protein *>::iterator iprot;
    vector<Promoter *>::iterator ipromo;
    vector<Reaction *>::iterator ireact;
    vector<Recepteur *>::iterator irecept;
    Especeeff dummy = Especeeff(0, 1);
    esps.push_back(dummy);
    for (igene = cell->genes.begin(); igene != cell->genes.end(); igene++) {
        Gene & gen = *(*igene);
        Geneeff gn = Geneeff(cell->especenum(gen.arn), gen.transrate, *igene, *cell);
        genes.push_back(gn);
    }
    for (unsigned int i = 0; i < nb_cells; i++) {
        for (iarn = cell->arns.begin(); iarn != cell->arns.end(); iarn++) {
            Especeeff esp = Especeeff((*iarn)->cste, (*iarn)->qtite);
            esps.push_back(esp);
        }
        for (iprot = cell->proteins.begin(); iprot != cell->proteins.end(); iprot++) {
            Especeeff esp = Especeeff((*iprot)->cste, (*iprot)->qtite);
            esps.push_back(esp);
        }
    }
    for (ireact = cell->reactions.begin(); ireact != cell->reactions.end(); ireact++) {
        double c = (*ireact)->cste;
        unsigned int nr1 = cell->especenum((*ireact)->reactif1);
        unsigned int nr2 = cell->especenum((*ireact)->reactif2);
        unsigned int np1 = cell->especenum((*ireact)->produit1);
        unsigned int np2 = cell->especenum((*ireact)->produit2);
        for (unsigned int i = 0; i < nb_cells; i++) {
            unsigned int r1;
            unsigned int r2;
            unsigned int p1;
            unsigned int p2;
            if (nr1) r1 = nr1 + i * nbtot;
            else r1 = 0;
            if (nr2) r2 = nr2 + i * nbtot;
            else r2 = 0;
            if (np1) p1 = np1 + i * nbtot;
            else p1 = 0;
            if (np2) p2 = np2 + i * nbtot;
            else p2 = 0;
            reactions.push_back(Reacteff(c, r1, r2, p1, p2));
        }
    }
    for (irecept = system.recepteurs.begin(); irecept != system.recepteurs.end(); irecept++) {
        double c = (*irecept)->cste;
        double eqcste = (*irecept)->a0;
        int nps = cell->especenum((*irecept)->protsource);
        int npc = cell->especenum((*irecept)->protcible);
        int npp = cell->especenum((*irecept)->protphospho);
        if (behavior == 2) {
            // Statut différent de la cellule centrale!
            //
            unsigned int ps;
            unsigned int pc;
            unsigned int pp;
            ps = nps;
            pc = npc + nbtot;
            pp = npp + nbtot;
            //         reactions.push_back(Reacteff(c,ps,pc,ps,pp));
            recepts.push_back(Recepteff(ps, pc, pp, c, eqcste));
            ps = nps + nbtot;
            pc = npc;
            pp = npp;
            //         reactions.push_back(Reacteff(2*c,ps,pc,ps,pp));
            recepts.push_back(Recepteff(ps, pc, pp, 2 * c, eqcste));
            for (unsigned int i = 1; i < nb_cells - 1; i++) {
                unsigned int ps2;
                unsigned int pc2;
                unsigned int pp2;
                ps2 = nps + i * nbtot;
                pc2 = npc + (i + 1) * nbtot;
                pp2 = npp + (i + 1) * nbtot;
                //            reactions.push_back(Reacteff(c,ps2,pc2,ps2,pp2));
                recepts.push_back(Recepteff(ps2, pc2, pp2, c, eqcste));
                ps2 = nps + (i + 1) * nbtot;
                pc2 = npc + i * nbtot;
                pp2 = npp + i * nbtot;
                //            reactions.push_back(Reacteff(c,ps2,pc2,ps2,pp2));
                recepts.push_back(Recepteff(ps2, pc2, pp2, c, eqcste));
            }
        } else if (behavior == 22) {
            for (unsigned int i = 0; i < nb_cells - 1; i++) {
                unsigned int ps2;
                unsigned int pc2;
                unsigned int pp2;
                ps2 = nps + i * nbtot;
                pc2 = npc + (i + 1) * nbtot;
                pp2 = npp + (i + 1) * nbtot;
                //            reactions.push_back(Reacteff(c,ps2,pc2,ps2,pp2));
                recepts.push_back(Recepteff(ps2, pc2, pp2, c, eqcste));
                ps2 = nps + (i + 1) * nbtot;
                pc2 = npc + i * nbtot;
                pp2 = npp + i * nbtot;
                //            reactions.push_back(Reacteff(c,ps2,pc2,ps2,pp2));
                recepts.push_back(Recepteff(ps2, pc2, pp2, c, eqcste));
            }
        }
    }
}

int
func(double t, const double y[], double f[], void * params)
{
    Celleff & cell = *((Celleff *)params);
    const int nbesp = cell.esps.size();
    for (int i = 1; i < nbesp; i++) {
        f[i] = -cell.esps[i].cste * y[i];
    }
    for (vector<Reacteff>::iterator ireact = cell.reactions.begin(); ireact != cell.reactions.end(); ireact++) {
        const Reacteff & react = *ireact;
        double deriv = react.cste * y[react.r1] * y[react.r2];
        f[react.r1] -= deriv;
        f[react.r2] -= deriv;
        f[react.p1] += deriv;
        f[react.p2] += deriv;
    }
    for (vector<Geneeff>::iterator igene = cell.genes.begin(); igene != cell.genes.end(); igene++) {
        Geneeff & gene = *igene;
        double sum = 1.0;
        double rate = gene.transrate;
        for (vector<Promoeff>::iterator ipromo = gene.promos.begin(); ipromo != gene.promos.end(); ipromo++) {
            double conc = y[(*ipromo).prot];
            sum += ipromo->eqcte * pow(conc, hill);
            rate += ipromo->transrate * (*ipromo).eqcte * pow(conc, hill);
        }
        f[gene.arn] += rate / sum;
    }
    if (behavior == 2 || behavior == 22) {
        if (recept_satur) {
            for (ivrecepte ire = cell.recepts.begin(); ire != cell.recepts.end(); ire++) {
                //            double cons2=y[ire->protc]*y[ire->protc];
                //            double a02=ire->eqcte*ire->eqcte;
                //            double deriv=ire->kphospho*cons2/a02/(1+cons2/a02)*y[ire->protc];
                //            f[ire->protc]-=deriv;
                //            f[ire->protp]+=deriv;
                double cons2 = y[ire->protc];
                double deriv = ire->kphospho * ire->eqcte * y[ire->prots] * y[ire->prots] / (1 + ire->eqcte * y[ire->prots] * y[ire->prots]) * y[ire->protc];
                f[ire->protc] -= deriv;
                f[ire->protp] += deriv;
            }
        } else {
            for (ivrecepte ire = cell.recepts.begin(); ire != cell.recepts.end(); ire++) {
                double deriv = ire->kphospho * y[ire->prots] * y[ire->protc];
                f[ire->protc] -= deriv;
                f[ire->protp] += deriv;
            }
        }
    }
    /* morphogen gradient */
    if (behavior == 7) {
        f[2 * nb_genes + 1] = 0;
    }
    /* logic gates input */
    if (behavior == 3) {
        f[2 * nb_genes + 1] = 0;
        f[2 * nb_genes + 2] = 0;
    }
    f[0] = 0;
    return GSL_SUCCESS;
}

int
funcsys(double t, const double y[], double f[], void * params)
{
    Celleff & cell = *((Celleff *)params);
    const int nbesp = cell.esps.size();
    const int nbtot = (nbesp - 1) / nb_cells;
    for (int i = 1; i < nbesp; i++) {
        f[i] = -cell.esps[i].cste * y[i];
    }
    for (vector<Reacteff>::iterator ireact = cell.reactions.begin(); ireact != cell.reactions.end(); ireact++) {
        const Reacteff & react = *ireact;
        double deriv = react.cste * y[react.r1] * y[react.r2];
        f[react.r1] -= deriv;
        f[react.r2] -= deriv;
        f[react.p1] += deriv;
        f[react.p2] += deriv;
    }
    for (vector<Geneeff>::iterator igene = cell.genes.begin(); igene != cell.genes.end(); igene++) {
        for (unsigned int i = 0; i < nb_cells; i++) {
            Geneeff & gene = *igene;
            double sum = 1.0;
            double rate = gene.transrate;
            for (vector<Promoeff>::iterator ipromo = gene.promos.begin(); ipromo != gene.promos.end(); ipromo++) {
                double conc = y[(*ipromo).prot + i * nbtot];
                sum += (*ipromo).eqcte * pow(conc, hill);
                rate += (*ipromo).transrate * (*ipromo).eqcte * pow(conc, hill);
            }
            f[gene.arn + i * nbtot] += rate / sum;
        }
    }
    if (behavior == 2 || behavior == 22) {
        if (recept_satur) {
            for (ivrecepte ire = cell.recepts.begin(); ire != cell.recepts.end(); ire++) {
                //            double cons2=y[ire->protc]*y[ire->protc];
                //            double a02=ire->eqcte*ire->eqcte;
                //            double deriv=ire->kphospho*cons2/a02/(1+cons2/a02)*y[ire->protc];
                //            f[ire->protc]-=deriv;
                //            f[ire->protp]+=deriv;
                double cons2 = y[ire->protc];
                double deriv = ire->kphospho * ire->eqcte * y[ire->prots] * y[ire->prots] / (1 + ire->eqcte * y[ire->prots] * y[ire->prots]) * y[ire->protc];
                f[ire->protc] -= deriv;
                f[ire->protp] += deriv;
            }
        } else {
            for (ivrecepte ire = cell.recepts.begin(); ire != cell.recepts.end(); ire++) {
                double deriv = ire->kphospho * y[ire->prots] * y[ire->protc];
                f[ire->protc] -= deriv;
                f[ire->protp] += deriv;
            }
        }
    }
    f[0] = 0;
    return GSL_SUCCESS;
}

int
funcx3c(double t, const double y[], double f[], void * params)
{
    Celleff & cell = *((Celleff *)params);
    int nbtot = (cell.esps.size() - 1) / nb_cells;
    funcsys(t, y, f, params);
    f[2 * nb_genes + 1] = 0;
    f[nbtot + 2 * nb_genes + 1] = 0;
    if (nb_cells > 2) {
        f[2 * nbtot + 2 * nb_genes + 1] = 0;
        if (nb_cells > 3) {
            f[3 * nbtot + 2 * nb_genes + 1] = 0;
        }
    }
    return GSL_SUCCESS;
}

int
jac(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
    cout << "calc jac" << endl;
    Celleff & cell = *((Celleff *) params);
    const int nbesp = cell.esps.size();
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, nbesp, nbesp);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set_zero(m);
    for (int i = 0; i < nbesp; i++) {
        gsl_matrix_set(m, i, i, -cell.esps[i].cste);
    }
    for (vector<Reacteff>::iterator ireact = cell.reactions.begin(); ireact != cell.reactions.end(); ireact++) {
        const Reacteff & react = *ireact;
        double mval = gsl_matrix_get(m, react.r1, react.r1);
        gsl_matrix_set(m, react.r1, react.r1, mval - react.cste * y[react.r2]);
        mval = gsl_matrix_get(m, react.r1, react.r2);
        gsl_matrix_set(m, react.r1, react.r2, mval - react.cste * y[react.r1]);
        mval = gsl_matrix_get(m, react.r2, react.r1);
        gsl_matrix_set(m, react.r2, react.r1, mval - react.cste * y[react.r2]);
        mval = gsl_matrix_get(m, react.r2, react.r2);
        gsl_matrix_set(m, react.r2, react.r2, mval - react.cste * y[react.r1]);
        mval = gsl_matrix_get(m, react.p1, react.r1);
        gsl_matrix_set(m, react.p1, react.r1, mval + react.cste * y[react.r2]);
        mval = gsl_matrix_get(m, react.p1, react.r2);
        gsl_matrix_set(m, react.p1, react.r2, mval + react.cste * y[react.r1]);
        mval = gsl_matrix_get(m, react.p2, react.r1);
        gsl_matrix_set(m, react.p2, react.r1, mval + react.cste * y[react.r2]);
        mval = gsl_matrix_get(m, react.p2, react.r2);
        gsl_matrix_set(m, react.p2, react.r2, mval + react.cste * y[react.r1]);
    }
    for (vector<Geneeff>::iterator igene = cell.genes.begin(); igene != cell.genes.end(); igene++) {
        Geneeff & gene = *igene;
        double sum = 1.0;
        double rate = gene.transrate;
        for (vector<Promoeff>::iterator ipromo = gene.promos.begin(); ipromo != gene.promos.end(); ipromo++) {
            double conc = y[(*ipromo).prot];
            sum += (*ipromo).eqcte * pow(conc, hill);
            rate += (*ipromo).transrate * (*ipromo).eqcte * pow(conc, hill);
        }
        for (vector<Promoeff>::iterator ipromo = gene.promos.begin(); ipromo != gene.promos.end(); ipromo++) {
            double mval = gsl_matrix_get(m, gene.arn, (*ipromo).prot);
            gsl_matrix_set(m, gene.arn, (*ipromo).prot, mval + ((*ipromo).eqcte * hill * pow(y[(*ipromo).prot], (hill - 1))) / sum * ((*ipromo).transrate - rate / sum));
        }
    }
    for (int i = 0; i < nbesp; i++) {
        gsl_matrix_set(m, 0, i, 0);
    }
    for (int i = 0; i < nbesp; i++) {
        gsl_matrix_set(m, i, 0, 0);
    }
    for (int i = 0; i < nbesp; i++) {
        dfdt[i] = 0;
    }
    if (behavior == 7) {
        int ngenes = cell.genes.size();
        for (int i = 0; i < nbesp; i++) {
            gsl_matrix_set(m, 2 * ngenes + 1, i, 0);
        }
    }
    if (behavior == 3) {
        int ngenes = cell.genes.size();
        for (int i = 0; i < nbesp; i++) {
            gsl_matrix_set(m, 2 * ngenes + 1, i, 0);
            gsl_matrix_set(m, 2 * ngenes + 2, i, 0);
        }
    }
    return GSL_SUCCESS;
}

int
jacsys(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
    Celleff & cell = *((Celleff *) params);
    const int nbesp = cell.esps.size();
    const int nbtot = (nbesp - 1) / nb_cells;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, nbesp, nbesp);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set_zero(m);
    for (int i = 0; i < nbesp; i++) {
        gsl_matrix_set(m, i, i, -cell.esps[i].cste);
    }
    for (vector<Reacteff>::iterator ireact = cell.reactions.begin(); ireact != cell.reactions.end(); ireact++) {
        const Reacteff & react = *ireact;
        double mval = gsl_matrix_get(m, react.r1, react.r1);
        gsl_matrix_set(m, react.r1, react.r1, mval - react.cste * y[react.r2]);
        mval = gsl_matrix_get(m, react.r1, react.r2);
        gsl_matrix_set(m, react.r1, react.r2, mval - react.cste * y[react.r1]);
        mval = gsl_matrix_get(m, react.r2, react.r1);
        gsl_matrix_set(m, react.r2, react.r1, mval - react.cste * y[react.r2]);
        mval = gsl_matrix_get(m, react.r2, react.r2);
        gsl_matrix_set(m, react.r2, react.r2, mval - react.cste * y[react.r1]);
        mval = gsl_matrix_get(m, react.p1, react.r1);
        gsl_matrix_set(m, react.p1, react.r1, mval + react.cste * y[react.r2]);
        mval = gsl_matrix_get(m, react.p1, react.r2);
        gsl_matrix_set(m, react.p1, react.r2, mval + react.cste * y[react.r1]);
        mval = gsl_matrix_get(m, react.p2, react.r1);
        gsl_matrix_set(m, react.p2, react.r1, mval + react.cste * y[react.r2]);
        mval = gsl_matrix_get(m, react.p2, react.r2);
        gsl_matrix_set(m, react.p2, react.r2, mval + react.cste * y[react.r1]);
    }
    for (vector<Geneeff>::iterator igene = cell.genes.begin(); igene != cell.genes.end(); igene++) {
        for (unsigned int i = 0; i < nb_cells; i++) {
            Geneeff & gene = *igene;
            double sum = 1.0;
            double rate = gene.transrate;
            for (vector<Promoeff>::iterator ipromo = gene.promos.begin(); ipromo != gene.promos.end(); ipromo++) {
                double conc = y[(*ipromo).prot + i * nbtot];
                sum += (*ipromo).eqcte * pow(conc, hill);
                rate += (*ipromo).transrate * (*ipromo).eqcte * pow(conc, hill);
            }
            for (vector<Promoeff>::iterator ipromo = gene.promos.begin(); ipromo != gene.promos.end(); ipromo++) {
                double mval = gsl_matrix_get(m, gene.arn + i * nbtot, (*ipromo).prot + i * nbtot);
                gsl_matrix_set(m, gene.arn + i * nbtot, (*ipromo).prot + i * nbtot, mval + ((*ipromo).eqcte * hill * pow(y[(*ipromo).prot], (hill - 1))) / sum * ((*ipromo).transrate - rate / sum));
            }
        }
    }
    for (int i = 0; i < nbesp; i++) {
        gsl_matrix_set(m, 0, i, 0);
    }
    for (int i = 0; i < nbesp; i++) {
        gsl_matrix_set(m, i, 0, 0);
    }
    for (int i = 0; i < nbesp; i++) {
        dfdt[i] = 0;
    }
    return GSL_SUCCESS;
}

int
jacx(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
    jac(t, y, dfdy, dfdt, params);
    //double tanht=tanh(t-100.0);
    dfdt[2 * nb_genes + 1] = 0;
    //dfdt[7]=2*tanht*(-1.0+gsl_pow_2(tanht));
    return GSL_SUCCESS;
}

int
jacxy(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
    jac(t, y, dfdy, dfdt, params);
    //double tanht=tanh(t-100.0);
    dfdt[2 * nb_genes + 1] = 0;
    dfdt[2 * nb_genes + 2] = 0;
    //dfdt[7]=2*tanht*(-1.0+gsl_pow_2(tanht));
    return GSL_SUCCESS;
}

int
jacx3c(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
    Celleff & cell = *((Celleff *)params);
    int nbtot = (cell.esps.size() - 1) / nb_cells;
    jacsys(t, y, dfdy, dfdt, params);
    dfdt[2 * nb_genes + 1] = 0;
    dfdt[nbtot + 2 * nb_genes + 1] = 0;
    if (nb_cells > 2) {
        dfdt[2 * nbtot + 2 * nb_genes + 1] = 0;
        if (nb_cells > 3) {
            dfdt[3 * nbtot + 2 * nb_genes + 1] = 0;
        }
    }
    cout << "comp jac!!\n";
    return GSL_SUCCESS;
}

double calcscore(const double functomatch[], double result[])
{
    double score = 0;
    for (unsigned int i = 0; i < nb_steps; i++) {
        score += gsl_pow_2(functomatch[i] - result[i]);
    }
    return score;
}

void makeoscill(double functomatch[], double finosc)
{
    double deltat = finosc / (double)nb_steps;
    for (unsigned int i = 0; i < nb_steps; i++) {
        functomatch[i] = 1 + 0.7 * cos(deltat * i * 1.0);
    }
}

void Celleff::recoverstate(double state[])
{
    const int nbesp = esps.size();
    for (int i = 0; i < nbesp; i++) {
        esps[i].qtite = state[i];
    }
}

void Celleff::savestate(double state[])
{
    const int nbesp = esps.size();
    for (int i = 0; i < nbesp; i++) {
        state[i] = esps[i].qtite;
    }
}


double Celleff::scorebiooscill(double finosc, double valstate)
{
    double result[nb_steps];
    double functomatch[nb_steps];
    double score = 0;
    esps[5].qtite = 1.7;
    esps[7].qtite = 0.1;
    //integrfunc(result,finosc);
    makeoscill(functomatch, finosc);
    score += 4 * calcscore(functomatch, result);
    esps[7].qtite = 1.0;
    //integrfunc(result,final_time);
    for (unsigned int i = 0; i < nb_steps; i++) {
        functomatch[i] = valstate;
    }
    score += 1.0 * calcscore(functomatch, result);
    return score;
}

double distcust(double pos1, double pos2)
{
    double mean = (pos1 + pos2) / 2;
    double delta = 0;
    if (mean > 1e-2) delta = mean;
    else delta = 1e-2;
    delta = 2;
    double space = 1.0 - fabs(pos1 - pos2) / delta;
    if (space < 0) return 0;
    else return space;
    //   return 1.0-tanh(fabs(pos1-pos2)/delta)
    //   double space=1.0-fabs(pos1-pos2)/delta;
    if (space < 0) return 0;
    else return space;
}

double distcust(const vd & pos1, const vd & pos2)
{
    double res = 0;
    civd iv2 = pos2.begin() + nb_steps;
    for (civd iv1 = pos1.begin() + nb_steps; iv1 != pos1.end(); iv1++) {
        res += distcust(*iv1, *iv2);
        iv2++;
    }
    return res;
}

int
Celleff::integrgrad(double conc, vd & integrA, vd & integrB, double & score)
{
    initconc(0);
    int ngenes = genes.size();
    esps[2 * ngenes + 1].qtite = conc;
    int ctrl = integrfunc(integrA, integrB, 350.0);
    if (ctrl) return 1;
    esps[2 * ngenes + 1].qtite = 0.1;
    ctrl = integrfunc(integrA, integrB, 350.0);
    if (ctrl) return 1;
    vd dumA, dumB;
    dumA.reserve(2 * nb_steps);
    dumB.reserve(2 * nb_steps);
    initconc(0);
    esps[2 * ngenes + 1].qtite = conc * 1.1;
    ctrl = integrfunc(dumA, dumB, 350.0);
    if (ctrl) return 1;
    esps[2 * ngenes + 1].qtite = 0.1;
    ctrl = integrfunc(dumA, dumB, 350.0);
    if (ctrl) return 1;
    score += 20 * fabs(integrA - dumA);
    score += 20 * fabs(integrB - dumB);
    return 0;
}


double Celleff::scoretristagrad()
{
    int ctrl;
    double score = 0;
    vd integrAhigh, integrAmed, integrAlow, integrBhigh, integrBmed, integrBlow;
    integrAhigh.reserve(2 * nb_steps);
    integrAmed.reserve(2 * nb_steps);
    integrAlow.reserve(2 * nb_steps);
    integrBhigh.reserve(2 * nb_steps);
    integrBmed.reserve(2 * nb_steps);
    integrBlow.reserve(2 * nb_steps);
    ctrl = integrgrad(100.0, integrAhigh, integrBhigh, score);
    if (ctrl) return 1e10;
    ctrl = integrgrad(4.0, integrAmed, integrBmed, score);
    if (ctrl) return 1e10;
    ctrl = integrgrad(0.1, integrAlow, integrBlow, score);
    if (ctrl) return 1e10;
    double d1 = distcust(integrAhigh, integrAmed);
    double d2 = distcust(integrAhigh, integrAlow);
    double d3 = distcust(integrAmed, integrAlow);
    double d4 = distcust(integrBhigh, integrBmed);
    double d5 = distcust(integrBhigh, integrBlow);
    double d6 = distcust(integrBmed, integrBlow);
    if (print_mode) {
        ivd ivAm = integrAmed.begin();
        ivd ivAl = integrAlow.begin();
        ivd ivBh = integrBhigh.begin();
        ivd ivBm = integrBmed.begin();
        ivd ivBl = integrBlow.begin();
        for (ivd iv = integrAhigh.begin(); iv != integrAhigh.end(); iv++) {
            outintegr << *iv << " " <<
                      *ivAm << " " <<
                      *ivAl << " " <<
                      *ivBh << " " <<
                      *ivBm << " " <<
                      *ivBl << "\n";
            ivAm++;
            ivAl++;
            ivBh++;
            ivBm++;
            ivBl++;
        }
    }
    double scorestart = score;
    double lower = mini(d1, d2);
    lower = mini(lower, d3);
    lower = mini(lower, d4);
    lower = mini(lower, d5);
    lower = mini(lower, d6);
    score += lower;
    score += 1000;
    if (score < 1010) {
        score = scorestart + mini(d1, d4) + mini(d2, d5) + mini(d3, d6);
    } else {
        return score;
    }
    return score;
}

double Celleff::scoremultistable()
{
    double result1[nb_steps];
    double result2[nb_steps];
    double result3[nb_steps];
    double result4[nb_steps];
    double result5[nb_steps];
    double result6[nb_steps];
    double result20[nb_steps];
    double result21[nb_steps];
    double interm = 0;
    // double functomatch[nb_steps];
    double score = 0;
    int ctrl = 0;
    initconc(-1);
    int ngenes = genes.size();
    esps[1].qtite = 0.1;
    esps[2].qtite = 0.2;
    esps[ngenes + 1].qtite = 0.1;
    esps[ngenes + 2].qtite = 0.2;
    //   ctrl=integrfunc(result1,result2,350.0);
    if (ctrl) return 1e10;
    initconc(-1);
    esps[1].qtite = 0.12;
    esps[2].qtite = 0.25;
    esps[ngenes + 1].qtite = 0.12;
    esps[ngenes + 2].qtite = 0.25;
    //   ctrl=integrfunc(result20,result21,350.0);
    if (ctrl) return 1e10;
    for (unsigned int i = 0; i < nb_steps; i++) {
        interm += fabs(result1[i] - result20[i]);
        interm += fabs(result2[i] - result21[i]);
    }
    //  score += ;
    initconc(0);
    esps[1].qtite = 10.0;
    esps[2].qtite = 1.5;
    esps[ngenes + 1].qtite = 10.0;
    esps[ngenes + 2].qtite = 1.5;
    //   ctrl=integrfunc(result3,result4,350.0);
    if (ctrl) return 1e10;
    initconc(0);
    esps[1].qtite = 11.0;
    esps[2].qtite = 1.8;
    esps[ngenes + 1].qtite = 11.0;
    esps[ngenes + 2].qtite = 1.8;
    //   ctrl=integrfunc(result20,result21,350.0);
    if (ctrl) return 1e10;
    for (unsigned int i = 0; i < nb_steps; i++) {
        interm += fabs(result3[i] - result20[i]);
        interm += fabs(result4[i] - result21[i]);
    }
    /*
       initconc(0);
       esps[1].qtite=12.0;
       esps[2].qtite=2.5;
       esps[ngenes+1].qtite=12.0;
       esps[ngenes+2].qtite=2.5;
       ctrl=integrfunc(result7,result8,350.0);
       if (ctrl) return 1e10;
       */
    initconc(1);
    esps[1].qtite = 30.0;
    esps[2].qtite = 40.0;
    esps[ngenes + 1].qtite = 30.0;
    esps[ngenes + 2].qtite = 40.0;
    //   ctrl=integrfunc(result5,result6,350.0);
    if (ctrl) return 1e10;
    initconc(1);
    esps[1].qtite = 31.0;
    esps[2].qtite = 41.0;
    esps[ngenes + 1].qtite = 31.0;
    esps[ngenes + 2].qtite = 41.0;
    //   ctrl=integrfunc(result20,result21,350.0);
    if (ctrl) return 1e10;
    for (unsigned int i = 0; i < nb_steps; i++) {
        interm += fabs(result5[i] - result20[i]);
        interm += fabs(result6[i] - result21[i]);
    }
    //  for (int i=0;i<nb_steps;i++){
    //    functomatch[i]=50.0;
    //  }
    //  score += calcscore(functomatch,result);
    double d1, d2, d3, d4, d5, d6;
    d1 = d2 = d3 = d4 = d5 = d6 = 0;
    for (unsigned int i = 0; i < nb_steps; i++) {
        d1 += distcust(result1[i], result3[i]);
        d2 += distcust(result1[i], result5[i]);
        d3 += distcust(result3[i], result5[i]);
        d4 += distcust(result2[i], result4[i]);
        d5 += distcust(result2[i], result6[i]);
        d6 += distcust(result4[i], result6[i]);
    }
    score = mini(d1, d4) + mini(d2, d5) + mini(d3, d6); //+d7+d8+d9+d10;
    score += 20 * interm;
    return score;
}

double Celleff::scoremonostable()
{
    double result1[nb_steps];
    // double functomatch[nb_steps];
    double score = 0;
    int ctrl = 0;
    initconc(-1);
    esps[1].qtite = 0.1;
    esps[2].qtite = 0.2;
    ctrl = integrunifunc(result1, 350.0);
    if (ctrl) return 1e10;
    for (unsigned int i = 20; i < nb_steps; i++) {
        score += fabs(result1[i] - 10.0);
    }
    //  score += ;
    return score;
}

double Celleff::scorebistable()
{
    double result1[nb_steps];
    double result1p[nb_steps];
    //   double result2[nb_steps];
    // double functomatch[nb_steps];
    double score = 0;
    int ctrl = 0;
    initconc(-1);
    int ngenes = genes.size();
    esps[1].qtite = 0.1;
    esps[2].qtite = 15.2;
    esps[ngenes + 1].qtite = 0.1;
    esps[ngenes + 2].qtite = 15.1;
    //   ctrl=integrfunc(result1,result2,500.0);
    if (ctrl) return 1e10;
    for (unsigned int i = 0; i < nb_steps; i++) {
        score += fabs(result1[i] - 0.1);
    }
    //  score += ;
    initconc(1);
    ngenes = genes.size();
    esps[1].qtite = 40.1;
    esps[2].qtite = 0.2;
    esps[ngenes + 1].qtite = 60.0;
    esps[ngenes + 2].qtite = 0.2;
    //   ctrl=integrfunc(result1p,result2,500.0);
    if (ctrl) return 1e10;
    for (unsigned int i = 0; i < nb_steps; i++) {
        score += fabs(result1p[i] - 50.0);
    }
    //  score += ;
    if (score < 1000) {
        score = 0;
        for (unsigned int i = 50; i < nb_steps; i++) {
            score += fabs(result1[i] - 1.0);
        }
        for (unsigned int i = 50; i < nb_steps; i++) {
            score += fabs(result1p[i] - 50.0);
        }
    }
    return score;
}

vd Celleff::scoreoscill()
{
    double score = 0;
    double score_auxi = 0;
    double fintime = 350.0;
    vd scores;
    int ctrl = 0;
    vd integrA, integrB;
    integrA.reserve(2 * nb_steps);
    integrB.reserve(2 * nb_steps);
    initconc(0);
    ctrl = integrfunc(integrA, integrB, fintime);
    if (ctrl) {
        scores.push_back(1e10);
        scores.push_back(0);
        return scores;
    }
    for (unsigned int i = 1; i < nb_steps - 1; i++) {
        score += distcust(integrA[i - 1] + integrA[i + 1], 2 * integrA[i]);
    }
    scores.push_back(score);
    for (unsigned int i = nb_steps / 2; i < nb_steps - 1; i++) {
        score_auxi += distcust(integrA[i - 1] + integrA[i + 1], 2 * integrA[i]);
    }
    scores.push_back(score_auxi);
    if (print_mode) {
        for (unsigned int i = 0; i < nb_steps; i++) {
            outintegr << i + 1 << " " << integrA[i] << "\n";
        }
    }
    return scores;
}

double Celleff::scorebistable2()
{
    double result1[nb_steps];
    double result1p[nb_steps];
    double score = 0;
    int ctrl = 0;
    initconc(-1);
    int ngenes = genes.size();
    esps[1].qtite = 0.1;
    esps[2].qtite = 15.2;
    esps[ngenes + 1].qtite = 0.1;
    esps[ngenes + 2].qtite = 15.1;
    //   ctrl=integrfunc(result1,result2,500.0);
    if (ctrl) return 1e10;
    //  score += ;
    initconc(1);
    ngenes = genes.size();
    esps[1].qtite = 40.1;
    esps[2].qtite = 0.2;
    esps[ngenes + 1].qtite = 60.0;
    esps[ngenes + 2].qtite = 0.2;
    //   ctrl=integrfunc(result1p,result2,500.0);
    if (ctrl) return 1e10;
    //  score += ;
    score = 0;
    for (unsigned int i = 0; i < nb_steps; i++) {
        score += distcust(result1[i], result1p[i]);
    }
    return score;
}

int Celleff::integrunifunc(double result1[], double fintime)
{
    unsigned int espsize = esps.size();
    gsl_odeiv_step * s = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, espsize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new(1e-3, 0.0);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(espsize);
    gsl_odeiv_system sys = {func, jac, espsize, this};
    double deltat = fintime / (double)nb_steps;
    double t = 0.0;
    double t1 = 0.0;
    double h = 1e-3;
    int nbgene = genes.size();
    double * y = new double[espsize];
    for (unsigned int i = 0; i < espsize; i++) {
        y[i] = esps[i].qtite;
    }
    y[0] = 1.0;
    for (unsigned int i = 0; i < nb_steps; i++) {
        t1 += deltat;
        while (t < t1) {
            int status = gsl_odeiv_evolve_apply(e, c, s,
                                                &sys,
                                                &t, t1,
                                                &h, y);
            y[0] = 1.0;
            for (unsigned int j = 0; j < espsize; j++) {
                if (y[j] > 1000 || y[j] < 0.001) {
                    return -1;
                }
            }
            if (status != GSL_SUCCESS) {
                cout << "error while integrating" << endl;
                break;
            }
        }
        result1[i] = y[nbgene + 1];
    }
    for (unsigned int i = 0; i < espsize; i++) {
        esps[i].qtite = y[i];
    }
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    recoverstate(y);
    delete [] y;
    return 0;
}


int Celleff::integrfunc(vd & integrA, vd & integrB, double fintime)
{
    unsigned int espsize = esps.size();
    gsl_odeiv_step * s = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, espsize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new(1e-3, 0.0);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(espsize);
    gsl_odeiv_system sys = {func, jac, espsize, this};
    double deltat = fintime / (double)nb_steps;
    double t = 0.0;
    double t1 = 0.0;
    double h = 1e-3;
    int nbgene = genes.size();
    double * y = new double[espsize];
    for (unsigned int i = 0; i < espsize; i++) {
        y[i] = esps[i].qtite;
    }
    y[0] = 1.0;
    for (unsigned int i = 0; i < nb_steps; i++) {
        t1 += deltat;
        while (t < t1) {
            int status = gsl_odeiv_evolve_apply(e, c, s,
                                                &sys,
                                                &t, t1,
                                                &h, y);
            y[0] = 1.0;
            for (unsigned int j = 0; j < espsize; j++) {
                if (y[j] > 500 || y[j] < 0.002) {
                    gsl_odeiv_evolve_free(e);
                    gsl_odeiv_control_free(c);
                    gsl_odeiv_step_free(s);
                    delete [] y;
                    return -1;
                }
            }
            if (status != GSL_SUCCESS) {
                cout << "error while integrating" << endl;
                break;
            }
        }
        //      cout << y[2*nbgene+1] << endl;
        integrA.push_back(y[nbgene + 1]);
        if (behavior != 1) {
            integrB.push_back(y[nbgene + 2]);
        }
    }
    for (unsigned int i = 0; i < espsize; i++) {
        esps[i].qtite = y[i];
    }
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    recoverstate(y);
    delete [] y;
    return 0;
}


vector<double> Celleff::integrporte(int nprotA, int nprotX, int nprotY, bool protX1, bool protY1, bool protX2, bool protY2)
{
    vector<double> resultat;
    return resultat;
}

vector<double> Celleff::integrbistable(int nprot)
{
    vector<double> resultat;
    return resultat;
}

vector< vector<double> > Celleff::integrmultistable(int prot0, int prot1, int prot2)
{
    vector< vector<double> > resultat;
    return resultat;
}

Especeeff::Especeeff(double c, double q)
    : cste(c), qtite(q)
{
    dirinit = 0;
}

Reacteff::Reacteff(double c, unsigned int re1,
                   unsigned int re2, unsigned int pr1, unsigned int pr2)
    : cste(c), r1(re1), r2(re2), p1(pr1), p2(pr2)
{
}

Geneeff::Geneeff(unsigned int narn, double dtransrate, Gene * gn, Cellule & cell)
    : arn(narn), transrate(dtransrate)
{
    for (vector<Promoter *>::iterator ipromo = cell.promoters.begin(); ipromo != cell.promoters.end(); ipromo++) {
        Promoter & promo = *(*ipromo);
        if (promo.gene == gn) {
            promos.push_back(Promoeff(cell.especenum(promo.protein), promo.transrate, promo.eqcte));
        }
    }
}

Promoeff::Promoeff(unsigned int p, double t, double e)
    : eqcte(e), transrate(t), prot(p)
{
}

Celleff::Celleff(Cellule & cell)
{
    vector<Gene *>::iterator igene;
    vector<Arn *>::iterator iarn;
    vector<Protein *>::iterator iprot;
    vector<Promoter *>::iterator ipromo;
    vector<Reaction *>::iterator ireact;
    Especeeff dummy = Especeeff(0, 1);
    esps.push_back(dummy);
    for (igene = cell.genes.begin(); igene != cell.genes.end(); igene++) {
        Gene & gen = *(*igene);
        Geneeff gn = Geneeff(cell.especenum(gen.arn), gen.transrate, *igene, cell);
        genes.push_back(gn);
    }
    if (print_mode) {
        cout << genes.size() << "\n";
    }
    for (iarn = cell.arns.begin(); iarn != cell.arns.end(); iarn++) {
        Especeeff esp = Especeeff((*iarn)->cste, (*iarn)->qtite);
        esp.dirinit = 1;
        esp.meanc = (*iarn)->qtite;
        esps.push_back(esp);
        if (print_mode) {
            cout << (*iarn)->label << "\n";
        }
    }
    for (iprot = cell.proteins.begin(); iprot != cell.proteins.end(); iprot++) {
        Especeeff esp = Especeeff((*iprot)->cste, (*iprot)->qtite);
        esp.meanc = (*iprot)->qtite;
        esp.dirinit = (*iprot)->dirinit;
        esps.push_back(esp);
        if (print_mode) {
            cout << (*iprot)->label << "\n";
        }
    }
    for (ireact = cell.reactions.begin(); ireact != cell.reactions.end(); ireact++) {
        double c = (*ireact)->cste;
        unsigned int nr1 = cell.especenum((*ireact)->reactif1);
        unsigned int nr2 = cell.especenum((*ireact)->reactif2);
        unsigned int np1 = cell.especenum((*ireact)->produit1);
        unsigned int np2 = cell.especenum((*ireact)->produit2);
        reactions.push_back(Reacteff(c, nr1, nr2, np1, np2));
    }
    initconc(0);
    //double result[nb_steps];
    //integrfunc(result,result,200.0);
    //printlong();
}


void Celleff::printlong()
{
    for (ivgeff igene = genes.begin(); igene != genes.end(); igene++) {
        cout << "gene" << endl;
    }
    for (ivreaceff ireact = reactions.begin(); ireact != reactions.end(); ireact++) {
        cout << "react " << (*ireact).r1 << " " << (*ireact).r2 << " " << (*ireact).p1 << " " << (*ireact).p2 << endl;
    }
    for (ivespeff iesp = esps.begin(); iesp != esps.end(); iesp++) {
        cout << "esp\n";
    }
}


void Celleff::initconc(int state)
{
    for (ivespeff iesp = esps.begin(); iesp != esps.end(); iesp++) {
        Especeeff & esp = *iesp;
        //      if (esp.dirinit!=0){
        esp.qtite = esp.meanc * pow(2.0, state * esp.dirinit);
        //      }
    }
}

void Celleff::initintegr(Integrsys & intsys)
{
    intsys.step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, esps.size());
    intsys.control = gsl_odeiv_control_y_new(1e-6, 1e-3);
    intsys.evolve = gsl_odeiv_evolve_alloc(esps.size());
}

void freeintegr(Integrsys & intsys, double * & y)
{
    gsl_odeiv_evolve_free(intsys.evolve);
    gsl_odeiv_control_free(intsys.control);
    gsl_odeiv_step_free(intsys.step);
    delete [] y;
}

void Celleff::initstarty(double * & y)
{
    unsigned int espsize = esps.size();
    for (unsigned int i = 0; i < espsize; i++) {
        y[i] = esps[i].qtite;
    }
    y[0] = 1.0;
}

void inityA(double * & y, unsigned int & ysize, double cell1, double cell2, double cell3, double cell4)
{
    int nbtot = (ysize - 1) / nb_cells;
    if (anchorsignal) {
        y[2 * nb_genes + 1] = cell1;
        y[nbtot + 2 * nb_genes + 1] = cell2;
        if (nb_cells > 2) {
            y[2 * nbtot + 2 * nb_genes + 1] = cell3;
            if (nb_cells > 3) {
                y[3 * nbtot + 2 * nb_genes + 1] = cell4;
            }
        }
    } else {
        y[nb_genes + 1] = cell1;
        y[nbtot + nb_genes + 1] = cell2;
        if (nb_cells > 2) {
            y[2 * nbtot + nb_genes + 1] = cell3;
            if (nb_cells > 3) {
                y[3 * nbtot + nb_genes + 1] = cell4;
            }
        }
    }
}

void respush(double * & y, unsigned int & ysize, Resultatsys & res)
{
    const int nbtot = (ysize - 1) / nb_cells;
    res.result1a.push_back(y[nb_genes + 1]);
    res.result1b.push_back(y[nb_genes + 2]);
    res.result2a.push_back(y[nbtot + nb_genes + 1]);
    res.result2b.push_back(y[nbtot + nb_genes + 2]);
    if (nb_cells > 2) {
        res.result3a.push_back(y[2 * nbtot + nb_genes + 1]);
        res.result3b.push_back(y[2 * nbtot + nb_genes + 2]);
        if (nb_cells > 3) {
            res.result4a.push_back(y[3 * nbtot + nb_genes + 1]);
            res.result4b.push_back(y[3 * nbtot + nb_genes + 2]);
        }
    }
}

void printres(double & time, double * & y, unsigned int & ysize)
{
    const int nbtot = (ysize - 1) / nb_cells;
    if (nb_cells < 3) {
        outintegr << time << " " << y[nb_genes + 1] << " " << y[nbtot + nb_genes + 1] << " " << y[nb_genes + 2] << " " << y[nbtot + nb_genes + 2] << endl;
    } else {
        outintegr << time << " " << y[nb_genes + 1] << " " << y[nbtot + nb_genes + 1] << " " << y[2 * nbtot + nb_genes + 1] << " " << y[3 * nbtot + nb_genes + 1] << " " << y[nb_genes + 2] << " " << y[nbtot + nb_genes + 2] << " " << y[2 * nbtot + nb_genes + 2] << " " << y[3 * nbtot + nb_genes + 2] << endl;
    }
}

int integrcells(Integrsys & intsys, double & t, double & t1, double * & y, unsigned int & ysize)
{
    static double h = 1e-3;
    while (t < t1) {
        int status = gsl_odeiv_evolve_apply(intsys.evolve, intsys.control, intsys.step,
                                            &intsys.system,
                                            &t, t1,
                                            &h, y);
        y[0] = 1.0;
        unsigned int nb_var = (ysize - 1) / nb_cells;
        for (unsigned int j = 0; j < ysize; j++) {
            if (isnan(y[j])) {
                freeintegr(intsys, y);
                return -1;
            }
            //         if (y[j]>500 || y[j]<0.002){
            if ((j % nb_var) != ((2 * nb_genes + 1) % nb_var) && (y[j] > 100 || y[j] < 0.01)) {
                //         if (y[j]>500 || isnan(y[j])){
                freeintegr(intsys, y);
                return -1;
            }
            if (y[j] < 1e-8) y[j] = 1e-8;
        }
        if (status != GSL_SUCCESS) {
            cout << "error while integrating" << endl;
            break;
        }
    }
    return 0;
}

int Celleff::integrsys3c(System & sys, Resultatsys & res, vd initc)
{
    unsigned int espsize = esps.size();
    Integrsys intsys;
    initintegr(intsys);
    intsys.system = (gsl_odeiv_system) {
        funcx3c, jacx3c, espsize, this
    };
    double deltat = final_time / (double)nb_steps;
    double t = 0.0;
    double t1 = 0.0;
    double * y = new double[espsize];
    initstarty(y);
    for (unsigned int i = 0; i < nb_steps; i++) {
        t1 += deltat;
        if (recept_test) {
            double mean_init = (initc[0] + initc[1]) / 2.0;
            inityA(y, espsize, mean_init * 1.05, mean_init * 0.95, 0.98, 0.96);
        } else {
            if (nb_cells > 2) {
                //            inityA(y, espsize, 4.0,1.0,0.25,0.25);
                inityA(y, espsize, 10.0, 1.0, 0.1, 0.1);
                //            inityA(y, espsize, 20.0,20.0,0.05,0.05);
            } else {
                inityA(y, espsize, initc[0], initc[1], 0.05, 0.05);
            }
        }
        if (integrcells(intsys, t, t1, y, espsize)) return -1;
        respush(y, espsize, res);
        if (print_mode) printres(t1, y, espsize);
    }
    deltat = 2 * final_time / (double)nb_steps;
    if (print_mode) {
        deltat *= 10;
        cout << initc[0] << " " << initc[1] << endl;
        cout << recept_test << endl;
    }
    for (unsigned int i = 0; i < nb_steps; i++) {
        t1 += deltat;
        if (anchorsignal) inityA(y, espsize, 1e-4, 1e-4, 1e-4, 1e-4);
        if (integrcells(intsys, t, t1, y, espsize)) return -1;
        respush(y, espsize, res);
        if (print_mode) printres(t1, y, espsize);
    }
    // Re - integration to inspect stability of solutions
    if (nb_cells > 2) {
        initstarty(y);
        deltat = final_time / (double)nb_steps;
        for (unsigned int i = 0; i < nb_steps; i++) {
            t1 += deltat;
            if (nb_cells > 2) {
                inityA(y, espsize, 10.0 * 1.1, 1.0 * 1.1, 0.1 * 1.1, 0.1 * 1.1);
                //            inityA(y, espsize, 4.0*1.1,1.0*1.1,0.25*1.1,0.25*1.1);
                //         inityA(y, espsize, 22.0,22.0,0.055,0.055);
            } else {
                inityA(y, espsize, initc[0] * 1.1, initc[0] * 1.1, 0.1, 0.1);
            }
            if (integrcells(intsys, t, t1, y, espsize)) return -1;
            respush(y, espsize, res);
            if (print_mode) printres(t1, y, espsize);
        }
        deltat = 2 * final_time / (double)nb_steps;
        for (unsigned int i = 0; i < nb_steps; i++) {
            t1 += deltat;
            if (anchorsignal) inityA(y, espsize, 1e-4, 1e-4, 1e-4, 1e-4);
            if (integrcells(intsys, t, t1, y, espsize)) return -1;
            respush(y, espsize, res);
            if (print_mode) printres(t1, y, espsize);
        }
    }
    //   recoverstate(y);
    freeintegr(intsys, y);
    return 0;
}

double
scoreres(vector<double> vectdb, double value)
{
    int i = 0;
    double s = 0;
    vector<double>::iterator idb;
    for (idb = vectdb.begin(); idb != vectdb.end(); idb++) {
        s += (*idb - value) * (*idb - value) * 2 * i / nb_steps;
        i++;
    }
    return s;
}

int
Celleff::integrgate(vd & res)
{
    res.clear();
    unsigned int espsize = esps.size();
    double fintime = 500.0;
    gsl_odeiv_step * s = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, espsize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new(1e-3, 0.0);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc(espsize);
    gsl_odeiv_system sys = {func, jac, espsize, this};
    double deltat = fintime / (double)nb_steps;
    double t = 0.0;
    double t1 = 0.0;
    double h = 1e-3;
    int nbgene = genes.size();
    double * y = new double[espsize];
    for (unsigned int i = 0; i < espsize; i++) {
        y[i] = esps[i].qtite;
    }
    y[0] = 1.0;
    for (unsigned int i = 0; i < nb_steps; i++) {
        t1 += deltat;
        while (t < t1) {
            int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, y);
            y[0] = 1.0;
            for (unsigned int j = 0; j < espsize; j++) {
                if (y[j] > 200 || y[j] < 0.005) {
                    gsl_odeiv_evolve_free(e);
                    gsl_odeiv_control_free(c);
                    gsl_odeiv_step_free(s);
                    delete [] y;
                    return -1;
                }
            }
            if (status != GSL_SUCCESS) {
                cout << "error while integrating" << endl;
                break;
            }
        }
        res.push_back(y[nbgene + 1]);
        if (print_mode) {
            outintegr << t1 << " " << y[nbgene + 1] << "\n";
        }
    }
    outintegr << "\n";
    for (unsigned int i = 0; i < espsize; i++) {
        esps[i].qtite = y[i];
    }
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    recoverstate(y);
    delete [] y;
    return 0;
}

double
Celleff::scoregate(bool i, bool j)
{
    int ngenes = genes.size();
    if (i) {
        esps[2 * ngenes + 1].qtite = 1.0;
    } else {
        esps[2 * ngenes + 1].qtite = 0.01;
    }
    if (j) {
        esps[2 * ngenes + 2].qtite = 1.0;
    } else {
        esps[2 * ngenes + 2].qtite = 0.01;
    }
    vd res;
    int ctrl = integrgate(res);
    if (ctrl == -1) {
        return 1e6;
    }
    double s = 0;
    if (i ^ j) {
        for (ivd ires = res.begin(); ires != res.end(); ires++) {
            s += distcust(0, *ires);
        }
        if (print_mode) {
            cout << i << " " << j << " " << (i ^ j) << "\n";
            cout << "X " << esps[2 * ngenes + 1].qtite << " Y " << esps[2 * ngenes + 2].qtite << "\n";
        }
    } else {
        for (ivd ires = res.begin(); ires != res.end(); ires++) {
            s += *ires;
        }
        if (print_mode) {
            cout << i << " " << j << " " << (i ^ j) << "\n";
            cout << "X " << esps[2 * ngenes + 1].qtite << " Y " << esps[2 * ngenes + 2].qtite << "\n";
        }
    }
    return s;
}

double
Celleff::scorelogicgate()
{
    double s = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            initconc(0);
            s += scoregate(i, j);
        }
    }
    return s;
}

void
Celleff::down_recept()
{
    //   unsigned int nb_react_recept=sys.recepteurs.size()*2*(nb_cells-1);
    //   for (ivreac ir=reactions.end()-nb_react_recept;ir!=reactions.end();ir++){
    //      ir->cste=0;
    //   }
    for (ivrecepte ire = recepts.begin(); ire != recepts.end(); ire++) {
        ire->kphospho = 0.0;
    }
}

Recepteff::Recepteff(unsigned int ps, unsigned int pc, unsigned int pp , double k, double e)
    : eqcte(e), kphospho(k), prots(ps), protc(pc), protp(pp)
{
}
