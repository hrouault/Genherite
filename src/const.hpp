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

#ifndef Const_H
#define Const_H

// Contraintes générales du système
extern unsigned int nb_promoters_max;
extern unsigned int nb_proteins_max;
const unsigned int nb_reactions_max = 80;
const unsigned int nb_recepteurs_max = 3;
//const unsigned int nb_recepteurs_max = 8;

// Taux relatifs de modifications aléatoires

const double t_new_r = 1.2;     //0.1
const double t_mod_r = 2.0;
const double t_mod_pro = 2.0;
const double t_mod_qpromo = 2.0;
const double t_mod_qtite = 3.0;
const double t_mod_recept = 2.0;
const double t_phospho = 1.2;   //0.1
const double t_new_promo = 0.3; //0.6;//0.05
const double t_rmreact = 0.2;   //0.4;//0.05
const double t_mutdir = 1.8;

extern double t_new_recept;     //0.4;

// Taux relatifs des modifications aléatoires lors de l'optimisation
const double t_mod_r_opti = 3.0;
const double t_mod_pro_opti = 3.0;
const double t_mod_qtite_opti = 2.0;
const double t_rmreact_opti = 0.2;
const double t_mod_recept_opti = 2.0;

//Amplitude de la modification des quantités lors de l'optimisation
const double ampli_opti = 1.2;

// Seuils des constantes, qtites, etc
const double seuil_degrad_inf = 0.1;
const double seuil_degrad_sup = 20.0;
const double seuil_react_inf = 0.1;
const double seuil_react_sup = 20.0;
const double seuil_qtite_inf = 0.10;
const double seuil_qtite_sup = 20.0;

//Nombre de pas d'intégration
extern unsigned int nb_steps;
const int nb_steps_init = 40000;        // Nb de pas d'intégration initiaux dans la procédure porte
const double final_time = 150.0;

//Nombre de générations
const int nb_generations_opti = 300;

// Valeurs du comportement bistable
const double val_haute = 1;
const double val_basse = 0;
const double val_init = 0.1;
const double val_switch = 1.6;

// Constantes des schémas spatiaux
extern unsigned int nb_cells;

//Sore limite avant optimisation
const double score_limite = 1;

//Pulsation et amplitude de l'oscillation
const double omega = 0.3;
const double am = 0.7;

//Pénalisations du nb de réactions et nb especes, etc
const double coef_react = 10;
const double coef_prot = 10;
const double coef_promo = 10;
const double coef_recept = 10;

//amplitude des modifications
const double defampli = 2;

//valeur des portes (deux protéines)
const bool val_0_0 = false;
const bool val_a_a = false;
const bool val_a_b = false;
const bool val_a_ab = false;
const bool val_b_0 = false;
const bool val_b_a = false;
const bool val_b_b = false;
const bool val_b_ab = false;
const bool val_ab_0 = false;
const bool val_ab_a = false;
const bool val_ab_b = false;
const bool val_ab_ab = false;

// Méthode d'intégration
// 0 pour Euler
// 1 pour Runge-Kutta
const int methode_int = 0;

// Nb de réseaux à générer

// Utilisation d'un réseau préexistant
const bool use_net = false;     //true;

// Exposant de Hill pour les promoteurs
const unsigned int hill = 2;

#endif                          /* Const_H */
