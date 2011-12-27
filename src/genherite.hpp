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

#ifndef _GenHerite_H_
#define _GenHerite_H_

#include <gsl/gsl_rng.h>
#include <fstream>
#include <vector>

using namespace std;

typedef vector<double> vd;
typedef vd::iterator ivd;
typedef vd::const_iterator civd;

vd operator+(const vd & vec1, const vd & vec2);
vd operator-(const vd & vec1, const vd & vec2);

double fabs(const vd & vec);

double frand();
double frand2();

extern struct gengetopt_args_info args_info;
extern unsigned int pop_size;
extern unsigned int nb_genes;

extern bool print_mode;

extern unsigned int behavior;
extern unsigned int recept_satur;
extern unsigned int recept_test;
extern unsigned int anchorsignal;

extern unsigned int opti_stage;

extern std::ofstream outintegr;

#endif /* _GenHerite_H_ */
