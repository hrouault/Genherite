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

#ifndef _Evolution_H_
#define _Evolution_H_

#include "cellule.hpp"

class Milieu {
  public:
    vector < Cellule * >cellules;

    Milieu();
    ~Milieu();
    void sort();
    void selection();
    void selection_temp(double temperature);
    void evolution();
    void optimisation();
    void optiselection();
};

bool compcellules(const Cellule * cell1, const Cellule * cell2);

#endif                          /* _Evolution_hpp_ */
