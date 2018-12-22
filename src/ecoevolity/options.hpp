/******************************************************************************
 * Copyright (C) 2015-2016 Jamie R. Oaks.
 *
 * This file is part of Ecoevolity.
 *
 * Ecoevolity is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Ecoevolity is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Ecoevolity.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef ECOEVOLITY_OPTIONS_HPP
#define ECOEVOLITY_OPTIONS_HPP


class EcoevolityOptions {

    public:
        EcoevolityOptions() { }

		enum ModelPrior {
            fixed = 0,
            uniform = 1,
            dpp = 2,
            pyp = 3,
        };

        enum ModelOperator {
            none = 0,
            rj = 1,
            gibbs_dpp = 2,
            gibbs_pyp = 3,
        };

};

#endif
