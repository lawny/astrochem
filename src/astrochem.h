/*
   astrochem.h - astrochem Function prototypes, various constant and data
   structures for Astrochem.

   Copyright (c) 2006-2013 Sebastien Maret

   This file is part of Astrochem.

   Astrochem is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Astrochem is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Astrochem.  If not, see <http://www.gnu.org/licenses/>.
   */

#/**
  * @file astrochem.h
  * @author Sebastion Maret
  * @date 28 August 2014
  * @brief File containing astrochem related private declaration
  */


/* Various definitions and constants */

#ifndef _ASTROCHEM_H_
#define _ASTROCHEM_H_

#include "libastrochem.h"

#ifdef WITH_MPI

#ifdef HAVE_OPENMP
#define MPI_GRAIN_DEFAULT OMP_NUM_THREADS
#else
#define MPI_GRAIN_DEFAULT 10
#endif

#define MPI_MASTER_RANK 0
#define MPI_FIRST_WORKER 1

typedef struct
{
  cell_table_t *cells;
  int *cell_idxs;
} cell_block_t;

#endif

#endif // _ASTROCHEM_H_
