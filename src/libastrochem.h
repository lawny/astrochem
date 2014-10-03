/*
   astrochem.h - Function prototypes, various constant and data
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

/**
 * @mainpage Astrochem
 *
 * Astrochem is a code to study the chemistry of a variety of astromical objects (dense clouds, prestellar cores, protoplanetary disks, etc.).
 * It is provided in the form of a library ( libastrochem ) and a standalone program ( astrochem ).
 *
 * The program astrochem is made to be used with an input file : 
 *<a href="http://smaret.github.io/astrochem/">http://smaret.github.io/astrochem/</a> 
 *
 * The library libastrochem ( documentation : @link libastrochem.h @endlink ) provide fonctions to easily create a system and solve it.
 * 
 */

/**
 * @file libastrochem.h
 * @author Sebastion Maret
 * @date 28 August 2014
 * @brief libastrochem public api
 * @note This is the public api of libastrochem, One can use it to create a solvable system, solve it and recover output abundances. <br>
 * An example can be found here : @link apiuser.c @endlink <br>
 * How to use it : <br>
 * @code
  int verbose = 1;
  char *chem_file = "../networks/osu2009.chm";
  net_t network;
  read_network(chem_file, &network, verbose );

  phys_t phys;
  phys.cosmic = 1e-17;
  phys.chi = 0;
  phys.grain_size = GRAIN_SIZE_DEFAULT;
  phys.grain_abundance = 0;

  double abs_err, rel_err;
  abs_err = ABS_ERR_DEFAULT;
  rel_err = REL_ERR_DEFAULT;

  const char* species[]  = {"CO", "HCO(+)", "e(-)"};
  const double initial_abundances[] = {1e-4, 1e-9, 1e-9};

  double *abundances;
  alloc_abundances( &network, &abundances ); // Allocate the abundances array; it contains all species.
  set_initial_abundances(species, 3, initial_abundances, &network, abundances); // Set initial abundances

  double density = 1000;
  double av = 20;
  double temperature = 10;

  cell_t cell;
  cell.nh = &density;
  cell.av = &av;
  cell.tgas = &temperature;
  cell.tdust = &temperature; // Assume tgas = tdust in this specific case

  astrochem_mem_t astrochem_mem;

  if( solver_init( &cell, &network, &phys, abundances , density, abs_err, rel_err, &astrochem_mem ) != 0 )
    {
      return EXIT_FAILURE;
    }
  int i;
  double time = 0;
  for( i = 0; i< 1000000 ; i++)
    {
      time += 1e-6; // advance time
      solve( &astrochem_mem, &network, abundances, time, verbose);

      // Do something with the results of abundances computations 
    }
  solver_close( &astrochem_mem );
  free_abundances( abundances );
  free_network (&network);
 @endcode
 * @todo move define in specific headers TODO
 * @todo last for define are not documented TODO
 */


/* Various definitions and constants */
#ifndef _LIBASTROCHEM_H_
#define _LIBASTROCHEM_H_

#include <nvector/nvector_serial.h>

#define MAX_LINE 512      /*!< Maximum number of characters in each input file line */
#define CHI_DEFAULT 1     /*!< Default chi value */
#define COSMIC_DEFAULT 1.3e-17  /*!< Default cosmic value */
#define GRAIN_SIZE_DEFAULT 1e-5 /*!< Default Grain radius, in cm */
#define TI_DEFAULT 1e-6  /*!< Default initial time */
#define TF_DEFAULT 1e7   /*!< Default final time */
#define ABS_ERR_DEFAULT 1e-20   /*!< Default absolute error */
#define REL_ERR_DEFAULT 1e-3    /*!< Default relative error */
#define TIME_STEPS_DEFAULT 32   /*!< Default number of times steps */
#define TRACE_ROUTES_DEFAULT 0  /*!< Deactivate route tracing by default */
#define N_OUTPUT_ROUTES 16      /*!< Defaults number of output routes */

#ifndef M_PI
#define M_PI  3.14159265358979323846264338327950288 /*!< Our own value of pi */
#endif

#define MAX_CHAR_SPECIES 32     /*!< Maximum number of characters in a specie name */

#define CONST_MKSA_YEAR 3.1536e7                /*!< Number of seconds in a year */
#define CONST_CGSM_BOLTZMANN (1.3806503e-16)    /*!< Boltzmann constant */
#define CONST_CGSM_MASS_PROTON (1.67262158e-24) /*!< Proton Mass */

#define MIN_ABUNDANCE 1e-20     /*!< Minimum abundance to write in output files */

#define FRACTION_TIME_GRAIN_70K 3.16e-19 
#define GAS_DUST_NUMBER_RATIO 7.57e+11    
#define GRAIN_SITES_PER_CM2 3.00e+15    /*!< cm-2 */
#define AVERAGE_UV_IRSF 1e8     /*!< photons cm-2 */

/* Data structures */

/**
 * @brief type of source
 * A Source can be static or dynamic.
 * With dynamic source, parameters can changed over time
 */
typedef enum
{ 
  STATIC = 0,  /*!< Static source */
  DYNAMIC = 1  /*!< Dynamic source */
} SOURCE_MODE;

/**
 * @brief type for species name
 */
typedef char species_name_t[MAX_CHAR_SPECIES];

/**
 * @brief struct containing an abundances and it's related index species in the network
 */
typedef struct
{
  int species_idx; /*!< Index of specie in network */
  double abundance; /*!< Abundance of specie */
} abund_t;

/**
 * @brief struct containing file names of chem file and source file
 */
typedef struct
{
  char chem_file[MAX_LINE]; /*!< Path to chem file */
  char source_file[MAX_LINE]; /*!< Path to source file */
} files_t;

/**
 * @brief struct containing physics parameters
 */
typedef struct
{
  double chi;             /*!< chi */
  double cosmic;          /*!< cosmic */
  double grain_size;      /*!< grain size */
  double grain_abundance; /*!< grain abundances */
} phys_t;

/**
 * @brief struct containing solver parameters
 */
typedef struct
{
  double ti;       /*!< initial time */
  double tf;       /*!< final time */
  double abs_err;  /*!< absolute error */
  double rel_err;  /*!< relative error */
} solver_t;

/**
 * @brief struct containing array of abundances 
 */
typedef struct
{
  abund_t *initial_abundances; /*!< Array of abundances */
  int n_initial_abundances;    /*!< Number of abundances in array */
} abundances_t;

/**
 * @brief struct containing output parameters
 */
typedef struct
{
  int *output_species_idx; /*!< array of output species idx */
  int n_output_species;    /*!< number of output species */
  int time_steps;          /*!< time steps used */
  int trace_routes;        /*!< If routes have been traced */
  char suffix[MAX_LINE];   /*!< Suffix to sue in output files */
} output_t;

/**
 * @brief struct containing input parametrs
 */
typedef struct
{
  files_t files;    /*!< Input files */
  phys_t phys;      /*!< Physics parameters */
  solver_t solver;  /*!< Solver parameters */
  abundances_t abundances; /*!< abundances */
  output_t output; /*!< Output parameters */
} inp_t;

/**
 * @brief struct containing cell parameters 
 */
typedef struct
{
  double *av;    /*!< av */
  double *nh;    /*!< density */
  double *tgas;  /*!< gas temperature */
  double *tdust; /*!< dust temperature */
} cell_table_t;

/**
 * @brief struct containing cell parameters
 */
typedef struct
{
  double av;    /*!< av */
  double nh;    /*!< density */
  double tgas;  /*!< gas temperature */
  double tdust; /*!< dust temperature */
} cell_t;


/**
 * @brief struct containing array of time steps
 */
typedef struct
{
  double *time_steps; /*!< Time steps */
  int n_time_steps;   /*!< Number of time steps */
} time_steps_t;

/**
 * @brief struct containing a source model 
 */
typedef struct
{
  cell_table_t *cell;    /*!< Array of cells */
  time_steps_t ts; /*!< Time steps */
  int n_cells;     /*!< Number of cells */
  SOURCE_MODE mode; /*!< Source mode */
} mdl_t;

/**
 * @brief struct containing a reaction
 */
typedef struct
{
  int reactant1; /*!< reactant 1*/
  int reactant2; /*!< reactant 2*/
  int reactant3; /*!< reactant 3*/
  int product1;  /*!< product 1*/
  int product2;  /*!< product 2*/
  int product3;  /*!< product 3*/
  int product4;  /*!< product 4*/
  double alpha;  /*!< reaction alpha*/
  double beta;   /*!< reaction beta*/
  double gamma;  /*!< reaction gamma*/
  int reaction_type; /*!< reaction type*/
  int reaction_no;   /*!< reaction number*/
} react_t;

/**
 * @brief struct containing a network
 */
typedef struct
{
  int n_species; /*!< number of species */
  int n_alloc_species; /*!< number of actully allocated species */
  species_name_t *species_names; /*!< species name */
  int n_reactions; /*!< number of reactions */
  react_t *reactions; /*!< array of reactions */
} net_t;

/**
 * @brief struct containing a rates of a reaction 
 */
typedef struct
{
  int reaction_no; /*!< number of concerned reaction */
  double rate;     /*!< rate value */
} r_t;

/**
 * @brief struct containing a route
 */
typedef struct
{
  r_t destruction; /*!< rate of destruction */
  r_t formation;   /*!< rate of formation */
} rout_t;

/**
 * @brief struct containing a result 
 */
typedef struct
{
  double *abundances; /*!< abundance result */
  rout_t *routes;     /*!< routes */
  int n_cells;        /*!< number of cells */
  int n_time_steps;   /*!< number of time steps */
  int n_output_abundances; /*!< number of output abundances */
} res_t;

/**
 * @brief bool enum
 */
typedef enum { false, true } bool;

/**
 * @brief solver specific struct
 */
typedef struct
{
  double *reac_rates; /*!< reaction rate */
  const react_t *reactions; /*!< reactions */
  int n_reactions;   /*!< number of reactions */
  int n_species;     /*!< number of species */
  double nh;         /*!< density */
  double av;         /*!< av */
  double tgas;       /*!< gas temperature */
  double tdust;      /*!< dust temparature */
  double chi;        /*!< chi */
  double cosmic;     /*!< cosmis */
  double grain_size; /*!< grain size */
  double grain_abundance; /*!< grain abundance */
} params_t;

/**
 * @brief struct containing solver related memory
 */
typedef struct
{
  void* cvode_mem; /*!< CVODE specific pointer */
  N_Vector y;      /*!< output vector */
  params_t params; /*!< user solver params */
  double density;  /*!< density */
} astrochem_mem_t; 


/* Fonction prototypes */
int alloc_abundances( const net_t* network, double** abundances );

void free_abundances( double* abundances );

int set_initial_abundances( const char** species, int n_initialized_abundances,
                            const double* initial_abundances, const net_t* network, double* abundances );

int solver_init( const cell_t* cell, const net_t* network, const phys_t* phys,
                 const double* abundances , double density, double abs_err, double rel_err,
                 astrochem_mem_t* astrochem_mem );

int solve( astrochem_mem_t* astrochem_mem, const net_t* network,
           double* abundances, double time , const cell_t* new_cell, int verbose );

void solver_close( astrochem_mem_t* astrochem_mem );


void read_input (const char *input_file, inp_t * input_params,
                 const net_t * network, int verbose);

void read_input_file_names (const char *input_file, files_t * files,
                            int verbose);

void free_input (inp_t * input_params);

void read_source (const char *source_file, mdl_t * source_mdl,
                  const inp_t * input_params, const int verbose);

void free_mdl (mdl_t * source_mdl);

void read_network (const char *chem_file, net_t * network, const int verbose);

void free_network (net_t * network);

void alloc_results (res_t * results, int n_time_steps, int n_cells,
                    int n_output_abundances);

void free_results (res_t * results);

int full_solve (int cell_index, const inp_t * input_params, SOURCE_MODE mode,
                const cell_table_t * cell, const net_t * network,
                const time_steps_t * ts, res_t * results, int verbose);

int get_abundance_idx (const res_t * results, int cell_idx, int ts_idx,
                       int abund_idx);

int get_route_idx (const res_t * results, int cell_idx, int ts_idx,
                   int abund_idx, int route_idx);

void alloc_input (inp_t * input_params, int n_initial_abundances,
                  int n_output_abundances);

void alloc_network (net_t * network, int n_species, int n_reactions);

void output (int n_cells, const inp_t * input_params,
             const mdl_t * source_mdl, const net_t * network,
             const res_t * results, int verbose);


#endif // _LIBASTROCHEM_H_
