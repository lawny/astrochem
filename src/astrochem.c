/* 
   Astrochem - compute the abundances of chemical species in the
   interstellar medium as as function of time.

   Copyright (c) 2006-2013 Sebastien Maret

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <mpi.h>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "astrochem.h"

void usage (void);
void version (void);

int
main (int argc, char *argv[])
{
  int cell_index;
#ifdef WITH_MPI
  int numprocs, rank, namelen;
  int mpi_grain = MPI_GRAIN_DEFAULT;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name (processor_name, &namelen);
  int nb_worker = numprocs - MPI_FIRST_WORKER;
  if (nb_worker < 1)
    {
      fprintf (stderr,
               "astrochem: %s:%d: There must be at least one worker, use mpirun\n",
               __FILE__, __LINE__);
      exit (1);
    }
  int tag = 1;
  MPI_Status status;

  if (rank == MPI_MASTER_RANK)
    {
#endif
      inp_t input_params;
      mdl_t source_mdl;
      net_t network;
      res_t results;

      int verbose = 1;
      char *input_file;

      /* Parse options and command line arguments. Diplay help 
         message if no (or more than one) argument is given. */

      int opt;

      static struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {"verbose", no_argument, NULL, 'v'},
        {"quiet", no_argument, NULL, 'q'},
#ifdef WITH_MPI
        {"mpi_grain", required_argument, 0, 'm' },
#endif
        {0, 0, 0, 0}
      };

      while ((opt = getopt_long (argc, argv, "hVvqm:", longopts, NULL)) != -1)
        {
          switch (opt)
            {
            case 'h':
              usage ();
              exit (0);
              break;
            case 'V':
              version ();
              exit (0);
              break;
            case 'v':
              verbose = 2;
              break;
            case 'q':
              verbose = 0;
              break;
#ifdef WITH_MPI
            case 'm':
              mpi_grain = atoi( optarg );
              break;
#endif
            default:
              usage ();
              exit (1);
            }
        };
      argc -= optind;
      argv += optind;
      if (argc != 1)
        {
          usage ();
          exit (1);
        }
      input_file = argv[0];

      /* Read the input file just to get file names */
      read_input_file_names (input_file, &input_params.files, verbose);

      /* Read the chemical network file */
      read_network (input_params.files.chem_file, &network, verbose);

      /* Read the input file */
      read_input (input_file, &input_params, &network, verbose);

      /* Read the source model file */
      read_source (input_params.files.source_file, &source_mdl, &input_params,
                   verbose);

      /* Allocate results */
      alloc_results (&results, input_params.output.time_steps,
                     source_mdl.n_cells,
                     input_params.output.n_output_species);
#ifdef WITH_MPI
      //Sending data
      int n_work;
      int cnt;
        
      // Checking mpi_grain
      if( mpi_grain <= 0 || mpi_grain >  source_mdl.n_cells )
        {
          if( verbose == 1 )
            {
              fprintf (stdout,
                       "astrochem: %s:%d: mpi_grain %i invalid, setting to %i\n",
                       __FILE__, __LINE__, mpi_grain, source_mdl.n_cells / nb_worker );
            }
          mpi_grain = source_mdl.n_cells / nb_worker;
        }

      for (cnt = 0; cnt < nb_worker; cnt++)
        {
          n_work = cnt + MPI_FIRST_WORKER;

          //Send all structures and array pointed by structures to each worker
          MPI_Send (&input_params, sizeof (inp_t), MPI_BYTE, n_work, tag,
                    MPI_COMM_WORLD);
          MPI_Send (input_params.abundances.initial_abundances,
                    input_params.abundances.n_initial_abundances *
                    sizeof (abund_t), MPI_BYTE, n_work, tag, MPI_COMM_WORLD);
          MPI_Send (input_params.output.output_species_idx,
                    input_params.output.n_output_species * sizeof (int),
                    MPI_BYTE, n_work, tag, MPI_COMM_WORLD);
          MPI_Send (&network, sizeof (net_t), MPI_BYTE, n_work, tag,
                    MPI_COMM_WORLD);
          MPI_Send (network.species_names,
                    network.n_alloc_species * sizeof (species_name_t),
                    MPI_BYTE, n_work, tag, MPI_COMM_WORLD);
          MPI_Send (network.reactions, network.n_reactions * sizeof (react_t),
                    MPI_BYTE, n_work, tag, MPI_COMM_WORLD);
          MPI_Send (&source_mdl.ts.n_time_steps, 1, MPI_INT, n_work, tag,
                    MPI_COMM_WORLD);
          MPI_Send (source_mdl.ts.time_steps, source_mdl.ts.n_time_steps,
                    MPI_DOUBLE, n_work, tag, MPI_COMM_WORLD);
          MPI_Send (&source_mdl.mode, sizeof (SOURCE_MODE), MPI_BYTE, n_work,
                    tag, MPI_COMM_WORLD);
          MPI_Send (&mpi_grain, 1, MPI_INT, n_work, tag, MPI_COMM_WORLD);
        }

      int current_cell = 0;
      int idx_array[mpi_grain];
      char tmp_msg;
      int i, n_cell, stopped_worker = 0;
      //MPI Master Loop
      while (stopped_worker < nb_worker)
        {
          //Wait for worker
          MPI_Recv (&tmp_msg, 1, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
                    MPI_COMM_WORLD, &status);
          // a worker is ready
          if (tmp_msg == 'R')
            {
              // no more cell to send
              if (current_cell == source_mdl.n_cells)
                {
                  //stop worker
                  tmp_msg = 'S';
                  MPI_Send (&tmp_msg, 1, MPI_CHAR, status.MPI_SOURCE, tag,
                            MPI_COMM_WORLD);
                  stopped_worker++;
                }
              //still got cell to send
              else
                {
                  tmp_msg = 'W';
                  //Inform worker that master gonna send cells
                  MPI_Send (&tmp_msg, 1, MPI_CHAR, status.MPI_SOURCE, tag,
                            MPI_COMM_WORLD);
                  //Prepare cells
                  n_cell = 0;
                  for (i = 0; i < mpi_grain; i++)
                    {
                      if (current_cell == source_mdl.n_cells)
                        {
                          idx_array[i] = -1;
                        }
                      else
                        {
                          idx_array[i] = current_cell;
                          current_cell++;
                          n_cell++;
                        }
                    }
                  if (verbose >= 1)
                    {
                      printf ("Sending %i to %i cell to worker %i\n",
                              idx_array[0], idx_array[mpi_grain - 1],
                              status.MPI_SOURCE);
                    }
                  //Send cell idxs
                  MPI_Send (idx_array, mpi_grain, MPI_INT,
                            status.MPI_SOURCE, tag, MPI_COMM_WORLD);
                  //Send cell data
                  MPI_Send (source_mdl.cell[idx_array[0]].nh,
                            sizeof (double) * source_mdl.ts.n_time_steps * 4 *
                            n_cell, MPI_BYTE, status.MPI_SOURCE, tag,
                            MPI_COMM_WORLD);
                }
            }
          //A worker got results
          else
            {
              n_work = status.MPI_SOURCE;
              //Receive cell idxs
              MPI_Recv (idx_array, mpi_grain, MPI_INT, n_work, MPI_ANY_TAG,
                        MPI_COMM_WORLD, &status);
              //Check cell idxs
              int n_cell = 0;
              for (i = 0; i < mpi_grain; i++)
                {
                  if (idx_array[i] != -1)
                    {
                      n_cell++;
                      if (i > 0 && idx_array[i] != idx_array[i - 1] + 1)
                        {
                          fprintf (stderr,
                                   "astrochem: %s:%d: index received are not consecutive\n",
                                   __FILE__, __LINE__);
                          exit (1);
                        }
                    }
                }
              //Receive results
              MPI_Recv (&results.abundances[get_abundance_idx
                                            (&results, idx_array[0], 0, 0)],
                        n_cell * results.n_time_steps *
                        results.n_output_abundances, MPI_DOUBLE, n_work,
                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              MPI_Recv (&results.routes[get_route_idx
                                        (&results, idx_array[0], 0, 0, 0)],
                        N_OUTPUT_ROUTES * n_cell * results.n_time_steps *
                        results.n_output_abundances * sizeof (rout_t),
                        MPI_BYTE, n_work, tag, MPI_COMM_WORLD, &status);
            }
        }

      // Write ouput
      output (source_mdl.n_cells, &input_params, &source_mdl, &network,
              &results, verbose);
      free_input (&input_params);
      free_mdl (&source_mdl);
      free_network (&network);
      free_results (&results);
    }
  else                          // Worker
    {
      inp_t input;
      net_t network;
      res_t results;
      SOURCE_MODE mode;

      //Receive all structures and pointed array

      //Input
      MPI_Recv (&input, sizeof (inp_t), MPI_BYTE, MPI_MASTER_RANK, tag,
                MPI_COMM_WORLD, &status);
      //Allocate input arrays before reception
      alloc_input (&input, input.abundances.n_initial_abundances,
                   input.output.n_output_species);
      MPI_Recv (input.abundances.initial_abundances,
                input.abundances.n_initial_abundances * sizeof (abund_t),
                MPI_BYTE, MPI_MASTER_RANK, tag, MPI_COMM_WORLD, &status);
      MPI_Recv (input.output.output_species_idx,
                input.output.n_output_species * sizeof (int), MPI_BYTE,
                MPI_MASTER_RANK, tag, MPI_COMM_WORLD, &status);

      //Network
      MPI_Recv (&network, sizeof (net_t), MPI_BYTE, MPI_MASTER_RANK, tag,
                MPI_COMM_WORLD, &status);
      //Allocate network arrays before reception
      alloc_network (&network, network.n_alloc_species, network.n_reactions);
      MPI_Recv (network.species_names,
                network.n_alloc_species * sizeof (species_name_t), MPI_BYTE,
                MPI_MASTER_RANK, tag, MPI_COMM_WORLD, &status);
      MPI_Recv (network.reactions, network.n_reactions * sizeof (react_t),
                MPI_BYTE, MPI_MASTER_RANK, tag, MPI_COMM_WORLD, &status);

      //Time steps
      time_steps_t ts;
      MPI_Recv (&ts.n_time_steps, 1, MPI_INT, MPI_MASTER_RANK, tag,
                MPI_COMM_WORLD, &status);
      //Allocate time steps before reception
      if ((ts.time_steps =
           malloc (sizeof (double) * ts.n_time_steps)) == NULL)
        {
          fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
                   __FILE__, __LINE__);
          exit (1);
        }
      MPI_Recv (ts.time_steps, ts.n_time_steps, MPI_DOUBLE, MPI_MASTER_RANK,
                tag, MPI_COMM_WORLD, &status);

      //source mode
      MPI_Recv (&mode, sizeof (SOURCE_MODE), MPI_BYTE, MPI_MASTER_RANK, tag,
                MPI_COMM_WORLD, &status);

      //mpi grain
      MPI_Recv (&mpi_grain, 1, MPI_INT, MPI_MASTER_RANK, tag,
                MPI_COMM_WORLD, &status);


      /*Allocate a standard memory aligned block of cells 
         See alloc_source comments for memory aligned allocation explanation */
      cell_block_t block;
      double *data;
      if ((block.cell_idxs = malloc (mpi_grain * sizeof (int))) == NULL)
        {
          fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
                   __FILE__, __LINE__);
          exit (1);
        }
      if ((data =
           malloc (4 * mpi_grain * ts.n_time_steps * sizeof (double))) ==
          NULL)
        {
          fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
                   __FILE__, __LINE__);
          exit (1);
        }
      if ((block.cells = malloc (sizeof (cell_t) * mpi_grain)) == NULL)
        {
          fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
                   __FILE__, __LINE__);
          exit (1);
        }
      int i;
      for (i = 0; i < mpi_grain; i++)
        {
          block.cells[i].nh = &(data[4 * i * ts.n_time_steps]);
          block.cells[i].av = &(data[(4 * i + 1) * ts.n_time_steps]);
          block.cells[i].tgas = &(data[(4 * i + 2) * ts.n_time_steps]);
          block.cells[i].tdust = &(data[(4 * i + 3) * ts.n_time_steps]);
        }

      //Allocate standard results block
      alloc_results (&results, ts.n_time_steps, mpi_grain, input.output.n_output_species);    //Dynamic solve incorrect alloc

      int n_block = 0;
      while (1)                 //MPI Worker Loop
        {
          char msg = 'R';
          //Send Ready flag to master
          MPI_Send (&msg, 1, MPI_CHAR, MPI_MASTER_RANK, tag, MPI_COMM_WORLD);
          //Receive master response
          MPI_Recv (&msg, 1, MPI_CHAR, MPI_MASTER_RANK, MPI_ANY_TAG,
                    MPI_COMM_WORLD, &status);
          //No more work, need to stop
          if (msg == 'S')
            {
              break;
            }
          //Still got work, receive cell indexes
          MPI_Recv (block.cell_idxs, mpi_grain, MPI_INT, MPI_MASTER_RANK,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          int n_cell = 0;
          for (i = 0; i < mpi_grain; i++)
            {
              if (block.cell_idxs[i] != -1)
                {
                  n_cell++;
                }
            }
          //Receive cells data
          MPI_Recv (block.cells[0].nh,
                    sizeof (double) * ts.n_time_steps * 4 * n_cell, MPI_BYTE,
                    MPI_MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule (dynamic, 1)
#endif
          // Solve the ODE system for each cell. 
          for (cell_index = 0; cell_index < n_cell; cell_index++)
            {
              if (block.cell_idxs[cell_index] != -1)
                {
                  full_solve (cell_index, &input, mode, &block.cells[cell_index],
                         &network, &ts, &results, 0);
                }
            }
          msg = 'F';
          // Send Finished flag to master
          MPI_Send (&msg, 1, MPI_CHAR, MPI_MASTER_RANK, tag, MPI_COMM_WORLD);
          //Send cells indexes to master
          MPI_Send (block.cell_idxs, mpi_grain, MPI_INT, MPI_MASTER_RANK,
                    tag, MPI_COMM_WORLD);
          //Send results to master
          MPI_Send (results.abundances,
                    n_cell * results.n_time_steps *
                    results.n_output_abundances, MPI_DOUBLE, MPI_MASTER_RANK,
                    tag, MPI_COMM_WORLD);
          MPI_Send (results.routes,
                    N_OUTPUT_ROUTES * n_cell * results.n_time_steps *
                    results.n_output_abundances * sizeof (rout_t), MPI_BYTE,
                    MPI_MASTER_RANK, tag, MPI_COMM_WORLD);
          n_block++;
        }
      //Free structures
      free_results (&results);
      free (block.cells[0].nh);
      free (block.cells);
      free (block.cell_idxs);
      free (ts.time_steps);
      free_network (&network);
      free_input (&input);
    }
  MPI_Finalize ();
#else //WITH_MPI

      /* Solve the ODE system for each cell. */
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule (dynamic, 1)
#endif
      for (cell_index = 0; cell_index < source_mdl.n_cells; cell_index++)
        {
          if (verbose >= 1)
            fprintf (stdout, "Computing abundances in cell %d...\n",
                     cell_index);
            full_solve (cell_index, &input_params, source_mdl.mode,
                 &source_mdl.cell[cell_index], &network, &source_mdl.ts,
                 &results, verbose);
          if (verbose >= 1)
            fprintf (stdout, "Done with cell %d.\n", cell_index);
        }

      /* Write the abundances in output files */
      output (source_mdl.n_cells, &input_params, &source_mdl, &network,
              &results, verbose);
      free_input (&input_params);
      free_mdl (&source_mdl);
      free_network (&network);
      free_results (&results);

#endif // WITH_MPI
      return (EXIT_SUCCESS);
    }

/*
   Display help message.
 */

void
usage (void)
{
  fprintf (stdout, "Usage: astrochem [option...] [file]\n\n");
  fprintf (stdout, "Options:\n");
  fprintf (stdout, "   -h, --help         Display this help\n");
  fprintf (stdout, "   -V, --version      Print program version\n");
  fprintf (stdout, "   -v, --verbose      Verbose mode\n");
  fprintf (stdout, "   -q, --quiet        Suppress all messages\n");
#ifdef WITH_MPI
  fprintf (stdout, "   -m, --mpi_grain value    mpi grain to use\n");
#endif
  fprintf (stdout, "\n");
  fprintf (stdout,
           "See the astrochem(1) manual page for more information.\n");
  fprintf (stdout, "Report bugs to <%s>.\n", PACKAGE_BUGREPORT);
}

/*
   Display version.
 */

void
version (void)
{
  fprintf (stdout, "This is astrochem, version %s\n", PACKAGE_VERSION);
#ifdef HAVE_OPENMP
  fprintf (stdout, "OpenMP support enabled, ");
#else
  fprintf (stdout, "OpenMP support disabled, ");
#endif
#ifdef USE_LAPACK
  fprintf (stdout, "LAPACK support enabled.\n");
#else
  fprintf (stdout, "LAPACK support disabled.\n");
#endif
  fprintf (stdout, "Copyright (c) 2006-2013 Sebastien Maret\n");
  fprintf (stdout, "\n");
  fprintf (stdout,
           "This is free software. You may redistribute copies of it under the terms\n");
  fprintf (stdout,
           "of the GNU General Public License. There is NO WARRANTY, to the extent\n");
  fprintf (stdout, "permitted by law.\n");
}
