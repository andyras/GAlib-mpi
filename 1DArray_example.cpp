#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ga-mpi/ga.h>
#include <ga-mpi/std_stream.h>
#include "mpi.h"
#include <sys/wait.h>
#include <string>

// This are the declaration of the objective functions which are defined later.
float objective(GAGenome &);
float dynamixObjective(GAGenome &);

// declare Initializer also
void Initializer(GAGenome &);

int mpi_tasks, mpi_rank;

int main(int argc, char **argv)
{
  // MPI init
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // See if we've been given a seed to use (for testing purposes).  When you
  // specify a random seed, the evolution will be exactly the same each time
  // you use that seed number
  unsigned int seed = 0;
  for(int i=1 ; i<argc ; i++)
    if(strcmp(argv[i++],"seed") == 0)
      seed = atoi(argv[i]);

  // Declare variables for the GA parameters and set them to some default values.
  int popsize  = 2; // Population
  int ngen     = 2; // Generations
  float pmut   = 0.03;
  float pcross = 0.65;

  // popsize / mpi_tasks must be an integer
  popsize = mpi_tasks * int((double)popsize/(double)mpi_tasks+0.999);

  // Create the phenotype for two variables.  The number of bits you can use to
  // represent any number is limited by the type of computer you are using.
  // For this case we use 10 bits for each var, ranging the square domain [0,5*PI]x[0,5*PI]
  ///GABin2DecPhenotype map;
  ///GABin2DecPhenotype map;
  ///map.add(10, 0.0, 5.0 * M_PI);
  ///map.add(10, 0.0, 5.0 * M_PI);

  // Create the template genome using the phenotype map we just made.
  ///GABin2DecGenome genome(map, objective);
  //GA1DArrayGenome<double> genome(2, objective);
  GA1DArrayGenome<double> genome(2, dynamixObjective);
  // define own initializer, can do the same for mutator and comparator
  genome.initializer(::Initializer);

  // Now create the GA using the genome and run it. We'll use sigma truncation
  // scaling so that we can handle negative objective scores.
  GASimpleGA ga(genome); // TODO change to steady-state
  GALinearScaling scaling;
  ga.minimize();		// by default we want to minimize the objective
  ga.populationSize(popsize);
  ga.nGenerations(ngen);
  ga.pMutation(pmut);
  ga.pCrossover(pcross);
  ga.scaling(scaling);
  if(mpi_rank == 0)
    ga.scoreFilename("evolution.txt");
  else
    ga.scoreFilename("/dev/null");
  ga.scoreFrequency(1);
  ga.flushFrequency(1);
  ga.selectScores(GAStatistics::AllScores);
  // Pass MPI data to the GA class
  ga.mpi_rank(mpi_rank);
  ga.mpi_tasks(mpi_tasks);
  ga.evolve(seed);

  // Dump the GA results to file
  if(mpi_rank == 0)
  {
    genome = ga.statistics().bestIndividual();
    printf("GA result:\n");
    printf("x = %f, y = %f\n",
	genome.gene(0), genome.gene(1));
  }

  MPI_Finalize();

  return 0;
}

float objective(GAGenome &c)
{
  ///GABin2DecGenome &genome = (GABin2DecGenome &)c;
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;
  float x, y, error;

  x = genome.gene(0);
  y = genome.gene(1);

  // Function with local minima. The lowest is located at (5/2*PI, 5/2*PI)
  error = ((1.-sin(x)*sin(y))+sqrt((x-M_PI*2.5)*(x-M_PI*2.5)+(y-M_PI*2.5)*(y-M_PI*2.5))/10.0)/2.5;

  return error;
}

float dynamixObjective(GAGenome &c) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;
  // get MPI rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // variables and output
  double g1 = genome.gene(0);
  double g2 = genome.gene(1);
  double g1_c = genome.gene(2);
  float output = 1.0;

  pid_t pid;
  int status;
  char ** args;
  std::string arg = "/extra/scratch/foo";

  // names for original directory and job directory
  // TODO pass as cmd line params
  pid = getpid();
  std::string tmpDir ("/extra/scratch/dynamix");
  // job directory name
  char jobFmt [] = "%s%.12e_%.12e_%.12e_%d_%d";
  char jobDir [1000];
  sprintf(jobDir, jobFmt, tmpDir.c_str(), g1, g2, g1_c, rank, size);
  std::cout << jobDir << std::endl;
  std::string dynDir ("/home/andyras/git/dynamix/dm/");
  std::string dynInsDir (dynDir + "ins/");

  // ---- set up job directory ---- //
 
  // -- create job directory -- //

  args = new char * [4];
  args[0] = new char [6];
  strncpy(args[0], "mkdir", 6);
  args[1] = new char [3];
  strncpy(args[1], "-p", 3);
  args[2] = new char [arg.length()+1];
  strncpy(args[2], arg.c_str(), arg.length()+1);
  args[3] = NULL;

  fprintf(stdout, "COMMAND:");
  for (int ii = 0; ii < 4; ii++) {
    fprintf(stdout, " %s", args[ii]);
  }
  fprintf(stdout, "\n");

  if ((pid = fork()) < 0) { // fork fails
    fprintf(stdout, "Fork bork\n");
    _exit(EXIT_FAILURE);
  }
  else if (pid == 0) { // child process
    // make directory
    execv("/bin/mkdir", args);
    // just in case
    _exit(EXIT_FAILURE);
  }
  else { // parent process
    waitpid(pid, &status, 0);
  }
  // clean up
  for (int ii = 0; ii < 4; ii++) {
    delete [] args[ii];
  }
  delete [] args;

  // -- copy inputs to job dir -- //

  args = new char * [5];
  args[0] = new char [strlen("cp") + 1];
  strncpy(args[0], "cp", strlen("cp") + 1);
  args[1] = new char [strlen("-r") + 1];
  strncpy(args[1], "-r", strlen("-r") + 1);
  args[2] = new char [dynInsDir.length() + 1];
  strncpy(args[2], dynInsDir.c_str(), dynInsDir.length() + 1);
  args[3] = new char [strlen(jobDir) + 1];
  strncpy(args[3], jobDir, strlen(jobDir) + 1);
  args[4] = NULL;

  fprintf(stdout, "COMMAND:");
  for (int ii = 0; ii < 5; ii++) {
    fprintf(stdout, " %s", args[ii]);
  }
  fprintf(stdout, "\n");

  if ((pid = fork()) < 0) { // fork fails
    fprintf(stdout, "Fork bork\n");
    _exit(EXIT_FAILURE);
  }
  else if (pid == 0) { // child
    execv("/bin/cp", args);
    _exit(EXIT_FAILURE);
  }
  else { // parent
    waitpid(pid, &status, 0);
  }

  // clean up
  for (int ii = 0; ii < 5; ii++) {
    delete [] args[ii];
  }
  delete [] args;

  // change parameters in input file

  // ---- run code ---- //
  //
  // ---- check for success ---- //
  //
  // ---- read in outputs ---- //
  //
  // ---- calculate objective ---- //
  //
  // ---- remove job directory ---- //

  return output;
}

void Initializer(GAGenome &g) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;

  // there are two genes
  genome.gene(0, GARandomFloat(0.0,5*M_PI));
  genome.gene(1, GARandomFloat(0.0,5*M_PI));

  return;
}
