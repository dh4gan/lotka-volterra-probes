/*
 * parFile.cpp
 *
 *  Created on: Sep 23, 2013
 *
 *      Author: davidharvey
 */

#include "parFile.h"
#include "Constants.h"
#include <iostream> // Included for debug lines only
#include <math.h>
#include <string>
#include <ctime>

parFile::parFile()
    {
      nVertices=0;
  	iseed=-1;
  	initialPrey=0.0;
  	initialPred=0.0;
  	tmax=0.0;
  	dt = 0.0;
  	tsnap = 0.0;

  	icChoice= "GHZ"; // 0 = GHZ, 1= cluster
  	parChoice= "constant"; // "constant" = constant parameters, "uniform" = random sampling, "gaussian" =Gaussian sampling

  	rmin=0.0;
  	rmax=0.0;
  	scaleLength=0.0;

  	range=0.0;


  	// Variable 1 is the constant or minimum value selected, or the mean for gaussian sampling
  	// Variable 2 is the maximum value for uniform sampling, or the standard deviation for Gaussian sampling

  	predGrow1=0.0;
  	predGrow2=0.0;

  	predDeath1=0.0;
  	predDeath2=0.0;

  	preyGrow1=0.0;
  	preyGrow2=0.0;

  	preyDeath1=0.0;
  	preyDeath2=0.0;

  	mutationRate1=0.0;
  	mutationRate2=0.0;

  	outflowRate1=0.0;
  	outflowRate2=0.0;
  	velocity1=0.0;
  	velocity2=0.0;


    }

parFile::parFile(string name):parFile()
    {
  parFileName = name;
    }

void parFile::readParams()
{
  /*
   * Written 22/2/18 by dh4gan
   * reads parameters from file
   *
   */

  string line, par;
  ifstream myfile(parFileName.c_str());

  // Then loop through each line using getline and then
  // assign to vectors


  while (getline(myfile, line))
  	{
  	istringstream iss(line);
  	iss >> par;


  	// Read domain specific variables

  	  if(par== "nVertices"){iss >> nVertices;}
  	  if(par== "range"){iss >> range;}
  	  if(par== "icChoice"){iss >> icChoice;}
  	  if(par== "rmin"){iss >> rmin;}
  	  if(par== "rmax"){iss >> rmax;}
  	  if(par== "scaleLength"){iss >> scaleLength;}

  	  if(par== "tmax"){iss >> tmax;}
  	  if(par== "dt"){iss >> dt;}
  	  if(par=="tsnap"){iss >> tsnap;}

  	  if(par== "iseed"){iss>>iseed;}

  	  // Read LK system variables

  	  if(par== "parChoice"){iss>>parChoice;}

  	  if(par=="initialPrey"){iss>>initialPrey;}
  	  if(par=="initialPred"){iss>>initialPred;}

  	  if(par=="preyGrow"){ iss>>preyGrow1 >> preyGrow2;}
  	  if(par=="preyDeath"){ iss>>preyDeath1 >> preyDeath2;}

  	  if(par=="predGrow"){ iss>>predGrow1 >> predGrow2;}
  	  if(par=="predDeath"){ iss>>predDeath1 >> predDeath2;}
  	  if(par=="mutationRate"){ iss>>mutationRate1 >> mutationRate2;}
  	  if(par=="outflowRate") { iss>>outflowRate1 >> outflowRate2;}
  	  if(par=="velocity"){ iss>>velocity1 >> velocity2;}

  	  if(velocity1>1.0)
  	    {
  	      printf("WARNING: attempting to set probe velocity > c; automatically fixing max=c \n");
  	      velocity1 = 1.0;
  	    }


  	  if(velocity2>1.0)
  	    {
  	      printf("WARNING: attempting to set probe velocity > c; automatically fixing max=c \n");
  	      velocity2 = 1.0;
  	    }

  	  velocity1 = velocity1*c_kpc_Myr;
  	  velocity2 = velocity2*c_kpc_Myr;

  	}
}

void parFile::writeParamsToFile()
{
/*
 * Written 22/2/18 by dh4gan
 * Writes input parameters to file
 *
 */

  string outputFileString = parFileName+".info";
  FILE* outputFile = fopen(outputFileString.c_str(), "w");

  time_t rawtime;
  struct tm* timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  string divider = "---------------------";

  fprintf(outputFile,"Parameters for lotka_volterra_probes run \n");
  fprintf(outputFile,"Read from %s\n",parFileName.c_str());
  fprintf(outputFile, "Time: %s \n", asctime(timeinfo));
  fprintf(outputFile, "%s \n",divider.c_str());
  fprintf(outputFile, "Domain Parameters\n");
  fprintf(outputFile, "%s \n",divider.c_str());
  fprintf(outputFile,"Number of stars: %i \n",nVertices);
  fprintf(outputFile,"Connection Range: %f kpc \n",range);
  fprintf(outputFile,"Initial Condition Choice: %s \n",icChoice.c_str());
  fprintf(outputFile,"Minimum Radius: %f kpc \n",rmin);
  fprintf(outputFile,"Maximum Radius: %f kpc \n",rmax);

  if(icChoice=="GHZ")
    {
      fprintf(outputFile,"Scale Length: %f kpc \n",scaleLength);
    }

  fprintf(outputFile,"Maximum Time: %f Myr \n",tmax);
  fprintf(outputFile,"Timestep: %f Myr \n",dt);
  fprintf(outputFile,"Snapshot Interval: %f Myr \n",tsnap);

  fprintf(outputFile, "%s \n",divider.c_str());
  fprintf(outputFile, "Lotka Volterra Parameters\n");
  fprintf(outputFile, "%s \n",divider.c_str());

  fprintf(outputFile,"Parameter Choice: %s \n",parChoice.c_str());

  if(parChoice=="constant")
    {
      fprintf(outputFile, "Prey Growth Rate: %f \n",preyGrow1);
      fprintf(outputFile, "Prey Death Rate: %f \n",preyDeath1);

      fprintf(outputFile, "Predator Growth Rate: %f \n",predGrow1);
      fprintf(outputFile, "Predator Death Rate: %f \n",predDeath1);
    }
  else if(parChoice=="uniform")

    {
      fprintf(outputFile, "Minimum Prey Growth Rate: %f \n",preyGrow1);
      fprintf(outputFile, "Maximum Prey Growth Rate: %f \n",preyGrow2);

      fprintf(outputFile, "Minimum Prey Death Rate: %f \n",preyDeath1);
      fprintf(outputFile, "Maximum Prey Death Rate: %f \n",preyDeath2);

      fprintf(outputFile, "Minimum Predator Growth Rate: %f \n",predGrow1);
      fprintf(outputFile, "Maximum Predator Growth Rate: %f \n",predGrow2);

      fprintf(outputFile, "Minimum Predator Death Rate: %f \n",predDeath1);
      fprintf(outputFile, "Maximum Predator Death Rate: %f \n",predDeath2);
    }

  else if(parChoice=="gaussian")

      {
        fprintf(outputFile, "Mean Prey Growth Rate: %f \n",preyGrow1);
        fprintf(outputFile, "SD of Prey Growth Rate: %f \n",preyGrow2);

        fprintf(outputFile, "Mean Prey Death Rate: %f \n",preyDeath1);
        fprintf(outputFile, "SD of Prey Death Rate: %f \n",preyDeath2);

        fprintf(outputFile, "Mean Predator Growth Rate: %f \n",predGrow1);
        fprintf(outputFile, "SD of Predator Growth Rate: %f \n",predGrow2);

        fprintf(outputFile, "Mean Predator Death Rate: %f \n",predDeath1);
        fprintf(outputFile, "SD of Predator Death Rate: %f \n",predDeath2);
      }

  fclose(outputFile);
}


void parFile::writeParamsToScreen()

{
  /*
   * Written 22/2/18 by dh4gan
   * Writes parameter data to screen
   *
   */

    time_t rawtime;
    struct tm* timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    string divider = "---------------------";

   printf("Parameters for lotka_volterra_probes run \n");
   printf("Read from %s\n",parFileName.c_str());
   printf("Time: %s \n", asctime(timeinfo));
   printf("%s \n",divider.c_str());
   printf("Domain Parameters\n");
   printf("%s \n",divider.c_str());
   printf("Number of stars: %i \n",nVertices);
   printf("Connection Range: %f kpc \n",range);
   printf("Initial Condition Choice: %s \n",icChoice.c_str());
   printf("Minimum Radius: %f kpc \n",rmin);
   printf("Maximum Radius: %f kpc \n",rmax);

   if(icChoice=="GHZ")
     {
       printf("Scale Length: %f kpc \n",scaleLength);
     }

   printf("Maximum Time: %f Myr \n",tmax);
   printf("Timestep: %f Myr \n",dt);
   printf("Snapshot Interval: %f Myr \n",tsnap);

   printf("%s \n",divider.c_str());
   printf("Lotka Volterra Parameters\n");
   printf("%s \n",divider.c_str());

   printf("Parameter Choice: %s \n",parChoice.c_str());

   if(parChoice=="constant")
     {
       printf("Prey Growth Rate: %f \n",preyGrow1);
       printf("Prey Death Rate: %f \n",preyDeath1);

       printf("Predator Growth Rate: %f \n",predGrow1);
       printf("Predator Death Rate: %f \n",predDeath1);
     }
   else if(parChoice=="uniform")

     {
       printf("Minimum Prey Growth Rate: %f \n",preyGrow1);
       printf("Maximum Prey Growth Rate: %f \n",preyGrow2);

       printf("Minimum Prey Death Rate: %f \n",preyDeath1);
       printf("Maximum Prey Death Rate: %f \n",preyDeath2);

       printf("Minimum Predator Growth Rate: %f \n",predGrow1);
       printf("Maximum Predator Growth Rate: %f \n",predGrow2);

       printf("Minimum Predator Death Rate: %f \n",predDeath1);
       printf("Maximum Predator Death Rate: %f \n",predDeath2);
     }

   else if(parChoice=="gaussian")

       {
         printf("Mean Prey Growth Rate: %f \n",preyGrow1);
         printf("SD of Prey Growth Rate: %f \n",preyGrow2);

         printf("Mean Prey Death Rate: %f \n",preyDeath1);
         printf("SD of Prey Death Rate: %f \n",preyDeath2);

         printf("Mean Predator Growth Rate: %f \n",predGrow1);
         printf("SD of Predator Growth Rate: %f \n",predGrow2);

         printf("Mean Predator Death Rate: %f \n",predDeath1);
         printf("SD of Predator Death Rate: %f \n",predDeath2);
       }

printf("%s \n",divider.c_str());
}






