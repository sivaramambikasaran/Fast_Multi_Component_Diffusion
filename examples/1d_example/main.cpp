#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "domain.hpp"

int main()
{   
    int N_grid = 100;
    domain D(N_grid);
	// The following array stores the name of the files from which data is to be read:
	std::string filename[] = {"data/table_224.txt","data/table_448.txt","data/table_896.txt",
                              "data/table_1792.txt","data/table_3585.txt","data/table_7171.txt"};

    // Used for timing:
	double tic, toc;
    // Tolerance for approximate methods:
    double tolerance = 1e-12;

	// Iterating over each of the files:
	for (int k = 0; k < 6; ++k) 
	{
        domain D(N_grid);

        // Reading the species from the data file:
        D.readSpecies(filename[k].c_str());
        std::cout << "Number of Species             :" << D.N_species << std::endl;
        // Getting the mole, mass fraction and gradient profiles:
        D.generateSpeciesProfile();
        // Getting the temperature, pressure profiles and their gradients:
        D.generateTemperaturePressureProfile();

        // Get velocities using the fast method:
        tic = omp_get_wtime();
        D.computeSpeciesVelocityFast(tolerance);
        toc = omp_get_wtime();
        std::cout << "Time Taken By Fast Method     :" << (toc - tic) << std::endl;

        // Getting exact velocities using A / b:
        tic = omp_get_wtime();
        D.computeSpeciesVelocityExact();
        toc = omp_get_wtime();
        std::cout << "Time Taken By Direct Method   :" << (toc - tic) << std::endl;
        
        // Getting velocities using conjugate gradient:
        tic = omp_get_wtime();
        D.computeSpeciesVelocityIterative(tolerance);
        toc = omp_get_wtime();
        std::cout << "Time Taken By Iterative Method:" << (toc - tic) << std::endl;

        std::cout << "Error(Fast Method)            :" << D.getError() << std::endl;
        std::cout << "========================================" << std::endl;
    }

    return 0;
}
