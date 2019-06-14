#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "domain.hpp"

int main()
{   
    int N_grid = 10;
    domain D(N_grid);

    // Reading the species from the data file:
    D.readSpecies("data/table_224.txt");
	// Getting the mole, mass fraction and gradient profiles:
	D.generateSpeciesProfile();
	// Getting the temperature, pressure profiles and their gradients:
	D.generateTemperaturePressureProfile();

    // Get velocities using the fast method:
    D.computeSpeciesVelocityFast(1e-12);
    // Getting exact velocities using A / b:
    D.computeSpeciesVelocityExact();

    std::cout << "Error:" << D.getError() << std::endl;
    return 0;
}
