#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "domain.hpp"



// This object is used to store data and solve the diffusion velocities for all species at a single grid point:
class grid_point
{
    // Number of species
    int N_species;

public:

    // Constructor
    grid_point(std::vector<species> species)
    {
        this->x         = x;
        this->N_species = N_species;
    }

    // Species Velocity, where species_velocity(i) denotes the velocity of species 'i' at the grid_point
    Eigen::VectorXd species_velocity_fast, species_velocity_exact, species_velocity_iterative;

    // Computes the species velocity
    void computeSpeciesVelocity(double tolerance);
    // Computes the species velocity iteratively
    void computeSpeciesVelocityIterative(double tolerance);
    // Computes the exact species velocity
    void computeSpeciesVelocityExact();
};

