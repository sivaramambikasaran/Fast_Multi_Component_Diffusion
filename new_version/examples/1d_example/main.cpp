#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "domain.hpp"

class fast_solver : public FMDV
{
public:
    fast_solver(int N_species) : FMDV(N_species)
    {}

    double getInverseDiffusionCoefficient(int i, int j)
    {
        // 1/D = K * P / T**1.5 * ((d_i+d_j))^2 * sqrt(m_i*m_j/(m_i+m_j))
        double dia_sum              = diameters(i) + diameters(j);
        double elastic_matrix_entry =   (pressure / (temperature * sqrt(temperature)))
                                      * dia_sum * dia_sum * sqrt(molecular_mass(i) * molecular_mass(j) 
                                      / (molecular_mass(i) + molecular_mass(j)));

        // This segment computes collision integral based on collision energies:
        if(collision_energies(i) != 0 && collision_energies(j) != 0)
        {
            static const double m[] = {6.8728271691,  9.4122316321,  7.7442359037,
                                       0.23424661229, 1.45337701568, 5.2269794238,
                                       9.7108519575,  0.46539437353, 0.00041908394781};
            
            double T  = temperature * sqrt(collision_energies(i) * collision_energies(j));
            double Nr = m[0] + T * (m[1] + T * (m[2] + T *  m[3]));
            double Dr = m[4] + T * (m[5] + T * (m[6] + T * (m[7] + T * m[8])));
        
            return (Nr/Dr) * elastic_matrix_entry;
        }

        else
        {
            return elastic_matrix_entry;
        }
    }
};


// This object is used to store data and solve the diffusion velocities for all species at a single grid point:
class grid_point
{
    // Solver object used to quickly solve for the diffusion velocities:
    fast_solver FS;
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

