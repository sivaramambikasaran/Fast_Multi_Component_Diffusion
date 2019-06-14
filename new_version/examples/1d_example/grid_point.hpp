#ifndef __grid_point_hpp__
#define __grid_point_hpp__

#include <vector>
#include <set>
#include "Eigen/Dense"

class fast_solver : public FMDV
{
public:
    fast_solver(int N_species) : FMDV(N_species)

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

class grid_point
{
    friend class domain;
    // Solver object used to quickly solve for the diffusion velocities:
    fast_solver FS;

    // Density, Temperature and Pressure values at the grid:
    double density, pressure, pressure_gradient, temperature, temperature_gradient;
    // Species Profile:
    Eigen::VectorXd mole_fraction, mass_fraction, mole_fraction_gradient, body_forces;
    // Species Properties:
    Eigen::VectorXd molecular_mass, diameters, collision_energy, thermal_diffusivities; 

    // X-location of the grid_point
    double x;
    // Number of species
    int N_species;

    // Inverse Diffusion Coefficient matrix
    Eigen::MatrixXd V;

public:

    // Constructor
    grid_point(double x);

    // Species Velocity, where species_velocity(i) denotes the velocity of species 'i' at the grid_point
    Eigen::VectorXd species_velocity, exact_species_velocity, iterative_species_velocity;
    // Error in solution
    double error;

    // Number of iterations for the iterative solver
    int N_iterations;

    // Computes the species velocity
    void computeSpeciesVelocityFast(double tolerance);
    // Computes the exact species velocity
    void computeSpeciesVelocityExact();
    // Computes the species velocity iteratively
    void computeSpeciesVelocityIterative(double tolerance);
};

grid_point::grid_point(double x, int N_species) 
{
    this->x         = x;
    this->N_species = N_species;
    this->FS        = fast_solver(N_species);
}

void grid_point::computeSpeciesVelocityFast(double tolerance)
{
    Eigen::VectorXd species_velocity = FS.computeSpeciesVelocities(tolerance);
} 


void grid_point::computeSpeciesVelocityExact() 
{
    obtain_Inverse_Diffusion_Coefficients(my_species);
    Eigen::MatrixXd M       =   mole_fraction.asDiagonal()*V + Eigen::VectorXd::Random(N_species)*W.transpose() - Eigen::MatrixXd((V*mole_fraction).asDiagonal());
    exact_species_velocity    =   M.fullPivLu().solve(rhs);
}

void grid_point::computeSpeciesVelocityIterative(double tolerance) 
{
    // obtain_Inverse_Diffusion_Coefficients(my_species);
    Eigen::VectorXd temp    =   V*mole_fraction;
    Eigen::MatrixXd M       =   mole_fraction.asDiagonal()*V + Eigen::VectorXd::Random(N_species)*V.transpose();
    for (int j=0; j<N_species; ++j) {
        M(j,j)-=temp(j);
    }

    // solve Ax = b using CG with matrix-free version:
    Eigen::BiCGSTAB<Eigen::MatrixXd > cg;
    cg.setTolerance(100*tolerance);
    cg.compute(M);
    interative_species_velocity    =   cg.solve(rhs);
    this->iterative_error    =   interative_species_velocity.cwiseProduct(W).sum(); //cg.error();
    this->N_iterations       =   cg.iterations();
}

#endif /*__grid_point_hpp__*/
