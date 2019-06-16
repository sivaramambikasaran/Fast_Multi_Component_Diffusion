#ifndef __grid_point_hpp__
#define __grid_point_hpp__

#include <vector>
#include <set>
#include "FMDV.hpp"
#include "Eigen/IterativeLinearSolvers"

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

class grid_point
{
    friend class domain;

    // X-location of the grid_point
    double x;
    // Solver object used to quickly solve for the diffusion velocities:
    fast_solver* FS;

public:

    // Constructor
    grid_point(double x);

    // Species Velocity, where species_velocity(i) denotes the velocity of species 'i' at the grid_point
    Eigen::VectorXd species_velocity, exact_species_velocity, iterative_species_velocity;
    // Error in solution
    double error;

    // Computes the species velocity
    void computeSpeciesVelocityFast(double tolerance);
    // Computes the exact species velocity
    void computeSpeciesVelocityExact();
    // Computes the species velocity iteratively
    void computeSpeciesVelocityIterative(double tolerance);
};

grid_point::grid_point(double x) 
{
    this->x  = x;
}

void grid_point::computeSpeciesVelocityFast(double tolerance)
{
    species_velocity = FS->computeSpeciesVelocities(tolerance);
} 

void grid_point::computeSpeciesVelocityExact() 
{
    Eigen::MatrixXd V = FS->getInverseDiffusionCoefficientMatrix();
    Eigen::VectorXd X = FS->mole_fraction;
    Eigen::VectorXd W = FS->molecular_mass;
    Eigen::VectorXd D = FS->thermal_diffusivities;
    Eigen::VectorXd Y = FS->mass_fraction;
    Eigen::VectorXd S = FS->getRandomVector();

    // (diag(VX) - diag(X)V - S * Wt)z = (b - alpha * s)
    Eigen::MatrixXd M       = Eigen::MatrixXd((V * X).asDiagonal()) - X.asDiagonal() * V - S * W.transpose();
    Eigen::VectorXd z_exact = M.fullPivLu().solve(FS->getRHS());

    double prefactor = FS->temperature_gradient / (FS->density * FS->temperature);
    exact_species_velocity = Eigen::VectorXd::Zero(X.size()); 
    for(int i = 0; i < X.size(); i++)
    {
        if(X(i) == 0)
        {
            exact_species_velocity(i) = 0;
        }

        else
        {
            exact_species_velocity(i) = z_exact(i) / X(i) - prefactor * D(i) / Y(i);
        }
    }
}

void grid_point::computeSpeciesVelocityIterative(double tolerance) 
{
    Eigen::MatrixXd V = FS->getInverseDiffusionCoefficientMatrix();
    Eigen::VectorXd X = FS->mole_fraction;
    Eigen::VectorXd W = FS->molecular_mass;
    Eigen::VectorXd D = FS->thermal_diffusivities;
    Eigen::VectorXd Y = FS->mass_fraction;
    Eigen::VectorXd S = FS->getRandomVector();

    // (diag(VX) - diag(X)V - S * Wt)z = (b - alpha * s)
    Eigen::MatrixXd M = Eigen::MatrixXd((V * X).asDiagonal()) - X.asDiagonal() * V - S * W.transpose();
    // solve Ax = b using CG with matrix-free version:
    Eigen::BiCGSTAB<Eigen::MatrixXd> cg;
    cg.setTolerance(100 * tolerance);
    cg.compute(M);
    Eigen::VectorXd z = cg.solve(FS->getRHS());

    double prefactor = FS->temperature_gradient / (FS->density * FS->temperature);
    iterative_species_velocity = Eigen::VectorXd::Zero(X.size()); 
    for(int i = 0; i < X.size(); i++)
    {
        if(X(i) == 0)
        {
            iterative_species_velocity(i) = 0;
        }

        else
        {
            iterative_species_velocity(i) = z(i) / X(i) - prefactor * D(i) / Y(i);
        }
    }
}

#endif /*__grid_point_hpp__*/
