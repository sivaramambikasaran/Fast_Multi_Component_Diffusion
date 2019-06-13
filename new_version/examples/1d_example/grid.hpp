#ifndef __grid_point_hpp__
#define __grid_point_hpp__

#include <vector>
#include <set>
#include "Eigen/Dense"

class grid_point_point
{
    friend class domain;

public:

    //  Constructor
    grid_point(double x, int N_species);

    //  X-location of the grid_point
    double x;
    //  Number of species
    int N_species;

    //  Species Velocity, where species_velocity(i) denotes the velocity of species 'i' at the grid_point
    Eigen::VectorXd species_velocity, exact_species_velocity, interative_species_velocity;
    //  Error in solution
    double error, iterative_error;
    //  Low-rank form of the inverse diffusion coefficient matrix
    Eigen::MatrixXd L, R;
    //  Inverse Diffusion Coefficient matrix
    Eigen::MatrixXd V;
    //  Rank of the inverse diffusion coefficient matrix
    int rank;

    //  Number of iterations for the iterative solver
    int N_iterations;

    //  Computes the right hand side
    void computeRHS();

    //  Computes the species velocity
    void computeSpeciesVelocity(const std::vector<species>& my_species, const Eigen::VectorXd& W, double tolerance);
    //  Computes the species velocity iteratively
    void computeSpeciesVelocityIterative(const std::vector<species>& my_species, const Eigen::VectorXd& W, double tolerance);
    //  Computes the exact species velocity
    void computeSpeciesVelocityExact(const std::vector<species>& my_species, const Eigen::VectorXd& W);

    //  Compute diffusion coefficient using kinetic theory
    double computeInverseDiffusionCoefficientKineticTheory(int i, int j, const std::vector<species>& my_species);

    //  Compute the collision integral
    double computeCollisionIntegral(int i, int j, const std::vector<species>& my_species);

    //  Compute inverse diffusion coefficient using collision integral
    double computeInverseDiffusionCoefficientCollision(int i, int j, const std::vector<species>& my_species);

    //  Compute inverse diffusion coefficient returns either Kinetic theory result or Collision integral result
    double computeInverseDiffusionCoefficient(int i, int j, const std::vector<species>& my_species);

};

grid_point::grid_point(double x, int N_species) 
{
    this->x         = x;
    this->N_species = N_species;
}

double grid_point::computeInverseDiffusionCoefficientKineticTheory(int i, int j, const std::vector<species>& my_species) 
{
    //  1/D =   ((d_i+d_j)/(2d_{max}))^2*sqrt(2*m_i*m_j/(m_i+m_j)/m_{max})
    double dia_sum   =   (my_species[i].diameter+my_species[j].diameter);
    double temp     =   dia_sum*dia_sum*sqrt(2*my_species[i].mass*my_species[j].mass/(my_species[i].mass+my_species[j].mass));
    return  temp;
}

double grid_point::computeCollisionIntegral(int i, int j, const std::vector<species>& my_species) {
    static const double m[] =   {6.8728271691,9.4122316321,7.7442359037,0.23424661229,1.45337701568,5.2269794238,9.7108519575,0.46539437353,0.00041908394781};
    // double T =   kB*temperature/sqrt(my_species[i].energy*my_species[j].energy);
    double T    =   temperature*sqrt(my_species[i].energy*my_species[j].energy);
    double Nr   =   m[0]+T*(m[1]+T*(m[2]+T*m[3]));
    double Dr   =   m[4]+T*(m[5]+T*(m[6]+T*(m[7]+T*m[8])));
    return (Nr/Dr);
}

double grid_point::computeInverseDiffusionCoefficientCollision(int i, int j, const std::vector<species>& my_species) 
{
    return computeInverseDiffusionCoefficientKineticTheory(i,j,my_species)*computeCollisionIntegral(i,j,my_species);
}

double grid_point::computeInverseDiffusionCoefficient(int i, int j, const std::vector<species>& my_species) 
{
    return computeInverseDiffusionCoefficientKineticTheory(i,j,my_species);
}


void grid_point::computeSpeciesVelocityExact(const std::vector<species>& my_species, const Eigen::VectorXd& W) 
{
    obtain_Inverse_Diffusion_Coefficients(my_species);
    Eigen::MatrixXd M       =   mole_fraction.asDiagonal()*V + Eigen::VectorXd::Random(N_species)*W.transpose() - Eigen::MatrixXd((V*mole_fraction).asDiagonal());
    exact_species_velocity    =   M.fullPivLu().solve(rhs);
}

void grid_point::computeSpeciesVelocityIterative(const std::vector<species>& my_species, const Eigen::VectorXd& W, double tolerance) {
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