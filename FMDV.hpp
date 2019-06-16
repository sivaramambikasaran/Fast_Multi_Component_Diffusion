#ifndef __FMDV_HPP__
#define __FMDV_HPP__

#include <iostream>
#include <vector>
#include <set>
#include "Eigen/Dense"

// An instance of this class can be used to solve
// for all the diffusion velocities at a single grid point
class FMDV
{
// These methods are used internally
// Not needed for a user to access:
// private:
public:
	// Number of species
	int N_species;
    // Rank of inverse diffusion matrix:
    int rank;

    // Scalar α which is used in the computation:
    double alpha;
    // Average mass of all species:
    double average_molecular_mass;
    
    // RHS that needs to be solved for:
    Eigen::VectorXd rhs;
	// Low-rank form of the inverse diffusion coefficient matrix
	Eigen::MatrixXd L, R;
    // Diagonal matrix diag(LR^T X) and its inverse:
    Eigen::VectorXd diagonal, invdiagonal;
    // Random vector S to use in the computation:
    Eigen::VectorXd S;

    // The following methods are used to get the low rank decomposition:
    Eigen::VectorXd getRow(const int j);
    Eigen::VectorXd getCol(const int k);
    void maxAbsVector(const Eigen::VectorXd& v, 
                      const std::set<int>& allowed_indices,
                      double& max, int& index
                     );
    void compressInverseDiffusionCoefficients(double tolerance);

	// Computes the right hand side
    void computeAlpha();
    void computeDiagonal();
	void computeRHS();

// public:

    // The following variables are doubles:
    // Density at the grid point:
    double density;
    // Temperature at the grid point:
	double temperature;
    // Temperature gradient at the grid point:
	double temperature_gradient;
    // Pressure at the grid point:
    double pressure;
    // Pressure gradient at the grid point:
    double pressure_gradient;

    // The following variables are vectors of size N_species:
    // Molecular weights:
    Eigen::VectorXd molecular_mass;
    // Diameters:
    Eigen::VectorXd diameters;
    // Thermal diffusivities:
    Eigen::VectorXd thermal_diffusivities;
    // Collision energies:
    Eigen::VectorXd collision_energies;
    // Mole fraction:
    Eigen::VectorXd mole_fraction;
    // Mole fraction gradient:
    Eigen::VectorXd mole_fraction_gradient;
    // Mass fraction:
    Eigen::VectorXd mass_fraction;
    // Body forces:
    Eigen::VectorXd body_forces;

	// Constructor:
    // Initializing using number of grid points and the number of species:
	FMDV(int N_species);

    // Method that defines the matrix entries:
    virtual double getInverseDiffusionCoefficient(int i, int j)
    {
        return 0;
    }

    // This is the main method that the user utilizes to get the species velocities:
    Eigen::VectorXd computeSpeciesVelocities(double tolerance);

    // Provided as a means to compare against the naive computation:
    Eigen::MatrixXd getInverseDiffusionCoefficientMatrix();
    Eigen::VectorXd getRandomVector();
    Eigen::VectorXd getRHS();
};

// Used to set the number of species considered in the computation:
FMDV::FMDV(int N_species)
{
    srand(time(NULL));
    this->N_species = N_species;

    molecular_mass         = Eigen::VectorXd::Zero(N_species);
    diameters              = Eigen::VectorXd::Zero(N_species);
    thermal_diffusivities  = Eigen::VectorXd::Zero(N_species);
    collision_energies     = Eigen::VectorXd::Zero(N_species);
    mole_fraction          = Eigen::VectorXd::Zero(N_species);
    mole_fraction_gradient = Eigen::VectorXd::Zero(N_species);
    mass_fraction          = Eigen::VectorXd::Zero(N_species);
    body_forces            = Eigen::VectorXd::Zero(N_species);
    rhs                    = Eigen::VectorXd::Zero(N_species);
    S                      = Eigen::VectorXd::Random(N_species); // NOTE: DO NOT SET TO ZERO. THIS CAN LEAD TO FAILURE!

    return;
}

Eigen::VectorXd FMDV::getRow(const int j) 
{
    Eigen::VectorXd row(N_species);
    #pragma omp parallel for
    for(int k = 0; k < N_species; k++) 
    {   
        row(k) = this->getInverseDiffusionCoefficient(j,  k);
    }

    return row;
}

Eigen::VectorXd FMDV::getCol(const int k) 
{
    Eigen::VectorXd col(N_species);
    #pragma omp parallel for
    for (int j = 0; j < N_species; j++) 
    {
        col(j) = this->getInverseDiffusionCoefficient(j, k);
    }

    return col;
}

Eigen::MatrixXd FMDV::getInverseDiffusionCoefficientMatrix() 
{
    Eigen::MatrixXd mat(N_species, N_species);
    
    #pragma omp parallel for
    for (int j=0; j < N_species; ++j) 
    {
        #pragma omp parallel for
        for (int k=0; k < N_species; ++k) 
        {
            mat(j,k) = this->getInverseDiffusionCoefficient(j, k);
        }
    }

    return mat;
}

void FMDV::maxAbsVector(const Eigen::VectorXd& v, 
                        const std::set<int>& allowed_indices,
                        double& max, int& index
                       )
{
    std::set<int>::iterator it;
    index = *allowed_indices.begin();
    max   = v(index);

    for(it = allowed_indices.begin(); it != allowed_indices.end(); it++) 
    {
        if(fabs(v(*it))>fabs(max)) 
        {
            index   =   *it;
            max     =   v(index);
        }
    }
}

// This is the same optimized ACA routine that is used from HODLR
void FMDV::compressInverseDiffusionCoefficients(double tolerance)
{
    // Indices which have been used:
    std::vector<int> row_ind;      
    std::vector<int> col_ind;

    // Indices that are remaining:
    std::set<int> remaining_row_ind;
    std::set<int> remaining_col_ind;
    
    // Bases:
    std::vector<Eigen::VectorXd> u; 
    std::vector<Eigen::VectorXd> v;

    for(int k = 0; k < N_species; k++) 
    {
        remaining_row_ind.insert(k);
    }
    
    for(int k = 0; k < N_species; k++) 
    {
        remaining_col_ind.insert(k);
    }

    double max, gamma;

    // Initialize the matrix norm and the the first row index
    double matrix_norm = 0;
    row_ind.push_back(0);
    remaining_row_ind.erase(0);

    // Stores the pivot entry of the considered row / col:
    int pivot;

    int target_rank = 0;
    // This would get updated:
    int computed_rank = 0;
    Eigen::VectorXd row, col;

    // These quantities in finding the stopping criteria:
    double row_squared_norm, row_norm, col_squared_norm, col_norm;

    // So these would be particularly useful for poorly conditioned matrices:
    int max_tries = 10;
    int count;

    // Repeat till the desired tolerance / rank is obtained
    do 
    {
        // Generation of the row
        // Row of the residuum and the pivot column
        // By calling row_ind.back(), we are getting the last pushed number
        row = this->getRow(row_ind.back());
        
        for(int i = 0; i < computed_rank; i++) 
        {
            row = row - u[i](row_ind.back()) * v[i];
        }

        this->maxAbsVector(row, remaining_col_ind, max, pivot);
        count = 0;

        // Alternating upon each call:
        bool eval_at_end = false;
        // Toggling randomness
        bool use_randomization = true;

        // This while loop is needed if in the middle of the algorithm the 
        // row happens to be exactly the linear combination of the previous rows 
        // upto some tolerance. i.e. prevents from ACA throwing false positives
        while (fabs(max) < tolerance && 
               count < max_tries   && 
               remaining_col_ind.size() > 0 && 
               remaining_row_ind.size() > 0
              ) 
        {
            row_ind.pop_back();
            int new_row_ind;

            // When rank < 3, we will just choose entries from the ends of the matrix:
            if(computed_rank < 3)
            {
                if(eval_at_end == true)
                {
                    new_row_ind = *remaining_row_ind.end();
                }

                else
                {
                    new_row_ind = *remaining_row_ind.begin();
                }
    
                eval_at_end = !(eval_at_end);
            }

            // However, when we have rank >=3, we will choose the entries such that
            // the newly picked entry is at the mid-point of the already chosen ones:
            else
            {
                if(use_randomization == true)
                {
                    std::set<int>::const_iterator it(remaining_row_ind.begin());
                    std::advance(it, rand() % remaining_row_ind.size());
                    new_row_ind = *it;
                }

                else
                {
                    std::vector<int> row_ind_sort(row_ind);
                    std::sort(row_ind_sort.begin(), row_ind_sort.end());
                    std::vector<int> row_ind_diff(row_ind_sort.size() - 1);

                    int max = 0;
                    int idx = 0;

                    for(int i = 0; i < row_ind_sort.size() - 1; i++)
                    {
                        row_ind_diff[i] = row_ind_sort[i+1] - row_ind_sort[i];
                        if(row_ind_diff[i] > max)
                        {
                            idx = i;
                            max = row_ind_diff[i];
                        }
                    }

                    new_row_ind = row_ind_sort[idx] + max / 2;
                }

                use_randomization = !(use_randomization);
            }

            row_ind.push_back(new_row_ind);
            remaining_row_ind.erase(new_row_ind);
            // Generation of the row
            // Row of the residuum and the pivot column
            row = this->getRow(new_row_ind);
            for(int i = 0; i < computed_rank; i++) 
            {
                row = row - u[i](row_ind.back()) * v[i];
            }

            this->maxAbsVector(row, remaining_col_ind, max, pivot);
            count++;
        }

        // In case it failed to resolve in the previous step, 
        // we break out of the dowhile loop:
        if (count == max_tries || 
            remaining_col_ind.size() == 0 || 
            remaining_row_ind.size() == 0
           )
        {
            break;
        } 

        // Resetting count back to zero for columns:
        count = 0;
        
        col_ind.push_back(pivot);
        remaining_col_ind.erase(pivot);
        // Normalizing constant
        gamma = double(1.0) / max;

        // Generation of the column
        // Column of the residuum and the pivot row
        col = this->getCol(col_ind.back());
        for(int i = 0; i < computed_rank; i++) 
        {
            col = col - v[i](col_ind.back()) * u[i];
        }

        this->maxAbsVector(col, remaining_row_ind, max, pivot);
        // Repeating the same randomization we carried out for the rows, now for the columns:
        while (fabs(max)<tolerance && 
               count < max_tries && 
               remaining_col_ind.size() >0 && 
               remaining_row_ind.size() >0
              ) 
        {
            col_ind.pop_back();

            int new_col_ind;

            if(col_ind.size() < 3)
            {
                if(eval_at_end)
                {
                    new_col_ind = *remaining_col_ind.end();
                }

                else
                {
                    new_col_ind = *remaining_col_ind.begin();
                }
    
                eval_at_end = !eval_at_end;
            }

            else
            {                
                if(use_randomization == true)
                {
                    std::set<int>::const_iterator it(remaining_col_ind.begin());
                    std::advance(it, rand() % remaining_col_ind.size());
                    new_col_ind = *it;
                }

                else
                {
                    std::vector<int> col_ind_sort(col_ind);
                    std::sort(col_ind_sort.begin(), col_ind_sort.end());
                    std::vector<int> col_ind_diff(col_ind_sort.size() - 1);

                    int max = 0;
                    int idx = 0;

                    for(int i = 0; i < col_ind_sort.size() - 1; i++)
                    {
                        col_ind_diff[i] = col_ind_sort[i+1] - col_ind_sort[i];
                        if(col_ind_diff[i] > max)
                        {
                            idx = i;
                            max = col_ind_diff[i];
                        }
                    }

                    new_col_ind = col_ind_sort[idx] + max / 2;
                }

                use_randomization = !(use_randomization);
            }

            col_ind.push_back(new_col_ind);
            remaining_col_ind.erase(new_col_ind);

            // Generation of the column
            // Column of the residuum and the pivot row:
            col = this->getCol(new_col_ind);
            for(int i = 0; i < computed_rank; i++) 
            {
                col = col - v[i](col_ind.back()) * u[i];
            }

            this->maxAbsVector(col, remaining_row_ind, max, pivot);
            count++;
        }

        row_ind.push_back(pivot);
        remaining_row_ind.erase(pivot);

        // New vectors
        u.push_back(gamma * col);
        v.push_back(row);

        // New approximation of matrix norm
        row_squared_norm = row.squaredNorm();
        row_norm         = sqrt(row_squared_norm);

        col_squared_norm = col.squaredNorm();
        col_norm         = sqrt(col_squared_norm);

        // Updating the matrix norm:
        matrix_norm += std::abs(gamma * gamma * row_squared_norm * col_squared_norm);

        for(int j = 0; j < computed_rank; j++) 
        {
            matrix_norm += 2.0 * std::abs(u[j].dot(u.back())) 
                               * std::abs(v[j].dot(v.back()));
        }
        
        computed_rank++;
    }
    while(computed_rank * (N_species + N_species) * row_norm * col_norm > 
          fabs(max) * tolerance * matrix_norm 
          && 
          computed_rank < fmin(N_species, N_species)
         );

    // If the computed_rank is >= to full-rank
    // then return the trivial full-rank decomposition
    if(computed_rank >= N_species - 1) 
    {
        L = this->getInverseDiffusionCoefficientMatrix();
        R = Eigen::MatrixXd::Identity(N_species, N_species);
        computed_rank = N_species;
        this->rank = computed_rank;
    }
    
    // This is when ACA has succeeded:
    else 
    {
        this->rank = computed_rank;

        L = Eigen::MatrixXd(N_species, computed_rank);
        R = Eigen::MatrixXd(N_species, computed_rank);
        
        for (int j = 0; j < computed_rank; j++) 
        {
            L.col(j) = u[j];
            R.col(j) = v[j];
        }
    }
}

void FMDV::computeDiagonal()
{
    diagonal    = L * (R.transpose() * mole_fraction);
	invdiagonal = diagonal.cwiseInverse();
    return;
}

void FMDV::computeAlpha()
{
    average_molecular_mass = (molecular_mass.cwiseProduct(mole_fraction)).sum();
    // α  = W * ∇T / (ρ * T) Σ D
    alpha = (average_molecular_mass * temperature_gradient / (density * temperature)) * thermal_diffusivities.sum();
    return;
}

// This computes the RHS for the expression we want to solve:
// This is given by equations (9) and (20) from the paper
// RHS    : b - αS
// b      : -∇X + (Y - X)∇P/P + ρ/P Y (f - f_tilde) 
// f_tilde: Σ (Y * f) 
void FMDV::computeRHS() 
{
    double f_tilde = mass_fraction.dot(body_forces);
    rhs	= - (mole_fraction_gradient + (mass_fraction - mole_fraction) * pressure_gradient / pressure).array()
          + (density / pressure * mass_fraction.array() * (body_forces.array() - f_tilde).array()).array();
    rhs = rhs - alpha * S;
    return;
}

Eigen::VectorXd FMDV::getRandomVector()
{
    return S;
}

Eigen::VectorXd FMDV::getRHS()
{
    return rhs;
}

// This method solves the equation:
// (diag(VX) - diag(X)V - S * Wt)z = (b - alpha * s)
Eigen::VectorXd FMDV::computeSpeciesVelocities(double tolerance)
{
    // Step #1 of Algorithm 1: expresses V = L * Rt
    // This uses ACA:
    this->compressInverseDiffusionCoefficients(tolerance);
    // Step #2 of Algorithm 1: D = diag(L * Rt * X)
    this->computeDiagonal();
    // As stated in the paper, we declare variables P and Q:
	Eigen::MatrixXd P(N_species, rank + 1), Q(N_species, rank + 1);
    // Step #3 of Algorithm 1:
    // P = [diag(X)L, S] ∈ R^{N_species × (r+1)}
	P << mole_fraction.asDiagonal() * L, S;
    // Step #4 of Algorithm 1: computing α
    this->computeAlpha();
    // Q = [ R, W ] ∈ R^{N_species × (r+1)}
	Q << R, molecular_mass;
    // This computes the RHS of the equation (b - alpha * s) and stores it in the variable rhs:
	this->computeRHS();
    // Step #4 of Algorithm: computing b_tilde = inv(D) * (b - αS)
	Eigen::VectorXd b_tilde = invdiagonal.cwiseProduct(rhs);
    // // Step #5 of Algorithm: computing P_tilde = inv(D) * P
	Eigen::MatrixXd P_tilde = invdiagonal.asDiagonal() * P;
    // // Step #6 of Algorithm: computing b_bar = Q^T * b_tilde
	Eigen::VectorXd b_bar = Q.transpose() * b_tilde;
    // // Step #7 of Algorithm: computing P_bar = Q^T * P_tilde
	Eigen::MatrixXd P_bar = Q.transpose() * P_tilde;
    // // Step #8 of Algorithm: computing z = b_tilde + P_tilde * inv(I - P_bar) * b_bar
	Eigen::VectorXd z = b_tilde + P_tilde * (Eigen::MatrixXd::Identity(rank + 1, rank + 1) - P_bar).fullPivLu().solve(b_bar);
    // Step #9 of Algorithm: computing v
    double prefactor = temperature_gradient / (density * temperature);
    Eigen::VectorXd species_velocities(N_species);
    for(int i = 0; i < N_species; i++)
    {
        // If mole-fraction is zero, doesn't it mean that species doesn't exist and that velocity is zero:
        // TODO: Why is it assigned a finite value in the paper?
        if(mole_fraction(i) == 0)
        {
            species_velocities(i) = 0;
        }

        else
        {
            species_velocities(i) = z(i) / mole_fraction(i) - prefactor * thermal_diffusivities(i) / mass_fraction(i);
        }
    }

    return species_velocities;
}

#endif // __FMDV_HPP__
