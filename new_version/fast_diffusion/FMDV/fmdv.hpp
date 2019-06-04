#ifndef __fmdv_hpp__
#define __fmdv_hpp__

#include <iostream>
#include <vector>
#include "Eigen/Dense"

// Definition of various physical constants:
const double R        =	8.3144598;
const double AVAGADRO =	6.022140857e23;
const double kB       =	R / AVAGADRO;
const double PI       =	3.141592653589793238;

class fmdv 
{
// These methods are used internally
// Not needed for a user to access:
private:
	// Number of species
	int N_species;
    // Number of grid points:
    int N_grid;
    
    // The following variables are vectors of size N_grid:
    // Density at the grid points:
    Eigen::VectorXd density;
    // Temperature at the grid points:
	Eigen::VectorXd temperature;
    // Temperature gradient at the grid points:
	Eigen::VectorXd temperature_gradient;
    // Pressure at the grid points:
    Eigen::VectorXd pressure;
    // Pressure gradient at the grid points:
    Eigen::VectorXd pressure_gradients;

    // The following variables are vectors of size N_species:
    // Molecular weights:
    Eigen::VectorXd molecular_weights;
    // Diameters:
    Eigen::VectorXd diameters;
    // Thermal diffusivities:
    Eigen::VectorXd thermal_diffusivities;
    // Collision energies:
    Eigen::VectorXd collision_energies;

    // The following variables are matrices of size N_species X N_grid
    // Mole fraction:
    Eigen::MatrixXd mole_fraction;
    // Mole fraction gradient:
    Eigen::MatrixXd mole_fraction_gradient;
    // Mass fraction:
    Eigen::MatrixXd mass_fraction;
    // Body forces:
    Eigen::MatrixXd body_forces;

    // RHS that needs to be solved:
    Eigen::VectorXd rhs;
	// Low-rank form of the inverse diffusion coefficient matrix
	Eigen::MatrixXd U ,V;

	// Computes the right hand side
	void computeRHS();
    
	//	Compute diffusion coefficient using kinetic theory
	double compute_Inverse_Diffusion_Coefficient_Kinetic_Theory(int i, int j, const std::vector<species>& mySpecies);

	//	Compute the collision integral
	double compute_Collision_Integral(int i, int j, const std::vector<species>& mySpecies);

	//	Compute inverse diffusion coefficient using collision integral
	double compute_Inverse_Diffusion_Coefficient_Collision(int i, int j, const std::vector<species>& mySpecies);

	//	Compute inverse diffusion coefficient returns either Kinetic theory result or Collision integral result
	double compute_Inverse_Diffusion_Coefficient(int i, int j, const std::vector<species>& mySpecies);

	//	Obtain row of the inverse diffusion coefficient matrix
	Eigen::VectorXd get_Row(int i, int startCol, int colSize, const std::vector<species>& mySpecies);

	//	Obtain column of the inverse diffusion coefficient matrix
	Eigen::VectorXd get_Col(int i, int startRow, int rowSize, const std::vector<species>& mySpecies);

	//	Obtains a submatrix of the inverse diffusion coefficient matrix
	Eigen::MatrixXd get_Matrix(int startRow, int startCol, int nRows, int nCols, const std::vector<species>& mySpecies);

	//	Obtain the maximum entry in absolute value of a vector
	unsigned max_Abs_Vector(const Eigen::VectorXd& v, const std::set<int>& allowed_Indices, double& max);

	//	Compressed inverse diffusion coefficient matrix
	void compressed_Inverse_Diffusion_Coefficients(const std::vector<species>& mySpecies, double tolerance);

	//	Obtain inverse diffusion coefficient matrix
	void obtain_Inverse_Diffusion_Coefficients(const std::vector<species>& mySpecies);

	//	Obtain the diagonal matrix for solving the diffusion velocities
	void compute_Diagonal();

public:

	// Constructor:
    // Initializing using number of grid points and the number of species:
	fmdv(int N_grid, int N_species);

    // Setters:
    // The length of this vector would be N_grid
    void setPressure(Eigen::VectorXd pressure);
    // The length of this vector would be N_grid
    void setPressureGradient(Eigen::VectorXd pressure_gradient);
    // The length of this vector would be N_grid
    void setDensity(Eigen::VectorXd density);
    // The length of this vector would be N_grid
    void setTemperature(Eigen::VectorXd temperature);
    // The length of this vector would be N_grid
    void setTemperatureGradient(Eigen::VectorXd temperature_gradient);
    // The length of this vector would be N_species
    void setMolecularWeights(Eigen::VectorXd molecular_weights);
    // The length of this vector would be N_species
    void setDiameters(Eigen::VectorXd diameters);
    // The length of this vector would be N_species
    void setThermalDiffusivities(Eigen::VectorXd thermal_diffusivities);
    // The dimensions of this matrix would be N_species X N_grid:
    void setMoleFractions(Eigen::MatrixXd mole_fractions);
    // The dimensions of this matrix would be N_species X N_grid:
    void setMoleFractionGradients(Eigen::MatrixXd mole_fractions);
    // The dimensions of this matrix would be N_species X N_grid:
    void setBodyForces(Eigen::MatrixXd body_forces);

    // This is the main method that the user utilizes to get the Species Velocities:
    Eigen::VectorXd computeSpeciesVelocities(double tolerance);
};

fmdv::fmdv(int N_grid, int N_species)
{
    this->N_grid    = N_grid;
    this->N_species = N_species;

    return;
}

void fmdv::setPressure(Eigen::VectorXd pressure);
{
    this->pressure = pressure;
    return;
}

void fmdv::setPressureGradient(Eigen::VectorXd pressure_gradient);
{
    this->pressure_gradient = pressure_gradient;
    return;
}

void fmdv::setDensity(Eigen::VectorXd density);
{
    this->density = density;
    return;
}

void fmdv::setTemperature(Eigen::VectorXd temperature);
{
    this->temperature = temperature;
    return;
}

void fmdv::setTemperatureGradient(Eigen::VectorXd temperature_gradient);
{
    this->temperature_gradient = temperature_gradient;
    return;
}

void fmdv::setMolecularWeights(Eigen::VectorXd molecular_weights);
{
    this->molecular_weights = molecular_weights;
    return;
}

void fmdv::setDiameters(Eigen::VectorXd diameters);
{
    this->diameters = diameters;
    return;
}

void fmdv::setThermalDiffusivities(Eigen::VectorXd thermal_diffusivities);
{
    this->thermal_diffusivities = thermal_diffusivities;
    return;
}

void fmdv::setMoleFractions(Eigen::MatrixXd mole_fractions);
{
    this->mole_fractions = mole_fractions;
    return;
}

void fmdv::setBodyForces(Eigen::MatrixXd body_forces);
{
    this->body_forces = body_forces;
    return;
}

// This computes the RHS for the i-th grid cell:
void fmdv::computeRHS(int i) 
{
    rhs	= moleFractionGradient - (massFraction - moleFraction) * pressureGradient/pressure;
    rhs	= rhs * temperature / pressure * sqrt(temperature);
}

void fmdv::maxAbsVector(const Vec& v, 
                        const std::set<int>& allowed_indices,
                        dtype& max, int& index
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

void fmdv::rookPiv(Mat& L, Mat& R, double tolerance_or_rank,
                   int n_row_start, int n_col_start,
                   int n_rows, int n_cols
                  )
{
    // Indices which have been used:
    std::vector<int> row_ind;      
    std::vector<int> col_ind;

    // Indices that are remaining:
    std::set<int> remaining_row_ind;
    std::set<int> remaining_col_ind;
    
    // Bases:
    std::vector<Vec> u; 
    std::vector<Vec> v;

    for(int k = 0; k < n_rows; k++) 
    {
        remaining_row_ind.insert(k);
    }
    
    for(int k = 0; k < n_cols; k++) 
    {
        remaining_col_ind.insert(k);
    }

    dtype max, gamma;

    // Initialize the matrix norm and the the first row index
    dtype_base matrix_norm = 0;
    row_ind.push_back(0);
    remaining_row_ind.erase(0);

    // Stores the pivot entry of the considered row / col:
    int pivot;

    int target_rank = 0;
    // This would get updated:
    int computed_rank = 0;
    Vec row, col;

    double tolerance = 0;
    // These quantities in finding the stopping criteria:
    dtype_base row_squared_norm, row_norm, col_squared_norm, col_norm;

    if(tolerance_or_rank < 1)
        tolerance = tolerance_or_rank;
    else
        target_rank = tolerance_or_rank;

    // So these would be particularly useful for poorly conditioned matrices:
    int max_tries = 10;
    int count;

    // Repeat till the desired tolerance / rank is obtained
    do 
    {
        // Generation of the row
        // Row of the residuum and the pivot column
        // By calling row_ind.back(), we are getting the last pushed number
        row = this->A->getRow(n_row_start + row_ind.back(), n_col_start, n_cols);
        
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
            row = this->A->getRow(n_row_start + new_row_ind, n_col_start, n_cols);
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
        gamma = dtype_base(1.0) / max;

        // Generation of the column
        // Column of the residuum and the pivot row
        col = this->A->getCol(n_col_start + col_ind.back(), n_row_start, n_rows);
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
            col = this->A->getCol(n_col_start + new_col_ind, n_row_start, n_rows);
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
    while(((tolerance_or_rank < 1) ? 
          computed_rank * (n_rows + n_cols) * row_norm * col_norm > 
          fabs(max) * tolerance * matrix_norm 
          : computed_rank < target_rank)
          && 
          computed_rank < fmin(n_rows, n_cols)
         );

    // If the computed_rank is >= to full-rank
    // then return the trivial full-rank decomposition
    if (computed_rank >= fmin(n_rows, n_cols) - 1) 
    {
        if (n_rows < n_cols) 
        {
            L = Mat::Identity(n_rows, n_rows);
            R = this->A->getMatrix(n_row_start, n_col_start, n_rows, n_cols).transpose();
            computed_rank = n_rows;
        }

        else 
        {
            L = this->A->getMatrix(n_row_start, n_col_start, n_rows, n_cols);
            R = Mat::Identity(n_cols, n_cols);
            computed_rank = n_cols;
        }
    }
    
    // This is when ACA has succeeded:
    else 
    {
        L = Mat(n_rows, computed_rank);
        R = Mat(n_cols, computed_rank);
        
        for (int j = 0; j < computed_rank; j++) 
        {
            L.col(j) = u[j];
            R.col(j) = v[j];
        }
    }
}

// This method solves the equation:
// (diag(VX) - diag(X)V - S * Wt)z = (b - alpha * s)
Eigen::MatrixXd computeSpeciesVelocities(double tolerance)
{
	Eigen::VectorXd species_velocities = Eigen::VectorXd::Zero(this->N_species);
    Eigen::VectorXd S                  = Eigen::VectorXd::Random(N_species);
    // Step #1 of Algorithm 1: expresses V = L * Rt
    this->compressedInverseDiffusionCoefficients(tolerance);
    // Step #2 of Algorithm 1: D = diag(L * Rt * X)
    this->computeDiagonal();

    // This computes the RHS of the equation (b - alpha * s) and stores it in the variable rhs:
	this->computeRHS();
    // As stated in the paper, we declare variables P and Q:
    // P = [diag(X)L, S] ∈ R^{N_species × (r+1)}
    // Q = [ R, W ] ∈ R^{N_species × (r+1)}
	Eigen::MatrixXd P(N_species, rank + 1), Q(N_species, rank + 1);
	P << moleFraction.asDiagonal() * L, S;
	Q << R, W;

	Eigen::VectorXd b_tilde = -invdiagonal.cwiseProduct(rhs);
	Eigen::VectorXd P_tilde =  invdiagonal.asDiagonal() * P;

	Eigen::VectorXd b_bar = Q.transpose() * b_tilde;
	Eigen::VectorXd P_bar = Q.transpose() * P_tilde;

	Eigen::VectorXd = Eigen::MatrixXd::Identity(rank+1,rank+1) - Vtilde.transpose()*Utilde;
	speciesVelocity      = tempVelocity + Utilde*temp.fullPivLu().solve();
}
