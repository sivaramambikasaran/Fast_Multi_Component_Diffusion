#ifndef __grid_hpp__
#define __grid_hpp__

// const double R			=	8.3144598;
// const double Avagadro	=	6.022140857e23;
// const double kB			=	R/Avagadro;
// const double PI			=	3.141592653589793238;

#include<cassert>
#include<vector>
#include<set>
#include <Eigen/Dense>

class grid {
	friend class domain;
public:
	//	X-location of the grid
	double x;

	//	Number of species
	int nSpecies;

	//	Total mass at the grid point
	double mass;

	//	Mole fraction of the species, where moleFraction(i) denotes the mole fraction of species 'i' at the grid
	Eigen::VectorXd moleFraction;

	//	Gradient of the mole fraction of the species, where moleFractionGradient(i) denotes the gradient of the mole fraction of species 'i' at the grid
	Eigen::MatrixXd moleFractionGradient;

	//	Mass fraction of the species, where moleFraction(i) denotes the mole fraction of species 'i' at the grid
	Eigen::VectorXd massFraction;

	//	Gradient of the mass fraction of the species, where moleFractionGradient(i) denotes the gradient of the mole fraction of species 'i' at the grid
	Eigen::VectorXd massFractionGradient;

	//	Temperature at the grid
	double temperature;

	//	Pressure at each grid point in the domain
	double pressure;

	//	Pressure gradient at each grid point in the domain
	double pressureGradient;

	//	Species Velocity, where speciesVelocity(i) denotes the velocity of species 'i' at the grid
	Eigen::VectorXd speciesVelocity, exactSpeciesVelocity, iterativespeciesVelocity;

	//	Error in solution
	double error, iterativeerror;

	//	Low-rank form of the inverse diffusion coefficient matrix
	Eigen::MatrixXd U ,V;

	//	Inverse Diffusion Coefficient matrix
	Eigen::MatrixXd D;

	//	Rank of the inverse diffusion coefficient matrix
	int rank;

	//	Right hand side that needs to be solved
	Eigen::VectorXd rhs;

	//	Compute the diagonal matrix required for solving the system
	Eigen::VectorXd diagonal, invdiagonal;

	//	Number of iterations for the iterative solver
	int nIterations;

	//	Computes the right hand side
	void compute_RHS();

	//	Computes the species velocity
	void compute_Species_Velocity(const std::vector<species>& mySpecies, const Eigen::VectorXd& W, double tolerance);

	//	Computes the species velocity iteratively
	void compute_Species_Velocity_Iteratively(const std::vector<species>& mySpecies, const Eigen::VectorXd& W, double tolerance);

	//	Computes the exact species velocity
	void compute_Exact_Species_Velocity(const std::vector<species>& mySpecies, const Eigen::VectorXd& W);

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

	//	Constructor
	grid(double x, int nSpecies);
};

grid::grid(double x, int nSpecies) {
	this->x			=	x;
	this->nSpecies	=	nSpecies;
}

double grid::compute_Inverse_Diffusion_Coefficient_Kinetic_Theory(int i, int j, const std::vector<species>& mySpecies) {
	//	1/D	=	((d_i+d_j)/(2d_{max}))^2*sqrt(2*m_i*m_j/(m_i+m_j)/m_{max})
	double diaSum	=	(mySpecies[i].diameter+mySpecies[j].diameter);
	double temp		=	diaSum*diaSum*sqrt(2*mySpecies[i].mass*mySpecies[j].mass/(mySpecies[i].mass+mySpecies[j].mass));
	return	temp;
}

double grid::compute_Collision_Integral(int i, int j, const std::vector<species>& mySpecies) {
	static const double m[]	=	{6.8728271691,9.4122316321,7.7442359037,0.23424661229,1.45337701568,5.2269794238,9.7108519575,0.46539437353,0.00041908394781};
	// double T	=	kB*temperature/sqrt(mySpecies[i].energy*mySpecies[j].energy);
	double T	=	temperature*sqrt(mySpecies[i].energy*mySpecies[j].energy);
	double Nr	=	m[0]+T*(m[1]+T*(m[2]+T*m[3]));
	double Dr	=	m[4]+T*(m[5]+T*(m[6]+T*(m[7]+T*m[8])));
	return (Nr/Dr);
}

double grid::compute_Inverse_Diffusion_Coefficient_Collision(int i, int j, const std::vector<species>& mySpecies) {
	return compute_Inverse_Diffusion_Coefficient_Kinetic_Theory(i,j,mySpecies)*compute_Collision_Integral(i,j,mySpecies);
}

double grid::compute_Inverse_Diffusion_Coefficient(int i, int j, const std::vector<species>& mySpecies) {
	return compute_Inverse_Diffusion_Coefficient_Kinetic_Theory(i,j,mySpecies);
	// return compute_Inverse_Diffusion_Coefficient_Collision(i,j,mySpecies);
}

Eigen::VectorXd grid::get_Row(int i, int startCol, int colSize, const std::vector<species>& mySpecies) {
	Eigen::VectorXd row(colSize);
	// std::cout << i << "\t" << startCol << "\t" << colSize << "\n";
	#pragma omp parallel for
	for (int j=0; j<colSize; ++j) {
		row(j)	=	compute_Inverse_Diffusion_Coefficient(i,startCol+j,mySpecies);
	}
	return row;
}

Eigen::VectorXd grid::get_Col(int i, int startRow, int rowSize, const std::vector<species>& mySpecies) {
	Eigen::VectorXd col(rowSize);

	#pragma omp parallel for
	for (int j=0; j<rowSize; ++j) {
		col(j)	=	compute_Inverse_Diffusion_Coefficient(j+startRow, i, mySpecies);
	}
	return col;
}

Eigen::MatrixXd grid::get_Matrix(int startRow, int startCol, int nRows, int nCols, const std::vector<species>& mySpecies) {
	Eigen::MatrixXd invD(nRows, nCols);

	#pragma omp parallel for
	for (int j=0; j<nRows; ++j) {
		for (int k=0; k<nCols; ++k) {
			invD(j,k)	=	compute_Inverse_Diffusion_Coefficient(startRow+j,startCol+k, mySpecies);
		}
	}
	return invD;
}

unsigned grid::max_Abs_Vector(const Eigen::VectorXd& v, const std::set<int>& allowed_Indices, double& max) {
	std::set<int>::iterator it	=	allowed_Indices.begin();
	unsigned index	=	*allowed_Indices.begin();
	max				=	v(index);

	for (it=allowed_Indices.begin(); it!=allowed_Indices.end(); ++it) {
		if(fabs(v(*it))>fabs(max)) {
			index	=	*it;
			max		=	v(index);
		}
	}
	return index;
}

void grid::compressed_Inverse_Diffusion_Coefficients(const std::vector<species>& mySpecies, double tolerance) {
	int start_Row	=	0;
	int start_Col	=	0;
	int n_Rows		=	nSpecies;
	int n_Cols		=	nSpecies;

	std::vector<int> rowIndex;		///	This stores the row indices, which have already been used.
	std::vector<int> colIndex;		///	This stores the column indices, which have already been used.
	std::set<int> remainingRowIndex;/// Remaining row indicies
	std::set<int> remainingColIndex;/// Remaining row indicies
	std::vector<Eigen::VectorXd> u;	///	Stores the column basis.
	std::vector<Eigen::VectorXd> v;	///	Stores the row basis.

	for (int k=0; k<n_Rows; ++k) {
		remainingRowIndex.insert(k);
	}
	for (int k=0; k<n_Cols; ++k) {
		remainingColIndex.insert(k);
	}

	srand (time(NULL));
	double max, Gamma, unused_max;

	/*  INITIALIZATION  */

	/// Initialize the matrix norm and the the first row index
	double matrix_Norm  =   0;
	rowIndex.push_back(0);
	remainingRowIndex.erase(0);

	int pivot;

	rank   =   0;

	Eigen::VectorXd a, row, col;

	double row_Squared_Norm, row_Norm, col_Squared_Norm, col_Norm;

	/// Repeat till the desired tolerance is obtained
	do {
		/// Generation of the row
		/// Row of the residuum and the pivot column
		// row =   A.row(rowIndex.back());
		row	=	get_Row(start_Row+rowIndex.back(), start_Col, n_Cols, mySpecies);
		for (int l=0; l<rank; ++l) {
			row =   row-u[l](rowIndex.back())*v[l];
		}

		// std::cout << row << "\n";
		pivot   =   max_Abs_Vector(row, remainingColIndex, max);

		int max_tries  =   10;
		int count      =   0;

		/// This randomization is needed if in the middle of the algorithm the row happens to be exactly the linear combination of the previous rows upto some tolerance.
		while (fabs(max)<tolerance && count < max_tries) {
			rowIndex.pop_back();
			int new_rowIndex	=	*remainingRowIndex.begin();
			rowIndex.push_back(new_rowIndex);
			remainingRowIndex.erase(new_rowIndex);

			/// Generation of the row
			// a	=	A.row(new_rowIndex);
			a	=	get_Row(start_Row+new_rowIndex, start_Col, n_Cols, mySpecies);
			/// Row of the residuum and the pivot column
			row =   a;
			for (int l=0; l<rank; ++l) {
				row =   row-u[l](rowIndex.back())*v[l];
			}
			pivot   =   max_Abs_Vector(row, remainingColIndex, max);
			++count;
		}
		// std::cout << row << "\n";

		if (count == max_tries) break;

		count = 0;

		colIndex.push_back(pivot);
		remainingColIndex.erase(pivot);

		/// Normalizing constant
		Gamma   =   1.0/max;

		/// Generation of the column
		// a	=	A.col(colIndex.back());
		a	=	get_Col(start_Col+colIndex.back(), start_Row, n_Rows, mySpecies);
		/// Column of the residuum and the pivot row
		col =   a;
		for (int l=0; l<rank; ++l) {
			col =   col-v[l](colIndex.back())*u[l];
		}
		pivot   =   max_Abs_Vector(col, remainingRowIndex, unused_max);
		// std::cout << col << "\n";

		/// This randomization is needed if in the middle of the algorithm the columns happens to be exactly the linear combination of the previous columns.
		while (fabs(max)<tolerance && count < max_tries) {
			colIndex.pop_back();
			int new_colIndex	=	*remainingColIndex.begin();
			colIndex.push_back(new_colIndex);
			remainingColIndex.erase(new_colIndex);

			/// Generation of the column
			// a	=	A.col(new_colIndex);
			a	=	get_Col(start_Col+new_colIndex, start_Row, n_Rows, mySpecies);

			/// Column of the residuum and the pivot row
			col =   a;
			for (int l=0; l<rank; ++l) {
				col =   col-u[l](colIndex.back())*v[l];
			}
			pivot   =   max_Abs_Vector(col, remainingRowIndex, unused_max);
			++count;
		}

		if (count == max_tries) break;

		count = 0;

		rowIndex.push_back(pivot);
		remainingRowIndex.erase(pivot);

		/// New vectors
		u.push_back(Gamma*col);
		v.push_back(row);

		/// New approximation of matrix norm
		row_Squared_Norm    =   row.squaredNorm();
		row_Norm            =   sqrt(row_Squared_Norm);

		col_Squared_Norm    =   col.squaredNorm();
		col_Norm            =   sqrt(col_Squared_Norm);

		matrix_Norm         =   matrix_Norm +   Gamma*Gamma*row_Squared_Norm*col_Squared_Norm;

		for (int j=0; j<rank; ++j) {
			matrix_Norm     =   matrix_Norm +   2.0*(u[j].dot(u.back()))*(v[j].dot(v.back()));
		}
		++rank;
	} while (rank*row_Norm*col_Norm > fabs(max)*tolerance*matrix_Norm && rank < fmin(n_Rows, n_Cols)-1);

	U	=	Eigen::MatrixXd(n_Rows,rank);
	V	=	Eigen::MatrixXd(n_Cols,rank);

	#pragma omp parallel for
	for (int j=0; j<rank; ++j) {
		U.col(j)    =   u[j];
		V.col(j)    =   v[j];
	}
	// std::cout << "Computed rank: " << rank << "\n";
	return;
}

void grid::compute_RHS() {
	rhs	=	massFractionGradient + (massFraction-moleFraction)*pressureGradient/pressure;
	rhs	=	rhs*temperature/pressure*sqrt(temperature);
}

void grid::obtain_Inverse_Diffusion_Coefficients(const std::vector<species>& mySpecies) {
	D	=	get_Matrix(0, 0, nSpecies, nSpecies, mySpecies);
}

void grid::compute_Diagonal() {
	diagonal	=	U*(V.transpose()*moleFraction);
	invdiagonal	=	diagonal.cwiseInverse();
}

void grid::compute_Species_Velocity(const std::vector<species>& mySpecies, const Eigen::VectorXd& W, double tolerance) {
	this->speciesVelocity	=	Eigen::VectorXd::Zero(nSpecies);
	// std::cout << "\nComputing the compressed_Inverse_Diffusion_Coefficients...\n";
	compressed_Inverse_Diffusion_Coefficients(mySpecies, tolerance);
	compute_Diagonal();
	compute_RHS();
	Eigen::MatrixXd Utilde(nSpecies, rank+1), Vtilde(nSpecies, rank+1);
	Utilde << moleFraction.asDiagonal()*U, Eigen::VectorXd::Random(nSpecies);
	Vtilde << V, W;
	Eigen::VectorXd tempVelocity=	-invdiagonal.cwiseProduct(rhs);
	Utilde						=	invdiagonal.asDiagonal()*Utilde;
	Eigen::MatrixXd temp		=	Eigen::MatrixXd::Identity(rank+1,rank+1) - Vtilde.transpose()*Utilde;
	speciesVelocity				=	tempVelocity + Utilde*temp.fullPivLu().solve(Vtilde.transpose()*tempVelocity);
	speciesVelocity				=	speciesVelocity;
	// obtain_Inverse_Diffusion_Coefficients(mySpecies);
	// Eigen::MatrixXd Error	=	D-U*V.transpose();
	// std::cout << "Error is: " << Error.norm()/D.norm() << "\n";
	// std::cout << "Error is: " << Error.cwiseAbs().maxCoeff() << "\n";
	// std::cout << "Maximum element in D is: " << D.cwiseAbs().maxCoeff() << "\n";
	// std::cout << "Temperature is: " << temperature << "\n";
	// std::cout << "Grid is: " << x << "\n";
	// std::cout << "\nComputed the compressed_Inverse_Diffusion_Coefficients.\n";
	// std::cout << "\nComputing the right hand side...\n";
	// std::cout << "\nComputed the right hand side.\n";
}

void grid::compute_Exact_Species_Velocity(const std::vector<species>& mySpecies, const Eigen::VectorXd& W) {
	obtain_Inverse_Diffusion_Coefficients(mySpecies);
	Eigen::VectorXd temp	=	D*moleFraction;
	Eigen::MatrixXd M		=	moleFraction.asDiagonal()*D + Eigen::VectorXd::Random(nSpecies)*W.transpose();
	for (int j=0; j<nSpecies; ++j) {
		M(j,j)-=temp(j);
	}
	exactSpeciesVelocity	=	M.fullPivLu().solve(rhs);
}

void grid::compute_Species_Velocity_Iteratively(const std::vector<species>& mySpecies, const Eigen::VectorXd& W, double tolerance) {
	// obtain_Inverse_Diffusion_Coefficients(mySpecies);
	Eigen::VectorXd temp	=	D*moleFraction;
	Eigen::MatrixXd M		=	moleFraction.asDiagonal()*D + Eigen::VectorXd::Random(nSpecies)*W.transpose();
	for (int j=0; j<nSpecies; ++j) {
		M(j,j)-=temp(j);
	}

    // solve Ax = b using CG with matrix-free version:
	Eigen::BiCGSTAB<Eigen::MatrixXd > cg;
	cg.setTolerance(100*tolerance);
	cg.compute(M);
	iterativespeciesVelocity	=	cg.solve(rhs);
	this->iterativeerror	=	cg.error();
	this->nIterations		=	cg.iterations();
}
#endif /*__grid_hpp__*/