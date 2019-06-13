#ifndef __domain_hpp__
#define __domain_hpp__

#include <vector>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

const double R        = 8.3144598;
const double AVAGADRO = 6.022140857e23;
const double kB       = R / AVAGADRO;
const double PI       = 3.141592653589793238;

class domain
{
	// Number of grid points in the domain
	int N_grid;
	// Vector used to store grid content:
	std::vector<grid_point> grid;

public:
	// Constructor
	domain()
	{

	};
	// Destructor
	~domain(){};

	// Grid size
	double dx;
	// Generates the grid:
	void generateGrid(double start, double end, int N_grid);

	// Species
	std::vector<species> species_present;

	// Minimum, Maximum diameter and Minimum and Maximum Mass
	double diaMin, diaMax, massMin, massMax;

	// Number of species in the domain
	int nSpecies;

	// Generate species profile
	void generateSpeciesProfile(int N_nodes);
	// Generate temperature and pressure profile
	void generateTempPressureProfile(int N_nodes);

	// Compute species velocity
	void computeSpeciesVelocity(double tolerance);
	// Compute exact species velocity
	void computeExactSpeciesVelocity();
	// Compute species velocity with iterative solver
	void computeSpeciesVelocityIteratively(double tolerance);

	// Compute error
	void compute_Error();
	// Compute iterative error
	void compute_Iterative_Error();
};

void domain::generate_Grid(int N_grid, double left, double right) 
{
	this->N_grid		=	N_grid;
	assert(N_grid>1);
	this->dx		=	(right-left)/N_grid;

	for (int j=0; j<N_grid; ++j) {
		grid myGrid(left+(j+0.5)*dx, nSpecies);
		grid.push_back(myGrid);
	}
	// std::cout << "\nGenerating the species profile...\n";
	generateSpeciesProfile(20);
	// std::cout << "\nGenerated the species profile.\n";

	// std::cout << "\nGenerating the temperature profile...\n";
	generate_Temp_Pressure_Profile(20);
	// std::cout << "\nGenerated the temperature profile.\n";
}

void domain::read_Species(std::string fileName) 
{
	std::ifstream myfile;
	myfile.open(fileName.c_str(), std::ios::in);
	std::string tempString;
	double tempDouble;

	std::getline(myfile, tempString);

	while(!myfile.eof()) 
	{
		species	newSpecies;
		myfile >> tempString;
		newSpecies.name		=	tempString;
		myfile >> tempDouble;
		newSpecies.mass		=	tempDouble;
		myfile >> tempDouble;
		newSpecies.diameter	=	tempDouble;
		myfile >> tempDouble;
		newSpecies.energy	=	tempDouble;
		mySpecies.push_back(newSpecies);
	}
	myfile.close();

	this->nSpecies	=	mySpecies.size();
	this->diaMin	=	mySpecies[0].diameter;
	this->diaMax	=	mySpecies[0].diameter;
	this->massMin	=	mySpecies[0].mass;
	this->massMax	=	mySpecies[0].mass;

	for (int j=1; j<nSpecies; ++j) {
		if (mySpecies[j].mass < massMin) {
			massMin	=	mySpecies[j].mass;
		}
		if (mySpecies[j].mass > massMax) {
			massMax	=	mySpecies[j].mass;
		}
		if (mySpecies[j].diameter < diaMin) {
			diaMin	=	mySpecies[j].diameter;
		}
		if (mySpecies[j].diameter > diaMax) {
			diaMax	=	mySpecies[j].diameter;
		}
	}

	W	=	Eigen::VectorXd::Random(nSpecies).cwiseAbs();

	#pragma omp parallel for
	for (int j=0; j<nSpecies; ++j) {
		mySpecies[j].mass		=	W(j); //mySpecies[j].mass/massMax;
		mySpecies[j].diameter	=	(double)(rand()) / RAND_MAX;
		// W(j)					=	mySpecies[j].mass;
	}
}

// Generates species profile
void domain::generateSpeciesProfile(int N_nodes) 
{
	Eigen::MatrixXd T = Eigen::MatrixXd(N_grid, N_nodes);
	T.col(0)          = Eigen::VectorXd::Ones(N_grid);

	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		T(j,1) = grid[j].x;
	}

	#pragma omp parallel for
	for (int j=0; j<N_nodes-2; ++j) 
	{
		for (int k=0; k<N_grid; ++k) 
		{
			T(k,j+2)	=	2*grid[k].x*T(k,j+1)-T(k,j);
		}
	}
	T+=Eigen::MatrixXd::Ones(N_grid,N_nodes);
	Eigen::MatrixXd moleFraction	=	(0.5*T*(Eigen::MatrixXd::Ones(N_nodes,nSpecies)+Eigen::MatrixXd::Random(N_nodes,nSpecies))).transpose();
	// moleFraction	-	Matrix of species mole fraction at each grid point; Size is nSpecies by N_grid
	Eigen::VectorXd sum	=	moleFraction.colwise().sum();

	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		moleFraction.col(j)	=	moleFraction.col(j)/sum(j);
	}

	Eigen::MatrixXd moleFractionGradient	=	Eigen::MatrixXd::Zero(nSpecies, N_grid);
	moleFractionGradient.block(0,1,nSpecies,N_grid-2)	=	moleFraction.block(0,2,nSpecies,N_grid-2)-moleFraction.block(0,0,nSpecies,N_grid-2);
	moleFractionGradient	=	0.5/dx*moleFractionGradient;

	// Obtaining mass fractions and mass fraction gradients of size nSpecies by N_grid
	Eigen::MatrixXd massFraction		=	Eigen::MatrixXd(nSpecies, N_grid);
	Eigen::MatrixXd massFractionGradient=	Eigen::MatrixXd(nSpecies, N_grid);

	#pragma omp parallel for
	for (int j=0; j<nSpecies; ++j) {
		massFraction.row(j)	=	mySpecies[j].mass*moleFraction.row(j);
	}
	sum	=	massFraction.colwise().sum();
	
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		massFraction.col(j)	=	massFraction.col(j)/sum(j);
	}

	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].moleFraction			=	moleFraction.col(j);
		grid[j].moleFractionGradient	=	moleFractionGradient.col(j);
		grid[j].massFraction			=	massFraction.col(j);
		grid[j].massFractionGradient	=	massFractionGradient.col(j);
		grid[j].mass					=	0.0;
		for (int k=0; k<nSpecies; ++k) {
			grid[j].mass+=grid[j].moleFraction(k)*mySpecies[k].mass;
		}
	}
}

// Generates temperature and pressure profile
void domain::generate_Temp_Pressure_Profile(int N_nodes) {
	Eigen::MatrixXd T	=	Eigen::MatrixXd(N_grid, N_nodes);
	T.col(0)=	Eigen::VectorXd::Ones(N_grid);
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		T(j,1)	=	grid[j].x;
	}
	#pragma omp parallel for
	for (int j=0; j<N_nodes-2; ++j) {
		for (int k=0; k<N_grid; ++k) {
			T(k,j+2)	=	2*grid[k].x*T(k,j+1)-T(k,j);
		}
	}
	double tempMean				=	2000;
	double tempScale			=	500;
	Eigen::VectorXd temp		=	Eigen::VectorXd::Ones(N_nodes)+Eigen::VectorXd::Random(N_nodes);
	Eigen::VectorXd temperature	=	tempMean*Eigen::VectorXd::Ones(N_grid) + tempScale*(T*temp/temp.sum());

	double pressureMean			=	101325;
	Eigen::VectorXd pressure	=	pressureMean*Eigen::VectorXd::Ones(N_grid);

	Eigen::VectorXd pressureGradient	=	Eigen::VectorXd::Zero(N_grid);
	pressureGradient.segment(1,N_grid-2)	=	pressure.segment(2,N_grid-2)-pressure.segment(0,N_grid-2);
	pressureGradient	=	(0.5/dx)*pressureGradient;

	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		grid[j].temperature		=	temperature(j);
		grid[j].pressure			=	pressure(j);
		grid[j].pressureGradient	=	pressureGradient(j);
	}

}

// Computes the species velocity at every grid point
void domain::compute_Species_Velocity(double tolerance) {
	static const double preFactor	=	1;//2.0/(3.0*AVAGADRO*diaMax*diaMax*sqrt(massMax))*(R/PI)*sqrt(R/PI);
	// std::cout << "Prefactor is: " << preFactor << "\n";
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		grid[j].compute_Species_Velocity(mySpecies, W, tolerance);
		grid[j].speciesVelocity	=	preFactor*grid[j].speciesVelocity;
	}
	// std::cout << "Species velocity is: " << grid[N_grid/2].speciesVelocity.transpose() << "\n\n";
	// std::cout << "Molefraction gradient is: " << grid[N_grid-2].moleFractionGradient.transpose() << "\n\n";
}

// Computes the exact species velocity at every grid point
void domain::compute_Exact_Species_Velocity() {
	static const double preFactor	=	1;//2.0/(3.0*AVAGADRO*diaMax*diaMax*sqrt(massMax))*(R/PI)*sqrt(R/PI);
	// std::cout << "Prefactor is: " << preFactor << "\n";
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		grid[j].compute_Exact_Species_Velocity(mySpecies, W);
		grid[j].exactSpeciesVelocity	=	preFactor*grid[j].exactSpeciesVelocity;
	}
}

// Computes the error
void domain::compute_Error() {
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		grid[j].error	=	(grid[j].exactSpeciesVelocity-grid[j].speciesVelocity).norm()/(grid[j].exactSpeciesVelocity.norm());
		// std::cout << grid[j].error << "\n";
	}
	std::cout << "Rank of the inverse diffusion matrix is: " << grid[N_grid/2].rank << "\n";
}

// Solve the system iteratively using CG
void domain::compute_Species_Velocity_Iteratively(double tolerance, double& iterativeerror) {
	static const double preFactor	=	1;//2.0/(3.0*AVAGADRO*diaMax*diaMax*sqrt(massMax))*(R/PI)*sqrt(R/PI);
	// std::cout << "Prefactor is: " << preFactor << "\n";
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		grid[j].compute_Species_Velocity_Iteratively(mySpecies, W, tolerance);
		grid[j].iterativespeciesVelocity	=	preFactor*grid[j].iterativespeciesVelocity;
	}
	std::cout << "Number of iterations for the iterative solver is: " << grid[N_grid/2].nIterations << "\n";
	iterativeerror	=	grid[N_grid/2].iterativeerror;
}

#endif /*__domain_hpp__*/
