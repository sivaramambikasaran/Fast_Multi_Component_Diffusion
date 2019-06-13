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
	// Number of species in the domain
	int N_species;

	// Grid size
	double dx;
	// Maximum diameter and Maximum Mass
	double diaMax, massMax;

	// Vector used to store grid content:
	std::vector<grid_point> grid;
	// Species that are present:
	std::vector<species> species_present;

public:

	// Constructor
	domain(int N_grid)
	{	
		this->N_grid = N_grid;
		this->dx     = 2 / N_grid;
	};

	// Destructor
	~domain(){};

	// This function is used to get the information about the species in the domain:
	void readSpecies(std::string file_name); 
	// Generates the grid:
	void generateGrid(int N_grid);
	// Generate species profile: generates profile for mass 
	// and mole fraction and their gradients in addition to body forces 
	void generateSpeciesProfile();
	// Generate temperature and pressure profile
	void generateTempPressureProfile();

	// Compute species velocity
	void computeSpeciesVelocity(double tolerance);
	// Compute exact species velocity
	void computeExactSpeciesVelocity();
	// Compute species velocity with iterative solver
	void computeSpeciesVelocityIteratively(double tolerance);

	// Compute error
	void computeError();
	// Compute iterative error
	void computeIterative_Error();
};

void domain::generateGrid(int N_grid) 
{

	for (int j=0; j<N_grid; ++j) 
	{
		grid_point new_grid_point(left+(j+0.5)*dx, N_species);
		grid.push_back(new_grid_point);
	}

	// Getting the mole, mass fraction and gradient profiles:
	generateSpeciesProfile();
	// Getting the temperature, pressure profiles and their gradients:
	generateTemperaturePressureProfile();
}

// We will be using this method to read the data from our data files:
void domain::readSpecies(std::string file_name) 
{
    std::ifstream myfile;
    myfile.open(file_name.c_str(), std::ios::in);
    std::string temp_string;
    double temp_double;

    std::getline(myfile, temp_string);

    while(!myfile.eof()) 
    {
        species new_species;
        myfile >> temp_string;
        new_species.name = temp_string;
        myfile >> temp_double;
        new_species.mass = temp_double;
        myfile >> temp_double;
        new_species.diameter = temp_double;
        myfile >> temp_double;
        new_species.energy = temp_double;
        // Our data files contain data for the above attributes
        // However, they do not contain information about the
        // thermal diffusivities body forces. For the sake of
        // this example, we assign them using random values:
        new_species.thermal_diffusivity = double(rand()) / RAND_MAX;
        species_present.push_back(new_species);
    }
    
    myfile.close();

    this->N_species = species_present.size();
    this->dia_max   = species_present[0].diameter;
    this->mass_max  = species_present[0].mass;

    // Finding the maximum mass and diameter:
    for(int j = 1; j < N_species; ++j) 
    {
        if (species_present[j].mass > mass_max) 
        {
            mass_max = species_present[j].mass;
        }
        if (species_present[j].diameter > dia_max) 
        {
            dia_max = species_present[j].diameter;
        }
    }

	// Normalizing the quantities such that the maximum of the mass and diameter is 1:
    #pragma omp parallel for
    for(int j = 0; j < N_species; ++j) 
    {
        species_present[j].mass     = species_present[j].mass/mass_max;
        species_present[j].diameter = species_present[j].diameter/dia_max;
    }
}

// Generates a random species profile:
void domain::generateSpeciesProfile() 
{
	// Obtaining mole fraction matrix of size N_species by N_grid
	Eigen::MatrixXd mole_fraction =	Eigen::MatrixXd::Ones(N_species, N_grid) + Eigen::MatrixXd::Random(N_species, N_grid);
	Eigen::VectorXd sum           = mole_fraction.colwise().sum();

	#pragma omp parallel for
	for (int j = 0; j < N_grid; ++j) 
	{
		mole_fraction.col(j) = mole_fraction.col(j) / sum(j);
	}

	// Taking central difference assuming periodic BCs:
	Eigen::MatrixXd mole_fraction_gradient                  = Eigen::MatrixXd::Zero(N_species, N_grid);
	mole_fraction_gradient.block(0, 1, N_species, N_grid-2) = mole_fraction.block(0, 2, N_species, N_grid-2) - mole_fraction.block(0, 0, N_species, N_grid-2);
	#pragma omp parallel for
	for (int j = 0; j < N_species; ++j) 
	{
		mole_fraction_gradient(j, 0)          = mole_fraction(j, 1) - mole_fraction(j, N_grid - 1);
		mole_fraction_gradient(j, N_grid - 1) = mole_fraction(j, 0) - mole_fraction(j, N_grid - 2);
	}

	mole_fraction_gradient = 0.5 / dx * mole_fraction_gradient;

	// Obtaining mass fraction matrix of size N_species by N_grid
	Eigen::MatrixXd mass_fraction = Eigen::MatrixXd(N_species, N_grid);
	#pragma omp parallel for
	for(int j = 0; j < N_species; ++j) 
	{
		mass_fraction.row(j) = species_present[j].mass * mole_fraction.row(j);
	}
	
	sum	= mass_fraction.colwise().sum();
	#pragma omp parallel for
	for(int j = 0; j < N_grid; ++j) 
	{
		mass_fraction.col(j) = mass_fraction.col(j) / sum(j);
	}

	// Declaring matrix for body forces:
	Eigen::MatrixXd body_forces = Eigen::MatrixXd::Random(N_species, N_grid);
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].mole_fraction          = mole_fraction.col(j);
		grid[j].mole_fraction_gradient = mole_fraction_gradient.col(j);
		grid[j].mass_fraction          = mass_fraction.col(j);
		grid[j].body_forces            = body_forces.col(j);
	}
}

// Generates temperature and pressure profiles
void domain::generateTemperaturePressureProfile() 
{
	double tempMean				=	2000;
	double tempScale			=	500;
	Eigen::VectorXd temp		=	Eigen::VectorXd::Ones(N_nodes)+Eigen::VectorXd::Random(N_nodes);
	Eigen::VectorXd temperature	=	tempMean*Eigen::VectorXd::Ones(N_grid) + tempScale*(T*temp/temp.sum());

	Eigen::VectorXd pressureGradient	=	Eigen::VectorXd::Zero(N_grid);
	pressureGradient.segment(1,N_grid-2)	=	pressure.segment(2,N_grid-2)-pressure.segment(0,N_grid-2);
	pressureGradient	=	(0.5/dx)*pressureGradient;

	double pressureMean			=	101325;
	Eigen::VectorXd pressure	=	pressureMean*Eigen::VectorXd::Ones(N_grid);

	Eigen::VectorXd pressureGradient	=	Eigen::VectorXd::Zero(N_grid);
	pressureGradient.segment(1,N_grid-2)	=	pressure.segment(2,N_grid-2)-pressure.segment(0,N_grid-2);
	pressureGradient	=	(0.5/dx)*pressureGradient;

	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].temperature		=	temperature(j);
		grid[j].pressure			=	pressure(j);
		grid[j].pressureGradient	=	pressureGradient(j);
		grid[j].pressureGradient	=	pressureGradient(j);
	}
}

// Computes the species velocity at every grid point
void domain::compute_Species_Velocity(double tolerance) {
	static const double preFactor = 2.0/(3.0*AVAGADRO*diaMax*diaMax*sqrt(massMax))*(R/PI)*sqrt(R/PI);
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		grid[j].compute_Species_Velocity(species_present, W, tolerance);
		grid[j].speciesVelocity	=	preFactor*grid[j].speciesVelocity;
	}
}

// Computes the exact species velocity at every grid point
void domain::compute_Exact_Species_Velocity() {
	static const double preFactor	=	1;//2.0/(3.0*AVAGADRO*diaMax*diaMax*sqrt(massMax))*(R/PI)*sqrt(R/PI);
	// std::cout << "Prefactor is: " << preFactor << "\n";
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		grid[j].compute_Exact_Species_Velocity(species_present, W);
		grid[j].exactSpeciesVelocity	=	preFactor*grid[j].exactSpeciesVelocity;
	}
}

// Solve the system iteratively using CG
void domain::compute_Species_Velocity_Iteratively(double tolerance, double& iterativeerror) {
	static const double preFactor	=	1;//2.0/(3.0*AVAGADRO*diaMax*diaMax*sqrt(massMax))*(R/PI)*sqrt(R/PI);
	// std::cout << "Prefactor is: " << preFactor << "\n";
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) {
		grid[j].compute_Species_Velocity_Iteratively(species_present, W, tolerance);
		grid[j].iterativespeciesVelocity	=	preFactor*grid[j].iterativespeciesVelocity;
	}
	std::cout << "Number of iterations for the iterative solver is: " << grid[N_grid/2].nIterations << "\n";
	iterativeerror	=	grid[N_grid/2].iterativeerror;
}

// Computes the error
void domain::computeError() 
{
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].error	=	(grid[j].exactSpeciesVelocity-grid[j].speciesVelocity).norm()/(grid[j].exactSpeciesVelocity.norm());
	}
}

#endif /*__domain_hpp__*/
