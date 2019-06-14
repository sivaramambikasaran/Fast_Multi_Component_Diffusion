#ifndef __domain_hpp__
#define __domain_hpp__

#include <vector>
#include "grid_point.hpp"
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
	double dia_max, mass_max;

	// Vector used to store gridpoints:
	std::vector<grid_point> grid;
	// Species that are present:
	std::vector<species> species_present;

	// The following vectors store the species properties:
	Eigen::VectorXd molecular_mass, diameters, collision_energies, thermal_diffusivities;

	// Generates the grid:
	void generateGrid();
	// Generate species profile: generates profile for mass 
	// and mole fraction and their gradients in addition to body forces 
	void generateSpeciesProfile();
	// Generate temperature and pressure profile
	void generateTempPressureProfile();

public:

	// Constructor
	domain(int N_grid)
	{	
		this->N_grid = N_grid;
		// We assume standard interval of [-1, 1]
		this->dx     = 2 / N_grid;

		// Generates the grid:
		this->generateGrid();
	};

	// Destructor
	~domain(){};

	// This function is used to get the information about the species in the domain:
	void readSpecies(std::string file_name); 

	// Compute species velocity
	void computeSpeciesVelocityFast(double tolerance);
	// Compute exact species velocity
	void computeSpeciesVelocityExact();
	// Compute species velocity with iterative solver
	void computeSpeciesVelocityIterative(double tolerance);

	// Computes and returns the error
	void getError();
};

void domain::generateGrid() 
{

	for (int j=0; j<N_grid; ++j) 
	{
		grid_point new_grid_point(-1 + (j+0.5) * dx);
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

	// Setting these values obtained to Eigen vectors:
	this->molecular_mass        = Eigen::VectorXd::Zero(N_species);
	this->diameters             = Eigen::VectorXd::Zero(N_species);
	this->collision_energies    = Eigen::VectorXd::Zero(N_species);
	this->thermal_diffusivities = Eigen::VectorXd::Zero(N_species);

	#pragma omp parallel for
    for(int j = 0; j < N_species; ++j) 
    {
		this->molecular_mass(j)        = species_present[j].mass;
		this->diameters(j)             = species_present[j].diameter;
		this->collision_energies(j)    = species_present[j].energy;
		this->thermal_diffusivities(j) = species_present[j].thermal_diffusivity;
	}

	#pragma omp parallel for
	for (int j = 0; j < N_grid; ++j) 
	{
		grid[j].FS.molecular_mass        = this->molecular_mass;
		grid[j].FS.diameters             = this->diameters;
		grid[j].FS.collision_energies    = this->collision_energies;
		grid[j].FS.thermal_diffusivities = this->thermal_diffusivities;
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

	// Declaring matrix with random values for body forces:
	Eigen::MatrixXd body_forces = Eigen::MatrixXd::Random(N_species, N_grid);
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].FS.mole_fraction          = mole_fraction.col(j);
		grid[j].FS.mole_fraction_gradient = mole_fraction_gradient.col(j);
		grid[j].FS.mass_fraction          = mass_fraction.col(j);
		grid[j].FS.body_forces            = body_forces.col(j);
	}
}

// Generates temperature and pressure profiles
void domain::generateTemperaturePressureProfile() 
{
	// Initializing density with a constant value of 1:
	Eigen::VectorXd density = Eigen::VectorXd::Ones(N_grid);

	// Initializing with random temperature variation:
	double temp_mean            = 2000;
	double temp_scale           = 500;
	Eigen::VectorXd temperature	= temp_mean * Eigen::VectorXd::Ones(N_grid) + temp_scale * (Eigen::VectorXd::Random(N_grid));

	// Taking Central Difference:
	Eigen::VectorXd temperature_gradient     = Eigen::VectorXd::Zero(N_grid);
	temperature_gradient.segment(1,N_grid-2) = temperature.segment(2,N_grid-2) - temperature.segment(0,N_grid-2);
	temperature_gradient(0)                  = temperature(1) - temperature(N-1);
	temperature_gradient(N-1)                = temperature(0) - temperature(N-2);
	temperature_gradient                     = (0.5/dx) * temperature_gradient;

	// Assuming constant pressure throughout:
	double pressure_mean     = 101325;
	Eigen::VectorXd pressure = pressure_mean * Eigen::VectorXd::Ones(N_grid);

	// Due to constant pressure, pressure gradients would be zero:
	Eigen::VectorXd pressure_gradient = Eigen::VectorXd::Zero(N_grid);

	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].FS.density              = density(j);
		grid[j].FS.temperature          = temperature(j);
		grid[j].FS.pressure             = pressure(j);
		grid[j].FS.temperature_gradient = temperature_gradient(j);
		grid[j].FS.pressure_gradient    = pressure_gradient(j);
	}
}

// Computes the species velocity at every grid point
void domain::computeSpeciesVelocityFast(double tolerance) 
{
	static const double prefactor = 2.0/(3.0*AVAGADRO * dia_max * dia_max * sqrt(mass_max)) * (R/PI) * sqrt(R/PI);
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].computeSpeciesVelocityFast(tolerance);
		grid[j].species_velocity = prefactor * grid[j].species_velocity;
	}
}

// Computes the exact species velocity at every grid point
void domain::computeSpeciesVelocityExact() 
{
	static const double prefactor = 2.0/(3.0*AVAGADRO * dia_max * dia_max * sqrt(mass_max)) * (R/PI) * sqrt(R/PI);
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].computeSpeciesVelocityExact();
		grid[j].exact_species_velocity = prefactor * grid[j].exact_species_velocity;
	}
}

// Solve the system iteratively using CG
void domain::computeSpeciesVelocityIterative(double tolerance) 
{
	static const double prefactor = 2.0/(3.0*AVAGADRO * dia_max * dia_max * sqrt(mass_max)) * (R/PI) * sqrt(R/PI);
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].computeSpeciesVelocityIterative(tolerance);
		grid[j].iterative_species_velocity = prefactor * grid[j].iterative_species_velocity;
	}
}

// Computes and returns the error:
double domain::getError() 
{
	double error;
	#pragma omp parallel for
	for (int j=0; j<N_grid; ++j) 
	{
		grid[j].error = (grid[j].exact_species_velocity-grid[j].species_velocity).norm() / (grid[j].exact_species_velocity.norm());
		error        += grid[j].error;
	}

	return error / N_grid;
}

#endif /*__domain_hpp__*/
