#ifndef __domain_hpp__
#define __domain_hpp__


#include<cassert>
#include<vector>
#include<Eigen/Dense>
#include<Eigen/IterativeLinearSolvers>
#include"grid.hpp"

class domain
{
public:
	//	Constructor
	domain();

	//	Destructor
	~domain();

	//	Grids
	std::vector<grid> grids;

	//	Number of grid points in the domain
	int nGrid;

	//	Grid size
	double dx;

	//	Generates the grids
	void generate_Grid(int nGrid, double left, double right);

	//	Species
	std::vector<species> mySpecies;

	//	Vector of mass of the species
	Eigen::VectorXd W;

	//	Minimum, Maximum diameter and Minimum and Maximum Mass
	double diaMin, diaMax, massMin, massMax;

	//	Read the properties from the file
	void read_Species(std::string fileName);

	//	Number of species in the domain
	int nSpecies;

	//	Generate species profile
	void generate_Species_Profile(int nNodes);

	//	Generate temperature and pressure profile
	void generate_Temp_Pressure_Profile(int nNodes);

	//	Compute species velocity
	void compute_Species_Velocity(double tolerance);

	//	Compute exact species velocity
	void compute_Exact_Species_Velocity();

	//	Compute species velocity with iterative solver
	void compute_Species_Velocity_Iteratively(double tolerance, double& iterativeerror);

	//	Compute error
	void compute_Error();

	//	Compute iterative error
	void compute_Iterative_Error();
};

domain::domain() {};
domain::~domain() {};

void domain::generate_Grid(int nGrid, double left, double right) 
{
	this->nGrid		=	nGrid;
	assert(nGrid>1);
	this->dx		=	(right-left)/nGrid;

	for (int j=0; j<nGrid; ++j) {
		grid myGrid(left+(j+0.5)*dx, nSpecies);
		grids.push_back(myGrid);
	}
	// std::cout << "\nGenerating the species profile...\n";
	generate_Species_Profile(20);
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

	W	=	Eigen::VectorXd(nSpecies);

	#pragma omp parallel for
	for (int j=0; j<nSpecies; ++j) {
		mySpecies[j].mass		=	mySpecies[j].mass/massMax;
		mySpecies[j].diameter	=	0.5*mySpecies[j].diameter/diaMax;
		W(j)					=	mySpecies[j].mass;
	}
}

//	Generates species profile
void domain::generate_Species_Profile(int nNodes) 
{
	Eigen::MatrixXd T	=	Eigen::MatrixXd(nGrid, nNodes);
	T.col(0)=	Eigen::VectorXd::Ones(nGrid);

	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		T(j,1)	=	grids[j].x;
	}

	#pragma omp parallel for
	for (int j=0; j<nNodes-2; ++j) {
		for (int k=0; k<nGrid; ++k) {
			T(k,j+2)	=	2*grids[k].x*T(k,j+1)-T(k,j);
		}
	}
	T+=Eigen::MatrixXd::Ones(nGrid,nNodes);
	Eigen::MatrixXd moleFraction	=	(0.5*T*(Eigen::MatrixXd::Ones(nNodes,nSpecies)+Eigen::MatrixXd::Random(nNodes,nSpecies))).transpose();
	//	moleFraction	-	Matrix of species mole fraction at each grid point; Size is nSpecies by nGrid
	Eigen::VectorXd sum	=	moleFraction.colwise().sum();

	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		moleFraction.col(j)	=	moleFraction.col(j)/sum(j);
	}

	Eigen::MatrixXd moleFractionGradient	=	Eigen::MatrixXd::Zero(nSpecies, nGrid);
	moleFractionGradient.block(0,1,nSpecies,nGrid-2)	=	moleFraction.block(0,2,nSpecies,nGrid-2)-moleFraction.block(0,0,nSpecies,nGrid-2);
	moleFractionGradient	=	0.5/dx*moleFractionGradient;

	//	Obtaining mass fractions and mass fraction gradients of size nSpecies by nGrid
	Eigen::MatrixXd massFraction		=	Eigen::MatrixXd(nSpecies, nGrid);
	Eigen::MatrixXd massFractionGradient=	Eigen::MatrixXd(nSpecies, nGrid);

	#pragma omp parallel for
	for (int j=0; j<nSpecies; ++j) {
		massFraction.row(j)	=	mySpecies[j].mass*moleFraction.row(j);
	}
	sum	=	massFraction.colwise().sum();
	
	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		massFraction.col(j)	=	massFraction.col(j)/sum(j);
	}

	massFractionGradient	=	Eigen::MatrixXd::Zero(nSpecies, nGrid);
	massFractionGradient.block(0,1,nSpecies,nGrid-2)	=	massFraction.block(0,2,nSpecies,nGrid-2)-massFraction.block(0,0,nSpecies,nGrid-2);
	massFractionGradient	=	0.5/dx*massFractionGradient;

	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		grids[j].moleFraction			=	moleFraction.col(j);
		grids[j].moleFractionGradient	=	moleFractionGradient.col(j);
		grids[j].massFraction			=	massFraction.col(j);
		grids[j].massFractionGradient	=	massFractionGradient.col(j);
		grids[j].mass					=	0.0;
		for (int k=0; k<nSpecies; ++k) {
			grids[j].mass+=grids[j].moleFraction(k)*mySpecies[k].mass;
		}
	}
}

//	Generates temperature and pressure profile
void domain::generate_Temp_Pressure_Profile(int nNodes) {
	Eigen::MatrixXd T	=	Eigen::MatrixXd(nGrid, nNodes);
	T.col(0)=	Eigen::VectorXd::Ones(nGrid);
	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		T(j,1)	=	grids[j].x;
	}
	#pragma omp parallel for
	for (int j=0; j<nNodes-2; ++j) {
		for (int k=0; k<nGrid; ++k) {
			T(k,j+2)	=	2*grids[k].x*T(k,j+1)-T(k,j);
		}
	}
	double tempMean				=	2000;
	double tempScale			=	500;
	Eigen::VectorXd temp		=	Eigen::VectorXd::Ones(nNodes)+Eigen::VectorXd::Random(nNodes);
	Eigen::VectorXd temperature	=	tempMean*Eigen::VectorXd::Ones(nGrid) + tempScale*(T*temp/temp.sum());

	double pressureMean			=	101325;
	Eigen::VectorXd pressure	=	pressureMean*Eigen::VectorXd::Ones(nGrid);

	Eigen::VectorXd pressureGradient	=	Eigen::VectorXd::Zero(nGrid);
	pressureGradient.segment(1,nGrid-2)	=	pressure.segment(2,nGrid-2)-pressure.segment(0,nGrid-2);
	pressureGradient	=	(0.5/dx)*pressureGradient;

	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		grids[j].temperature		=	temperature(j);
		grids[j].pressure			=	pressure(j);
		grids[j].pressureGradient	=	pressureGradient(j);
	}

}

//	Computes the species velocity at every grid point
void domain::compute_Species_Velocity(double tolerance) {
	static const double preFactor	=	2.0/(3.0*Avagadro*diaMax*diaMax*sqrt(massMax))*(R/PI)*sqrt(R/PI);
	// std::cout << "Prefactor is: " << preFactor << "\n";
	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		grids[j].compute_Species_Velocity(mySpecies, W, tolerance);
		grids[j].speciesVelocity	=	preFactor*grids[j].speciesVelocity;
	}
	// std::cout << "Species velocity is: " << grids[nGrid/2].speciesVelocity.transpose() << "\n\n";
	// std::cout << "Molefraction gradient is: " << grids[nGrid-2].moleFractionGradient.transpose() << "\n\n";
}

//	Computes the exact species velocity at every grid point
void domain::compute_Exact_Species_Velocity() {
	static const double preFactor	=	2.0/(3.0*Avagadro*diaMax*diaMax*sqrt(massMax))*(R/PI)*sqrt(R/PI);
	// std::cout << "Prefactor is: " << preFactor << "\n";
	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		grids[j].compute_Exact_Species_Velocity(mySpecies, W);
		grids[j].exactSpeciesVelocity	=	preFactor*grids[j].exactSpeciesVelocity;
	}
}

//	Computes the error
void domain::compute_Error() {
	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		grids[j].error	=	(grids[j].exactSpeciesVelocity-grids[j].speciesVelocity).norm()/(grids[j].exactSpeciesVelocity.norm());
		// std::cout << grids[j].error << "\n";
	}
	std::cout << "Rank of the inverse diffusion matrix is: " << grids[nGrid/2].rank << "\n";
}

//	Solve the system iteratively using CG
void domain::compute_Species_Velocity_Iteratively(double tolerance, double& iterativeerror) {
	static const double preFactor	=	2.0/(3.0*Avagadro*diaMax*diaMax*sqrt(massMax))*(R/PI)*sqrt(R/PI);
	// std::cout << "Prefactor is: " << preFactor << "\n";
	#pragma omp parallel for
	for (int j=0; j<nGrid; ++j) {
		grids[j].compute_Species_Velocity_Iteratively(mySpecies, W, tolerance);
		grids[j].iterativespeciesVelocity	=	preFactor*grids[j].iterativespeciesVelocity;
	}
	std::cout << "Number of iterations for the iterative solver is: " << grids[nGrid/2].nIterations << "\n";
	iterativeerror	=	grids[nGrid/2].iterativeerror;
}


#endif /*__domain_hpp__*/