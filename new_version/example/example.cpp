#include <iostream>
#include "../FMDV.hpp"

struct species 
{
	//	Species name
	std::string name;
	//	Species mass
	double mass;
	//	Species diameter
	double diameter;
	//	Species Collision Energy
	double energy;
};

void readSpecies(std::string fileName, std::vector<species> &involved_species, int &N_species) 
{
	std::ifstream myfile;
	myfile.open(fileName.c_str(), std::ios::in);
	std::string temp_string;
	double temp_double;

	std::getline(myfile, temp_string);

	while(!myfile.eof()) 
	{
		species	newSpecies;
		myfile >> temp_string;
		myfile >> temp_double;
		newSpecies.mass     = temp_double;
		myfile >> temp_double;
		newSpecies.diameter	= temp_double;
		myfile >> temp_double;
		newSpecies.energy   = temp_double;
		involved_species.push_back(newSpecies);
	}

	myfile.close();

	this->nSpecies = involved_species.size();
	this->diaMin   = mySpecies[0].diameter;
	this->diaMax   = mySpecies[0].diameter;
	this->massMin  = mySpecies[0].mass;
	this->massMax  = mySpecies[0].mass;

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
	for (int j=0; j<nSpecies; ++j) 
    {
		mySpecies[j].mass		=	mySpecies[j].mass/massMax;
		mySpecies[j].diameter	=	0.5*mySpecies[j].diameter/diaMax;
		W(j)					=	mySpecies[j].mass;
	}
}

int main(int argc, char* argv[])
{
    srand(time(NULL));

    int N_grid = 100;
    // File names to read from:
	std::string filename[] = {"data/table_224.txt","data/table_448.txt","data/table_896.txt","data/table_1792.txt","data/table_3585.txt","data/table_7171.txt"};

	// Iterating over each of the files:
	for (int k = 0; k < 6; ++k) 
	{
        FMDV F(N);

    F.setPressure(double(rand()) / RAND_MAX);
    F.setPressureGradient(double(rand()) / RAND_MAX);
    F.setDensity(double(rand()) / RAND_MAX);
    F.setTemperature(double(rand()) / RAND_MAX);
    F.setTemperatureGradient(double(rand()) / RAND_MAX);

    Eigen::VectorXd X = Eigen::VectorXd(N).array() + 1;
    // Normalizing:
    X                 = X / X.sum();
    Eigen::VectorXd W = Eigen::VectorXd(N).array() + 1;

    F.setMolecularWeights(W);
    F.setDiameters(Eigen::VectorXd::Random(N).array() + 1);
    F.setThermalDiffusivities(Eigen::VectorXd::Random(N).array() + 1);
    F.setCollisionEnergies(Eigen::VectorXd::Random(N).array() + 1);
    F.setMoleFraction(X);
    F.setMoleFractionGradient(Eigen::VectorXd::Random(N).array() + 1);
    F.setBodyForces(Eigen::VectorXd::Random(N).array() + 1);

    Eigen::MatrixXd V    = F.getInverseDiffusionCoefficientMatrix();
    Eigen::VectorXd v_sp = F.computeSpeciesVelocities(1e-12);

    // (diag(VX) - diag(X)V - S * Wt)z = (b - alpha * s)
    Eigen::VectorXd S       = F.getRandomVector();
    Eigen::MatrixXd M       = Eigen::MatrixXd((V * X).asDiagonal()) - X.asDiagonal() * V - S * W.transpose();
    Eigen::VectorXd z_exact = M.fullPivLu().solve(F.getRHS());

    double prefactor  = F.getTemperatureGradient() / (F.getDensity() * F.getTemperature());
    Eigen::VectorXd D = F.getThermalDiffusivities();
    Eigen::VectorXd Y = F.getMassFraction();

    Eigen::VectorXd v_sp_exact(N);
    for(int i = 0; i < N; i++)
    {
        if(X(i) == 0)
        {
            v_sp_exact(i) = 0;
        }

        else
        {
            v_sp_exact(i) = z_exact(i) / X(i) - prefactor * D(i) / Y(i);
        }
    }

    std::cout << "Error:" << (v_sp - v_sp_exact).cwiseAbs().maxCoeff() << std::endl;
}
