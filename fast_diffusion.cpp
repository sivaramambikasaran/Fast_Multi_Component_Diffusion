#include<iostream>
#include<fstream>
#include<string>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include<Eigen/Dense>
#include"species.hpp"
#include"domain.hpp"

int main(int argc, char* argv[]) 
{
	srand(time(NULL));
	// The following array stores the name of the files from which data is to be read:
	std::string filename[] = {"data/table_224.txt","data/table_448.txt","data/table_896.txt","data/table_1792.txt","data/table_3585.txt","data/table_7171.txt"};
	// Used for timing:
	double CPS = CLOCKS_PER_SEC; // krithika
	clock_t start, end;

	// Iterating over each of the files:
	for (int k = 0; k < 6; ++k) 
	{
		domain myDomain;
		myDomain.read_Species(filename[k].c_str());
		std::cout << "----------------------------------------\n";
		std::cout << "Number of species is: " << myDomain.nSpecies << "\n";

		int nGrid			=	2;
		double left			=	-1.0;
		double right		=	1.0;
		double tolerance	=	1e-14;

		myDomain.generate_Grid(nGrid, left, right);

		start	=	clock(); //omp_get_wtime();  // krithika
		myDomain.compute_Species_Velocity(tolerance);
		end	=       clock(); //omp_get_wtime();  // krithika
		double fast_direct_time	=	(end-start)/CPS;


		start	=	clock(); //omp_get_wtime();  // krithika
		myDomain.compute_Exact_Species_Velocity();
		end	=       clock(); //omp_get_wtime();  // krithika
		double direct_time		=	(end-start)/CPS;

		myDomain.compute_Error();
		double fast_direct_error=	myDomain.grids[nGrid/2].error;

		double iterativeerror	=	0.0;
		start	=	clock(); //omp_get_wtime();  // krithika
		myDomain.compute_Species_Velocity_Iteratively(tolerance, iterativeerror);
		end	=       clock(); //omp_get_wtime();  // krithika
		double iterative_time	=	(end-start)/CPS;

		std::cout << "Time taken for: \n\n";
		std::cout << "\t Conv: \t" << direct_time << "\n";
		std::cout << "\t Fast: \t" << fast_direct_time << "\n";
		std::cout << "\t Iter: \t" << iterative_time << "\n";

		std::cout << "\nError in: \n\n";
		std::cout << "\t Fast: \t" << fast_direct_error << "\n";
		std::cout << "\t Iter: \t" << iterativeerror << "\n";

		std::cout << "\n----------------------------------------\n";
	}
}
