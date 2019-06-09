#include <iostream>
#include "FMDV.hpp"

int main(int argc, char* argv[])
{
    srand(time(NULL));

    int N = 100;
    FMDV F(N);

    F.setPressure(0.1);
    F.setPressureGradient(0.2);
    F.setDensity(0.3);
    F.setTemperature(0.4);
    F.setTemperatureGradient(0.5);

    Eigen::VectorXd X = Eigen::VectorXd(N).array() + 1;
    Eigen::VectorXd W = Eigen::VectorXd(N).array() + 1;

    F.setMolecularWeights(W);
    F.setDiameters(Eigen::VectorXd::Random(N).array() + 1);
    F.setThermalDiffusivities(Eigen::VectorXd::Random(N).array() + 1);
    F.setCollisionEnergies(Eigen::VectorXd::Zero(N));
    F.setMoleFraction(X);
    F.setMoleFractionGradient(Eigen::VectorXd::Random(N).array() + 1);
    F.setBodyForces(Eigen::VectorXd::Random(N).array() + 1);

    Eigen::MatrixXd V    = F.getInverseDiffusionCoefficientMatrix();
    Eigen::VectorXd v_sp = F.computeSpeciesVelocities(1e-12);

    // (diag(VX) - diag(X)V - S * Wt)z = (b - alpha * s)
    Eigen::VectorXd S       = F.getRandomVector();
	Eigen::MatrixXd M       = Eigen::MatrixXd((V * X).asDiagonal()) - X.asDiagonal() * V - S * W.transpose();
    Eigen::VectorXd z_exact = M.fullPivLu().solve(F.getRHS());

    double prefactor = 5 / (3 * 0.4);
    Eigen::VectorXd species_velocities_exact(N);
    for(int i = 0; i < N_species; i++)
    {
        if(X(i) == 0)
        {
            species_velocities_exact(i) = -prefactor * thermal_diffusivities(i) / mass_fraction(i);
        }

        else
        {
            species_velocities_exact(i) = z_exact(i) / mole_fraction(i) - prefactor * thermal_diffusivities(i) / mass_fraction(i);
        }
    }

    std::cout << "Error:" << (v_sp - v_sp_exact).cwiseAbs().maxCoeff() << std::endl;
}
