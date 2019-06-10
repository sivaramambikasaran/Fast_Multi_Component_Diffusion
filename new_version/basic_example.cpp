#include <iostream>
#include "FMDV.hpp"

int main(int argc, char* argv[])
{
    srand(time(NULL));

    int N = 50;
    FMDV F(N);

    F.setPressure(0.1);
    F.setPressureGradient(0.2);
    F.setDensity(0.3);
    F.setTemperature(0.4);
    F.setTemperatureGradient(0.5);

    // Vectors with random values in range [0, 1]
    F.setMolecularWeights(0.5 * (Eigen::VectorXd::Random(N).array() + 1));
    F.setDiameters(0.5 * (Eigen::VectorXd::Random(N).array() + 1));
    F.setThermalDiffusivities(0.5 * (Eigen::VectorXd::Random(N).array() + 1));
    F.setCollisionEnergies(Eigen::VectorXd::Zero(N));
    F.setMoleFraction(0.5 * (Eigen::VectorXd::Random(N).array() + 1));
    F.setMoleFractionGradient(0.5 * (Eigen::VectorXd::Random(N).array() + 1));
    F.setBodyForces(0.5 * (Eigen::VectorXd::Random(N).array() + 1));

    // Fast Method:
    Eigen::VectorXd v_sp = F.computeSpeciesVelocities(1e-12);

    // (diag(VX) - diag(X)V - S * Wt)z = (b - alpha * s)
    Eigen::VectorXd X       = F.getMoleFraction();
    Eigen::VectorXd W       = F.getMolecularWeights();
    Eigen::MatrixXd V       = F.getInverseDiffusionCoefficientMatrix();
    Eigen::VectorXd S       = F.getRandomVector();
	Eigen::MatrixXd M       = Eigen::MatrixXd((V * X).asDiagonal()) - X.asDiagonal() * V - S * W.transpose();
    Eigen::VectorXd z_exact = M.fullPivLu().solve(F.getRHS());

    // double prefactor = 5 / (3 * 0.4);
    // Eigen::VectorXd species_velocities_exact(N);
    // for(int i = 0; i < N_species; i++)
    // {
    //     if(X(i) == 0)
    //     {
    //         species_velocities_exact(i) = -prefactor * thermal_diffusivities(i) / mass_fraction(i);
    //     }

    //     else
    //     {
    //         species_velocities_exact(i) = z_exact(i) / mole_fraction(i) - prefactor * thermal_diffusivities(i) / mass_fraction(i);
    //     }
    // }

    std::cout << v_sp.sum() << std::endl;
    std::cout << z_exact.sum() << std::endl;

    std::cout << "Error:" << (v_sp - z_exact).cwiseAbs().maxCoeff() << std::endl;
}
