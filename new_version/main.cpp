#include <iostream>
#include "FMDV.hpp"

int main(int argc, char* argv[])
{
    int N = 1000;
    FMDV F(N);

    F.setPressure(0.1);
    F.setPressureGradient(0.2);
    F.setDensity(0.3);
    F.setTemperature(0.4);
    F.setTemperatureGradient(0.5);

    F.setMolecularWeights(Eigen::VectorXd::Random(N).array() + 1);
    F.setDiameters(Eigen::VectorXd::Random(N).array() + 1);
    F.setThermalDiffusivities(Eigen::VectorXd::Random(N).array() + 1);
    F.setCollisionEnergies(Eigen::VectorXd::Zero(N));
    F.setMoleFraction(Eigen::VectorXd::Random(N).array() + 1);
    F.setMoleFractionGradient(Eigen::VectorXd::Random(N).array() + 1);
    F.setBodyForces(Eigen::VectorXd::Random(N).array() + 1);

    Eigen::MatrixXd V    = F.getMatrix();
    Eigen::VectorXd v_sp = F.computeSpeciesVelocities(1e-12);

    // (diag(VX) - diag(X)V - S * Wt)z = (b - alpha * s)
    Eigen::VectorXd S = Eigen::VectorXd::Random(F.N_species);
    Eigen::VectorXd X = F.mole_fraction;
	Eigen::MatrixXd M = (V * X).asDiagonal() - X.asDiagonal() * V - S * F.mass.transpose();

    Eigen::VectorXd v_sp2 = M.fullPivLu().solve(F.rhs);
}
