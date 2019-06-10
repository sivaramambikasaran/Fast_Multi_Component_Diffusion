#include <iostream>
#include "../FMDV.hpp"

int main(int argc, char* argv[])
{
    srand(time(NULL));

    int N = 100;
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
    F.setThermalDiffusivities(Eigen::VectorXd::Zero(N).array() + 0);
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
