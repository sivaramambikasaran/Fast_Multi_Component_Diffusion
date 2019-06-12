#include <iostream>
#include "../FMDV.hpp"

int main(int argc, char* argv[])
{
    srand(time(NULL));

    int N = 100;
    FMDV F(N);


double FMDV::getMatrixDiffusionCoefficient(int i, int j)
{
	// 1/D = ((d_i+d_j)/(2d_{max}))^2 * sqrt(2*m_i*m_j/(m_i+m_j)/m_{max})
    // TODO: Currently the constants have not been included. Need to do that
	double dia_sum              = diameters(i) + diameters(j);
	double elastic_matrix_entry = (temperature * sqrt(temperature) / pressure) * 
                                  dia_sum * dia_sum * 
                                  sqrt(2 * mass(i) * mass(j) / (mass(i) + mass(j)));

    // This segment computes collision integral based on collision energies:
    // TODO: Find out on what this is based
    if(collision_energies(i) != 0 && collision_energies(j) != 0)
    {
        static const double m[] = {6.8728271691,  9.4122316321,  7.7442359037,
                                   0.23424661229, 1.45337701568, 5.2269794238,
                                   9.7108519575,  0.46539437353, 0.00041908394781};
        
        double T  = temperature * sqrt(collision_energies(i) * collision_energies(j));
        double Nr = m[0] + T * (m[1] + T * (m[2] + T *  m[3]));
        double Dr = m[4] + T * (m[5] + T * (m[6] + T * (m[7] + T * m[8])));
    
    	return (Nr/Dr) * elastic_matrix_entry;
    }

    else
    {
        return elastic_matrix_entry;
    }
}

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
