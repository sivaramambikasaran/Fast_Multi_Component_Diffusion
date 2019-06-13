#include "FMDV.hpp"

class fast_solver : public FMDV
{
public:
    fast_solver(int N_species) : FMDV(N_species)
    {
        // Density at the grid point:
        density = double(rand()) / RAND_MAX;
        // Temperature at the grid point:
        temperature = double(rand()) / RAND_MAX;
        // Temperature gradient at the grid point:
        temperature_gradient = double(rand()) / RAND_MAX;
        // Pressure at the grid point:
        pressure = double(rand()) / RAND_MAX;
        // Pressure gradient at the grid point:
        pressure_gradient = double(rand()) / RAND_MAX;

        // Molecular weights:
        molecular_mass = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
        molecular_mass = molecular_mass / molecular_mass.maxCoeff();

        // Diameters:
        diameters = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
        diameters = diameters / diameters.maxCoeff();

        // Thermal diffusivities:
        thermal_diffusivities = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
        // Collision energies:
        collision_energies = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
        // Mole fraction:
        mole_fraction  = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
        mole_fraction /= mole_fraction.sum();
        // Mole fraction gradient:
        mole_fraction_gradient = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
        // Ensuring that molefraction gradient sums up to zero:
        mole_fraction_gradient /= mole_fraction_gradient.sum();
        mole_fraction_gradient  = mole_fraction_gradient.array() - 1 / double(N_species);

        // Mass fraction:
        mass_fraction  = molecular_mass.cwiseProduct(mole_fraction);
        mass_fraction /= mass_fraction.sum();
        // Body forces:
        // The idea here is that body forces are proportional to mass of the species:
        body_forces = double(rand()) / RAND_MAX * mass_fraction;
    }

    double getInverseDiffusionCoefficient(int i, int j)
    {
        // 1/D = K * P / T**1.5 * ((d_i+d_j))^2 * sqrt(m_i*m_j/(m_i+m_j))
        double dia_sum              = diameters(i) + diameters(j);
        double elastic_matrix_entry =   (pressure / (temperature * sqrt(temperature)))
                                      * dia_sum * dia_sum * sqrt(molecular_mass(i) * molecular_mass(j) 
                                      / (molecular_mass(i) + molecular_mass(j)));

        // This segment computes collision integral based on collision energies:
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
};

int main()
{
    int N = 400;
    fast_solver FS(N);
    
    Eigen::MatrixXd V    = FS.getInverseDiffusionCoefficientMatrix();
    Eigen::VectorXd v_sp = FS.computeSpeciesVelocities(1e-14);

    Eigen::MatrixXd L = FS.L;
    Eigen::MatrixXd R = FS.R;

    Eigen::VectorXd X = FS.mole_fraction;
    Eigen::VectorXd W = FS.molecular_mass;
    Eigen::VectorXd D = FS.thermal_diffusivities;
    Eigen::VectorXd Y = FS.mass_fraction;

    // (diag(VX) - diag(X)V - S * Wt)z = (b - alpha * s)
    Eigen::VectorXd S = FS.getRandomVector();
    Eigen::MatrixXd M = Eigen::MatrixXd((V * X).asDiagonal()) - X.asDiagonal() * V - S * W.transpose();

    Eigen::VectorXd z_exact = M.fullPivLu().solve(FS.getRHS());

    double prefactor = FS.temperature_gradient / (FS.density * FS.temperature);
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
