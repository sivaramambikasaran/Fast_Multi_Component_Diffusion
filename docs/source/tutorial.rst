********
Tutorial
********

For this tutorial, we are going to be using the ``basic_example_usage/main.cpp`` file that is listed under ``examples/`` since it demonstrates all the features of this library. We go over the main features and functions that may be of interest to a user on this page.

This code can be used to solve for multicomponent diffusion velocities quickly from the Stefan-Maxwell's equation for a single grid point

.. math::
    \nabla X_p & = \overbrace{\displaystyle \sum_{k=1}^{N} \dfrac{X_pX_k}{\mathcal{D}_{pk}} \left( v_k - v_p \right)}^{\text{Difference in velocities}} + \underbrace{\displaystyle \sum_{k=1}^{N} \dfrac{X_pX_k}{\mathcal{D}_{pk}} \left( \dfrac{D^{(T)}_k}{Y_k} - \dfrac{D^{(T)}_p}{Y_p} \right) \dfrac{\nabla T}{\rho T}}_{\text{Temperature gradient (Soret effect)}} + \overbrace{\left(Y_p-X_p\right) \dfrac{\nabla P}{P}}^{\text{Pressure gradient}} + \underbrace{\dfrac{\rho}{P}\displaystyle \sum_{k=1}^N Y_pY_k \left( f_p - f_k \right)}_{\text{Difference in body force}}

This would mean that several input parameters such as mole-fraction, mass-fraction, pressure etc. are needed as an input to the solver.

Attributes Needed By Solver
---------------------------

The following attributes need to be a ``double`` since they describe data at the single grid point::

    // Density at the grid point:
    density
    // Temperature at the grid point:
    temperature
    // Temperature gradient at the grid point:
    temperature_gradient
    // Pressure at the grid point:
    pressure
    // Pressure gradient at the grid point:
    pressure_gradient

The following attributes need to be an ``Eigen::VectorXd`` of length ``N_species`` since it contains the values for each of the species used in the solver::
    
    // Molecular weights of the species:
    molecular_mass;
    // Diameters of the species:
    diameters;

    // Thermal diffusivities of the species:
    thermal_diffusivities;
    // Collision energies of the species:
    collision_energies;
    // Mole fraction of the species:
    // NOTE: A correct setup should ensure that this vector sums up to 1
    mole_fraction;
    // Mole fraction gradient of the species:
    // NOTE: A correct setup should ensure that this vector sums up to 0
    mole_fraction_gradient;
    // Mass fraction of the species:
    // NOTE: A correct setup should ensure that this vector sums up to 1
    mass_fraction;
    // Body forces of the species:
    body_forces;

Below, we go over how one can set the various parameters that are needed. The main solver object is a derived object of the main ``FMDV`` class. The various attributes needed for the solve along with the inverse diffusivity matrix is abstracted through this class. For the sake of the tutorial, we are calling this derived class ``fast_solver``. In this file, we have set all concerned parameters through the constructor. Alternatively, the parameters of interest can also be set / modified by changing these public attributes. The inverse diffusivity matrix is abstracted through the ``getInverseDiffusionCoefficient`` function which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix.

Declaring Derived Class of ``FMDV``
-----------------------------------

In ``main.cpp`` we have::

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

Here, we have set the attributes needed using the constructor itself. However, we can alternatively also just declare the object with the associated ``getInverseDiffusionCoefficient`` and set the attributes later before calling the solve function::

    class fast_solver : public FMDV
    {
    public:
        fast_solver(int N_species) : FMDV(N_species)
        {}

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

    fast_solver FS(N);

    // Density at the grid point:
    FS.density = double(rand()) / RAND_MAX;
    // Temperature at the grid point:
    FS.temperature = double(rand()) / RAND_MAX;
    // Temperature gradient at the grid point:
    FS.temperature_gradient = double(rand()) / RAND_MAX;
    // Pressure at the grid point:
    FS.pressure = double(rand()) / RAND_MAX;
    // Pressure gradient at the grid point:
    FS.pressure_gradient = double(rand()) / RAND_MAX;

    // Molecular weights:
    FS.molecular_mass = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
    FS.molecular_mass = FS.molecular_mass / FS.molecular_mass.maxCoeff();

    // Diameters:
    FS.diameters = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
    FS.diameters = FS.diameters / FS.diameters.maxCoeff();

    // Thermal diffusivities:
    FS.thermal_diffusivities = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
    // Collision energies:
    FS.collision_energies = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
    // Mole fraction:
    FS.mole_fraction  = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
    FS.mole_fraction /= FS.mole_fraction.sum();
    // Mole fraction gradient:
    mole_fraction_gradient = (Eigen::VectorXd::Random(N_species)).cwiseAbs();
    // Ensuring that molefraction gradient sums up to zero:
    FS.mole_fraction_gradient /= FS.mole_fraction_gradient.sum();
    FS.mole_fraction_gradient  = FS.mole_fraction_gradient.array() - 1 / double(N_species);

    // Mass fraction:
    FS.mass_fraction  = FS.molecular_mass.cwiseProduct(FS.mole_fraction);
    FS.mass_fraction /= FS.mass_fraction.sum();
    // Body forces:
    // The idea here is that body forces are proportional to mass of the species:
    FS.body_forces = double(rand()) / RAND_MAX * FS.mass_fraction;

Getting the Species Velocities
------------------------------

Once the solver object has been set with the required attributes, a single function call is all that is needed to get the diffusion velocities::

    Eigen::VectorXd v_sp = FS.computeSpeciesVelocities(1e-14);

**NOTE**: Currently we are performing the solve for a single grid point. However, when solving over several grid points, it is quite common for the considered physical problem to have a common matrix :math:`V` for all grid points. In such a case, further speedups can be acheived since the matrix doesn't have to be factorized at all grid points. Currently this interface doesn't support it in favour of simplicity, and would perform the factorization of the matrix for each grid point. If further speedups are required, this is worth looking into. We would be glad to aid in the development if needed :)
