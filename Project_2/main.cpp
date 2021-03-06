#define CATCH_CONFIG_RUNNER /* Configure Catch to use this main, and not its own */
#include "inc/variationalmontecarlo.h"
#include "inc/catch.hpp"
#include "inc/json.hpp"

#ifdef MPI_ON
#include <mpi.h>
#endif

using json = nlohmann::json; /* for convenience */

int RunCatchTests()
{
    return Catch::Session().run();
}


int main(int argumentCount, char *argumentVector[])
{
    /* Allocating input file for parameters */
    std::ifstream inputFile;

    /* Output argument number and argument to terminal */
    std::cout << "Arg. num " << "  Argument" << std::endl;
    for (int i = 0; i < argumentCount; i++)
    {
        std::cout << i << "          " << argumentVector[i] << std::endl;
    }


    /* Who am I? */
    std::string Iam = "Line";

    /* The whole project path is needed to read json file */
    std::string projectFolder;

    /* Check the number of command line arguments */
    if (argumentCount == 1) {
        if (Iam == "Trond")
        {
            std::cout << "Hi, Trond!" << std::endl;
            projectFolder = "/Users/trondwj/GitHub/FYS4411/Project_2/";         // Mac
            //projectFolder = "/home/twj/Documents/GitHub/FYS4411/Project_2/";  // Linux
            inputFile.open(projectFolder+"parameters.json", std::ios::in);
        } else if (Iam == "Line")
        {
            std::cout << "Hi, Line!" << std::endl;
            projectFolder = "/home/line/github/FYS4411/Project_2/";
            inputFile.open(projectFolder+"parameters.json", std::ios::in);
        } else
        {
            std::cerr << "You have to give the parameter file as command line argument." << std::endl;
            std::cerr << "Compile: g++ -std=c++11 -O3 -L '/armadillo/lib/' -I '/armadillo/include/' -I 'project folder' -o main.out main.cpp src/*.cpp -larmadillo" << std::endl;
            std::cerr << "Remember to include the paths of armadillo library and headers as stated above, and in addition the project folder." << std::endl;
            std::cerr << "Usage: ./main.out parameters.json" << std::endl;
            exit(1);
        }
    } else if (argumentCount == 2)
    {
        /* Input file from command line arguments */
        inputFile.open(argumentVector[1], std::ios::in);
    } else
    {
        /* Using MPI */
        /* Input file from command line arguments */
        //inputFile.open(argumentVector[5], std::ios::in);  /* Command line */

        if (Iam == "Trond")
        {
            projectFolder = "/Users/trondwj/GitHub/FYS4411/Project_2/";         // Mac
            //projectFolder = "/home/twj/Documents/GitHub/FYS4411/Project_2/";  // Linux
            inputFile.open(projectFolder+argumentVector[1], std::ios::in);
        } else if (Iam == "Line")
        {
            projectFolder = "/home/line/github/FYS4411/Project_2/";
            inputFile.open(projectFolder+argumentVector[5], std::ios::in);
        } else
        {
            inputFile.open(argumentVector[5], std::ios::in);
        }
    }

    /* READ PARAMETERS */
    json parameter;
    inputFile >> parameter;
    bool TEST = parameter["Run tests"];
    if (TEST) {
        std::cout << "Running tests..." << std::endl;
        return RunCatchTests();
    } else
    {
        int nDimensions = 2;
        bool UseAverageTiming = parameter["Use Average timing"];
        int nParticles = parameter["nParticles"];
        if (nParticles != 2 && nParticles != 6 && nParticles != 12 && nParticles != 20)
        {
            std::cout << std::endl;
            std::cerr << "The system is built to calculate the number of electrons in the quantity {2, 6, 12, 20}." << std::endl;
            std::cout << "Choose either of these particle numbers." << std::endl;
            std::cout << std::endl;
            exit(1);
        }
        int nCycles = parameter["nCycles"];
        double alpha = parameter["alpha"];
        double beta = parameter["beta"];
        double omega = parameter["omega"];
        double stepLength = parameter["stepLength"];
        bool UseImportanceSampling = parameter["Use Importance sampling"];
        double timeStep = parameter["timeStep"];
        bool UseJastrowFactor = parameter["Use Jastrow factor"];
        bool UseFermionInteraction = parameter["Use Fermion interaction"];
        bool UseAnalyticalExpressions = parameter["Use analytical expressions"];
        bool UseNumericalPotentialEnergy = parameter["Use numerical potential energy"];
        std::string cycleType = parameter["Cycle type"];
        int cycleStepToFile = parameter["Cycle step to file"];
        int terminalizationFactor = parameter["Terminalization factor"];


        /* Initialize VMC */
        VariationalMonteCarlo *VMC = new VariationalMonteCarlo();

        /* Allocation of information to be printed */
        arma::rowvec runVector;
        arma::mat runMatrix;

        /* For timing purposes */
        int trials;
        if (!UseAverageTiming)  /* No average timing */
        {
            trials = 1;
        } else
        {
            trials = 10;
        }

        /* MPI */
#ifdef MPI_ON
        {
            /* Initialize the MPI environment */
            MPI_Init(NULL, NULL);

            /* Get number of processes */
            int world_size;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);

            /* Get the rank of the process */
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

            /* Get the name of the processor */
            char processor_name[MPI_MAX_PROCESSOR_NAME];
            int name_len;
            MPI_Get_processor_name(processor_name, &name_len);

            /* Print off a hello world message */
            std::cout << "Hello world from processor " << processor_name << ", rank " << world_rank
                      << " out of " << world_size << " processors." << std::endl;
        }
#endif

        /* Run VMC */
        for (int i = 0; i < trials; i++)
        {
            runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                                    UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                                    UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
            runMatrix.insert_rows(i, runVector);
        }

        /* Used for running several trials */
        arma::rowvec columnSum = arma::sum(runMatrix, 0);

        /* Finding averages of trials */
        double runTime         = columnSum(0)/trials;
        double energy          = columnSum(1)/trials;
        double energySquared   = columnSum(2)/trials;
        double variance        = columnSum(3)/trials;
        double acceptanceRatio = columnSum(4)/trials;
        double kineticEnergy;
        double potentialEnergy;
        kineticEnergy   = columnSum(5)/trials;
        potentialEnergy = columnSum(6)/trials;

        /* Setup for printing to terminal */
        std::string Jastrow;
        UseJastrowFactor ? Jastrow = "On" : Jastrow = "Off";  /* (condition) ? (if_true) : (if_false) */

        std::string Importance;
        UseImportanceSampling ? Importance = "On" : Importance = "Off";

        std::string Interaction;
        UseFermionInteraction ? Interaction = "On" : Interaction = "Off";

        std::string Timing;
        UseAverageTiming ? Timing = "On" : Timing = "Off";

        std::string AnalyticExpressions;
        UseAnalyticalExpressions ? AnalyticExpressions = "On" : AnalyticExpressions = "Off";

        std::string NumericalPotentialEnergy;
        UseNumericalPotentialEnergy ? NumericalPotentialEnergy = "On" : NumericalPotentialEnergy = "Off";

        std::string WriteToFile;
        if (cycleStepToFile == 0) { WriteToFile = "Off"; }
        else                      { WriteToFile = "On"; }

        /* Write to terminal */
        std::cout << std::endl;
        std::cout << "Jastrow_factor "  << " Importance_sampling " << " Fermion_interaction " << " Average_timing "
                  << " Analytical_expressions " << " Numerical_potential_energy " << "      Cycle_type " << " Write_to_file "
                  << std::endl;

        std::cout << std::setw(14) << Jastrow
                  << std::setw(21) << Importance
                  << std::setw(21) << Interaction
                  << std::setw(16) << Timing
                  << std::setw(24) << AnalyticExpressions
                  << std::setw(28) << NumericalPotentialEnergy
                  << std::setw(17) << cycleType
                  << std::setw(15) << WriteToFile << std::endl;


        std::cout << "Particles "  << " Dimensions " << "    Cycles " << " Alpha " << " Beta " << " Omega "
                  << " Step_length " << " Time_step " << " Time_[sec] " << "      Energy "
                  << " Energy_squared " << " Variance "  << " Acceptance_ratio " << std::endl;

        std::cout << std::setw(9)  << std::setprecision(3) << nParticles
                  << std::setw(12) << std::setprecision(3) << nDimensions
                  << std::setw(11) << std::setprecision(8) << nCycles
                  << std::setw(7)  << std::setprecision(3) << alpha
                  << std::setw(6)  << std::setprecision(3) << beta
                  << std::setw(7)  << std::setprecision(3) << omega
                  << std::setw(13) << std::setprecision(3) << stepLength
                  << std::setw(11) << std::setprecision(6) << timeStep
                  << std::setw(12) << std::setprecision(6) << runTime
                  << std::setw(13) << std::setprecision(6) << energy
                  << std::setw(16) << std::setprecision(6) << energySquared
                  << std::setw(10) << std::setprecision(3) << variance
                  << std::setw(18) << std::setprecision(6) << acceptanceRatio << std::endl;


        std::cout << "Kinetic_energy "  << " Potential_energy " << std::endl;

        std::cout << std::setw(14) << std::setprecision(6) << kineticEnergy
                  << std::setw(18) << std::setprecision(6) << potentialEnergy << std::endl;

        /* Finalize the MPI environment */
#ifdef MPI_ON
        MPI_Finalize();
#endif

        return 0;
    }
}
