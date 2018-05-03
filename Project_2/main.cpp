#include <iostream>
#include <fstream>
#include <string>
#include "inc/json.hpp"

using json = nlohmann::json; // for convenience


int main()//int numberOfArguments, char *arguments[])
{
    /* The whole project path is needed to read json on mac */

    /* Who am I? */
    std::string Iam = "Line";
    std::string projectFolder;
    if (Iam == "Trond") {
        std::cout << "Hi, Trond!" << std::endl;
        projectFolder = "/Users/trondwj/GitHub/FYS4411/Project_2/";
    } else if (Iam == "Line") {
        std::cout << "Hi, Line!" << std::endl;
        projectFolder = "/home/line/github/FYS4411/Project_2/";
    } else {
        std::cerr << "I don't know who you are.." << std::endl;
        exit(1);
    }


    /* READ PARAMETERS */
    std::ifstream inputFile(projectFolder+"parameters.json", std::ios::in);
    json parameter;
    inputFile >> parameter;
    int nParticles = parameter["nParticles"];
    int nDimensions = parameter["nDimensions"];
    int nCycles = parameter["nCycles"];
    double alpha = parameter["alpha"];
    double beta = parameter["beta"];


    return 0;
}
