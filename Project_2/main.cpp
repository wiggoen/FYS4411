#include <iostream>
#include <fstream>
#include <string>
#include "inc/json.hpp"

using json = nlohmann::json; // for convenience


int main()//int numberOfArguments, char *arguments[])
{
    /* The whole project path is needed to read json on mac */
    std::string projectFolder = "/home/line/github/FYS4411/Project_2/";
    // --

    /* READ PARAMETERS */
    std::ifstream inputFile(projectFolder+"parameters.json", std::ios::in);
    json parameter;
    inputFile >> parameter;
    double alpha = parameter["alpha"];
    std::cout << alpha << std::endl;




    std::cout << "Hello World!" << std::endl;
    return 0;
}
