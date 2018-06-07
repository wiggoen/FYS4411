# Variational Monte Carlo of Fermionic systems

In collaboration with [Line](https://github.com/linegpe).

## Ways of running the program (you can change the parameters given by default in the 'parameters.json' file):
- Python script (you have to change to the library/include paths to fit your installation of [Armadillo](http://arma.sourceforge.net))
  - [Python 3](py3-run.py) [relative link]
- Manually by running it in [Qt Creator](https://www.qt.io) (you have to link [Armadillo](http://arma.sourceforge.net) to the .pro-file)
- Manually by using the command line with command line arguments (basically the same as the Python script provide)

Most of the parameter file is self-explanatory, but there are some hidden options.
"Cycle type" can take three different values
- "MonteCarlo" for running the code in normal mode calculating the energies.
- "SteepestDescent" for running the steepest descent method.
- "OneBodyDensity" for running the one body density calculations (this will write to a file called histogram.txt placed in the folder results inside the project folder).

"Cycle step to file" is an option to choose how much to write to file. The files written are the energies and the mean distances between the particles (this is only possible for two electrons at the moment).
- 0 is the option to not write to file.
- 1 is for writing every line to file.
- other values is for writing other kinds of lines to file (e.g. 1000 writes every 1000 energy and mean distance to file).
- When running "OneBodyDensity" it does not matter what the "Cycle step to file" is put on, one body density will write every line to its own file.

**NOTE!** If you don't have the folder 'results' in the project folder (Project_2), the files won't be written. Make sure to make the directory before you try to write to file.

If you run the script with <"UseMPI": true,> the program will include MPI for the calculations. But keep in mind that for now MPI only runs the same parameter file for all processors.

## Running automatic tests of the code
Change the first line in [parameters.json](parameters.json) [relative link] file to <"Run tests": true,>.
The program will run the tests when executed.
