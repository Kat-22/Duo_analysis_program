# Duo_analysis_program
Repository for Master's Project on the use of the stabilisation method to interpret Duo outputs, initially used for the molecule SH.
Before using the programs, check that the code analyses your file names correctly or the values will be strange. It splits them based on underscores.

## Full_code_command_line.py:
Main program for the generation and analysis of Duo files. Generate Duo input files, or use Duo input, Duo output, or grepped/zipped extracted energies from output files.

See full documentation here (add link).

Use command 'help' in command line to see the required inputs. 

## Full_code_duo_analysis.py:
Alternative version of Full_code_command_line.py designed to be run from within python rather than the command line for easier analysis. Inputs and navigation are handled differently but output is identical.

See full documentation here (add link).

No help given for inputs, but full details given in documentation.

## Single_file_analysis.py:
Program to analyse a single state using molecule, state, j, v, and omega parameters post extraction using Full_code_command_line.py. 

See full documentation here (add link).

Use command 'help' in command line to see the required inputs.

