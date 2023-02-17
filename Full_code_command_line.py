# -*- coding: utf-8 -*-
"""
Command Line Version of Code 

@author: Kathryn
"""
#Use conda activate base in terminal pre testing
#Assume already in correct folder in directory

#Import libraries
import os, sys, csv, numpy as np, matplotlib.pyplot as plt
import warnings
import subprocess
import math
from scipy.optimize import curve_fit 
from scipy import stats
from pylab import *

#Create Subfolder for extracted then sorted folders
def create_folder():
    CWD = os.getcwd()
    if not os.path.exists(f"{CWD}/Extracted_Files/Sorted"):
        os.makedirs(f"{CWD}/Extracted_Files/Sorted")

#Extract OUTput files and convert to list only (.txt file)
def extract_file_out():
    File_Directory= os.getcwd()
    for file in os.listdir(File_Directory):        
        splitfile = os.path.splitext(file)
        if splitfile[1]==".out":
              readingfile = open(file, "r")
              for line in readingfile:
                  if line.find("||") == -1:
                      continue
                  else: 
                      line2 = line.split(" ")
                      new_line = []
                      for entry in line2:
                          if entry != "":
                              new_line.append(entry)
                      with open(f"Extracted_Files/{splitfile[0]}.out", "a") as f:
                          writeline = " ".join(new_line)
                          f.write(writeline)
              readingfile.close()
              
#INPut to OUTput if INP files then to list only (.txt file): (extend to generation?)
def run_duo(): 
    File_Directory = os.getcwd()
    for file in os.listdir(File_Directory):
        splitfile = os.path.splitext(file)
        if splitfile[1]==".inp":
            print(file)
            outputfile = f"{splitfile[0]}.out"
            cmd = (f"duo_dos2.exe < {file} > {outputfile}")
            print(cmd)
            os.system(cmd)
            
def extract_file_inp():
   #Submit input files to duo
    run_duo()
    #Take resulting output files and extract
    extract_file_out()


def sort_input_params(atoms, molecule, nstates, Jmax, N, L, template_file):
    try:
        molecule = str(molecule)
        atoms = str(atoms.replace(":"," "))
        nstates = int(nstates)
        Jmax_sp = Jmax.split(":")
        if len(Jmax_sp) ==1:
            Jmax_new = Jmax_sp[0]
            Jmax = float(Jmax_new)
        elif len(Jmax_sp) ==3:
            trial_Jmax = np.arange(float(Jmax_sp[0]),float(Jmax_sp[1]),float(Jmax_sp[2]))
            if trial_Jmax == []:
                sys.exit("Please ensure your J range goes from smallest value to largest value followed by step size")
            else:
                new_Jmax = f"{Jmax_sp[0]}, {Jmax_sp[1]} {Jmax_sp[2]}"
                Jmax = new_Jmax
        N = int(N)
        L_sp = L.split(":")
        if len(L_sp) ==1:
            Res = 0
        elif len(L_sp) ==3:
            Res = L_sp[2]
            new_L = np.arange(float(L_sp[0]),float(L_sp[1]),float(L_sp[2]))
            L = []
            for i in new_L:
                i = f"{i:.6f}"
                i = float(i)
                L.append(i)
            if L == []:
                sys.exit("Check that your L range is smaller_number:larger_number:step_size. ")
            if L[0]<0.85:
                sys.exit("In order to have a reasonable box size, the minimum value for L is 0.85, please ensure your range is not below this. ")  
        else:
            sys.exit("Your L should take the form 'number' or 'start_number:end_number:step_size'. See 'help' command for more details. ")
        template_file = template_file.strip('"') 
        try:
            with open(template_file, "r"):
                template_file = template_file
        except FileNotFoundError:
            sys.exit("Please ensure you have entered a valid file path to your template file. ")
        except OSError:
            sys.exit("Please ensure you have entered a valid file path to your template file. ")
    except IndexError:#TypeError:
        sys.exit("One of your input variables was of the incorrect type. Please refer to the 'help' command for more details. ")
    return(atoms, molecule, nstates, Jmax, N, L, Res, template_file)

def get_input_files(atoms, molecule, nstates, Jmax, N, L, Res, template_file):
    with open(template_file, "r") as f:
        for line in f:
            if line.find("atoms") != -1:
                new_line = f"atoms {atoms} \n"
            elif line.find("molecule") != -1:
                new_line = f"molecule {molecule} \n"
            elif line.find("nstates") != -1:
                new_line = f"nstates {nstates} \n"
            elif line.find("jrot") != -1:
                new_line = f"jrot {Jmax} \n"
            elif line.find("npoints") != -1:
                new_line = f"  npoints {N}  (odd) \n"
            elif line.find("range") != -1:
                new_line = f"  range  0.85,{L} \n"
            elif line.find("vmax") != -1:
                new_line = f"  vmax 40 10  {N} {N} {N} \n"
            elif line.find("INTENSITY") != -1:
                new_line = "INTENSITY OFF \n"
            else:
                new_line = line
            J = Jmax.strip(",")
            K = J.strip(":")
            J = K.strip(" ")
            with open(f"{molecule}_intTEST_L{L}_J_Res{Res}.inp", "a") as new_file:
                    new_file.write(new_line)
    
def generate_input_files(atoms, molecule, nstates, Jmax, N, L, template_file):
    atoms, molecule, nstates, Jmax, N, L, Res, template_file = sort_input_params(atoms, molecule, nstates, Jmax, N, L, template_file)
    i = 1
    for L_val in L:
        get_input_files(atoms, molecule, nstates, Jmax, N, L_val, Res, template_file)
        print(f"{i} of {len(L)}")
        i = i+1
    create_folder()
    

#Extract other file to txt file
def extract_file_other(file_type):
    File_Directory= os.getcwd()
    for file in os.listdir(File_Directory):        
        splitfile = os.path.splitext(file)
        if splitfile[1]==file_type:
              readingfile = open(file, "r")
              i = 0
              for line in readingfile:
                  line2 = line.split(" ")
                  new_line = []
                  for entry in line2:
                      if entry != "":
                          new_line.append(entry)
                  new_out_file = new_line[0].strip(":")
                  with open(f"Extracted_Files/{new_out_file}.out", "a") as f:
                      writeline = " ".join(new_line[1:])
                      f.write(writeline)
                  i = i+1
                  print(i)
              readingfile.close()

#Sort files into split by state
def get_keys_sorted():
    File_Dir= os.getcwd()
    File_Directory = f"{File_Dir}/Extracted_Files"
    i = 0
    for file in os.listdir(File_Directory): 
        splitfile = os.path.splitext(file)
        if splitfile[1]==".out":
            reading_file = open(f"{File_Directory}/{file}","r")
            props = splitfile[0].split("_")
            molecule = props[0]
            L_extra = props[2]
            L = L_extra.strip("L")
            for line in reading_file:
                line_parts = line.split(" ")
                #mol, state, pol, j, v, omega
                new_file_name_list = [molecule]
                line_parts[-1] = line_parts[-1].strip("||")
                line_parts[-1] = line_parts[-1].strip("\n")
                new_file_name_list.append(line_parts[-1])
                new_file_name_list.append(line_parts[-2])
                new_file_name_list.append(line_parts[0])
                new_file_name_list.append(line_parts[4])
                new_file_name_list.append(line_parts[-3])
                new_file_name = "&".join(new_file_name_list)
                line_to_write = []
                line_to_write.append(L)
                line_to_write.append(line_parts[2])
                line_to_write.append("\n")
                line_write = " ".join(line_to_write)
                check_counter = 0
                with open(f"./Extracted_Files/Sorted/{new_file_name}.out", "a") as e:
                    check_counter = check_counter
                with open(f"./Extracted_Files/Sorted/{new_file_name}.out", "r") as f:
                    for entry in f:
                        if entry == line_write:
                            check_counter = 1
                if check_counter == 0:
                    with open(f"./Extracted_Files/Sorted/{new_file_name}.out", "a") as g:
                        g.write(line_write)
            i = i+1
            print(i)

#GET DATA AND MODELS
def get_x_y(file):
    L_list = []
    Energy_list = []
    with open(file, "r") as reading_file:
        for line in reading_file:
            splitline = line.split(" ")
            L_list.append(float(splitline[0]))
            Energy_list.append(float(splitline[1]))
            Title = os.path.splitext(file)[0].replace("&", " ")
    return(L_list, Energy_list, Title)

def lorentzian(energy, constant, energy_mean, energy_gamma):
    l_fit =  constant*(energy_gamma/((energy-energy_mean)**2 + energy_gamma**2))
    return(l_fit)

def get_loren_fit(Energy_list):
    energy_freq, energy_bins = np.histogram(Energy_list, int(len(Energy_list)/6)+1)
    energy_freq = energy_freq/sum(energy_freq)
    energy_values = []
    for i in range(len(energy_bins)-1):
        energy_val = (energy_bins[i]+energy_bins[i+1])/2
        energy_values.append(energy_val)
    energy_mean = sum(energy_values)/len(energy_values)  
    energy_gamma = 2/(np.pi*energy_freq.max())
    max_e = energy_freq.max()
    energy_freq = list(energy_freq)
    try:
        params2, matrix = curve_fit(lorentzian, energy_values, energy_freq, p0=[max_e,energy_mean,energy_gamma]) 
        matsd = np.sqrt(np.diag(matrix))
    except RuntimeError:
        params2 = [0,0,0]
        matsd = [0,0,0]
    except TypeError:
        params2 = [0,0,0]
        matsd = [0,0,0]
    #check for if the fit is bad (standard deviation of fits from oint >100)    
    if matsd[2] > 1:
        params2 = [0,0,0]
        matsd = [0,0,0]
    return energy_freq, energy_values, params2, matsd

def get_loren_fit_small(Energy_list):
    energy_freq, energy_bins = np.histogram(Energy_list, int(len(Energy_list)/0.4)+1)
    energy_freq = energy_freq/sum(energy_freq)
    energy_values = []
    for i in range(len(energy_bins)-1):
        energy_val = (energy_bins[i]+energy_bins[i+1])/2
        energy_values.append(energy_val)
    energy_mean = sum(energy_values)/len(energy_values)  
    energy_gamma = 2/(np.pi*energy_freq.max())
    max_e = energy_freq.max()
    energy_freq = list(energy_freq)
    #add weighting here
    weighting = []
    for i in range(len(energy_values)):
        if energy_freq[i]/max_e >0.3:
            weighting.append(0.1)
        else:
            weighting.append(1)
    try:
        params2, matrix = curve_fit(lorentzian, energy_values, energy_freq, p0=[max_e,energy_mean,energy_gamma], sigma=weighting, absolute_sigma=True) 
        matsd = np.sqrt(np.diag(matrix))
    except RuntimeError:
        params2 = [0,0,0]
        matsd = [0,0,0]
    except TypeError:
        params2 = [0,0,0]
        matsd = [0,0,0]
    if matsd[2] > 0.1:
        params2 = [0,0,0]
        matsd = [0,0,0]
    return energy_freq, energy_values, params2, matsd

def get_line_broadening_l(params2,matsd):
    HWHM = params2[2]
    HWHM = abs(HWHM)
    SD = matsd[2]
    return HWHM, SD

def plot_graph(Title, params2, matsd):
            if [params2[0], params2[1], params2[2]]== [0,0,0]:
                params2 = params2
            else:
                HWHM_l, SD = get_line_broadening_l(params2,matsd)
                print(f"{Title}: {HWHM_l:.6f} ± {SD:.6f}")
           
            
#Run Program
def extraction(file_type):
    if file_type == ".out":
        extract_file_out()
    else:
        extract_file_other(file_type)
    get_keys_sorted()
    
def extraction_inp():
    extract_file_inp()
    get_keys_sorted()
    

def analysis(state):
    os.chdir("./Extracted_Files/Sorted")
    File_Directory= os.getcwd()
    #file_name ="SH&A2Sigma+&+&0.5&4&0.5"
    print("Molecule State Polarity J v Omega: gamma(HWHM) ± standard_deviation")
    for file in os.listdir(File_Directory):
        splitfile = os.path.splitext(file)
        name_file = splitfile[0]
        name_file_split = name_file.split("&")
        try:
            if name_file_split[1] == state:
                L_list, Energy_list, Title = get_x_y(file)
                if int(name_file_split[4])== 0:
                    energy_freq, energy_values, params2, matsd = get_loren_fit_small(Energy_list)
                else:
                    energy_freq, energy_values, params2, matsd = get_loren_fit(Energy_list)
                plot_graph(Title, params2, matsd)
        except IndexError:
            state = state
    #plt.show(block = False)
    os.chdir("..")
    os.chdir("..")

def run_program():
    print("You can ask for help using the command 'help'. ")
    inputs = sys.argv[1:]
    try:
        if inputs[0] == "help":
            print("The required inputs are in order:")
            print("1: Whether you want file extraction and analysis (Y) or just analysis (N). This is case sensitive. ")
            print("")
            print("2: If you answered 'Y' to the previous question then what file type needs to be extracted, including the period e.g. '.txt', '.out', or '.inp' ")
            print("")
            print("3: If you answered '.inp' to the previous question then you need to answer if you would like the program to generate the input files for you (Y) or not (N) ")
            print("")
            print("4: If you answered 'Y' to the previous question then there are several variables needed after that: ")
            print("a: atom names with a : in between of the desired molecule e.g. S:H ")
            print("b: molecule name e.g SH")
            print("c: the number of states for the molecule e.g 6 ")
            print("d: the maximum j value you would like to analyse e.g. 1.5, or range - see notes")
            print("e: the number of points you would like for the analysis to use e.g. 1001, this should be odd")
            print("f: the range or size of the box to use to approximate the continuum in angstroms e.g. 5 - see notes")
            print('g: the file path to the template input file for that molecule WITH QUOTATION MARKS e.g. "C:/path/to/file.inp". You need the full path unless the file is in the current working directory, this does not include if it is in a folder of the current working directory. ')
            print("")
            print("5: the state you want analysed exactly as in the input file e.g A2Sigma+") 
            print("")
            print("")
            print("NOTE - ranges")
            print("inputs 4d and 4f are the only possible multi variable inputs")
            print("for both of them they can either be a single number or a range")
            print("if they are a range they MUST take the form lowerbound:upperbound:stepsize e.g. 5:6:0.001 for L or 0.5:10.5:1 for J with ':' separating the values")
            print("for your information, the upperbound will not be included in the L range. ")
        elif  inputs[0] == "Y":
            if inputs[1] == ".inp":
                if inputs[2] == "N":
                    extraction_inp()
                    analysis(inputs[3])
                elif inputs[2] == "Y":
                    generate_input_files(inputs[3],inputs[4],inputs[5],inputs[6],inputs[7],inputs[8],inputs[9])
                    extraction_inp() 
                    analysis(inputs[10])
                else:
                    sys.exit("If you are trying in use .inp files please ensure the third input is 'Y' or 'N'. Refer to the 'help' command for more details ") 
            elif inputs[1][0] != ".":
                sys.exit("Please ensure that the second input is a valid file extension that includes the period")
            else:
                extraction(inputs[1])
                analysis(inputs[2])
        elif inputs[0] == "N":
            analysis(inputs[1])
        else:
            sys.exit("Please ensure your first input is either 'Y' or 'N', or refer to help document")
    except IndexError:
        sys.exit("Please ensure you have entered the required number of inputs (2,3,4, or 11). Use 'help' command for more details ")

    
run_program()   





