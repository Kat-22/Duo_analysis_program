# -*- coding: utf-8 -*-
"""
Single state analysis
@author: Kathryn
"""
import os, sys, csv, numpy as np, matplotlib.pyplot as plt
import warnings
import math
from scipy.optimize import curve_fit 
from scipy import stats
from pylab import *

#print help if required
def helpme():
    print("To use this function the following inputs need to be used: ")
    print('1: file name IN QUOTATION MARKS EXCLUDING extension e.g "SH&A2Sigma+&+&0.5&3&0.5"')
    print("OR")
    print("1: molecule name e.g SH ")
    print("2: state name e.g. A2Sigma+ ")
    print("3: J value e.g. 0.5")
    print("4: v value e.g. 3")
    print("5: Omega value e.g. 0.5")

#Check inputs are correct
def check_inputs_file(inputs):
    file = "".join(inputs)
    new_file = file.strip('"')
    file = new_file
    try:
        os.chdir("./Extracted_Files")
        os.chdir("./Sorted")
        with open(f"{file}.out", "r"):
            file = file  
        os.chdir("..")
        os.chdir("..")
    except FileNotFoundError:
        sys.exit("Please ensure you have entered a valid file name and that the directory you are in has the Extracted_Files folder in it. ")
    except OSError:
        sys.exit("Please ensure you have entered a valid file name and that the directory you are in has the Extracted_Files folder in it. ")
    return file

def check_inputs_no_file(inputs):
    try:
        molecule = str(inputs[0])
        state = str(inputs[1])
        J_val = float(inputs[2])
        v_val = int(inputs[3])
        o_val = float(inputs[4])
        return(inputs)
    except TypeError:
         sys.exit("One of your input variables was of the incorrect type. Please refer to the 'help' command for more details. ")
    except IndexError:
        sys.exit("Please ensure you have entered the required number of inputs (1 or 5). Use 'help' command for more details ")

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
    if matsd[2] > 0.1:
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

def plot_graph(L_list, Energy_list, Title, energy_freq, energy_values, params2, matsd):
            Fig1 = plt.figure()
            plt.scatter(L_list, Energy_list)
            plt.xlabel("L (Å)")
            plt.ylabel("Energy ($\mathregular{cm^{-1}}$)")
            plt.title(Title)
            
            Fig2 =  plt.figure()
            plt.axis([min(energy_values)-(max(energy_values)-min(energy_values))/50, max(energy_values)+(max(energy_values)-min(energy_values))/50, None, None])
            plt.scatter(energy_values, energy_freq)
            plt.title(Title)
            plt.xlabel("Energy ($\mathregular{cm^{-1}}$)")
            plt.ylabel("Normalised Count")
            if [params2[0], params2[1], params2[2]]== [0,0,0]:
                params2 = params2
            else:
                loren_profile = lorentzian(energy_values, *params2)
                plt.plot(energy_values, loren_profile, "g")
                FWHM_l, SD = get_line_broadening_l(params2,matsd)
                print(f"{Title}: {FWHM_l:.6f} ± {SD:.6f}")
                
def single_file_analysis(inputs):
    os.chdir("./Extracted_Files/Sorted")
    File_Directory= os.getcwd()
    if len(inputs) == 1:
        file_name = inputs[0]
        try: 
            for file in os.listdir(File_Directory):
                splitfile = os.path.splitext(file)
                name_file = splitfile[0]
                name_file_split = name_file.split("&")
                if name_file == file_name:
                    L_list, Energy_list, Title = get_x_y(file)
                    if int(name_file_split[4]) == 0:
                        energy_freq, energy_values, params2, matsd = get_loren_fit_small(Energy_list)
                    else:
                        energy_freq, energy_values, params2, matsd = get_loren_fit(Energy_list)
                    plot_graph(L_list, Energy_list, Title, energy_freq, energy_values, params2, matsd)
        except IndexError:
            sys.exit("Index Error")
    else:
        for file in os.listdir(File_Directory):
            splitfile = os.path.splitext(file)
            name_file = splitfile[0]
            name_file_split = name_file.split("&")
            try:
                if name_file_split[0]== inputs[0]:
                    if name_file_split[1]== inputs[1]:
                        if name_file_split[3]== inputs[2]:
                            if name_file_split[4]== inputs[3]:
                                if name_file_split[5]== inputs[4]:
                                    L_list, Energy_list, Title = get_x_y(file)
                                    if int(name_file_split[4])== 0:
                                        energy_freq, energy_values, params2, matsd = get_loren_fit_small(Energy_list)
                                    else:
                                        energy_freq, energy_values, params2, matsd = get_loren_fit(Energy_list)
                                    plot_graph(L_list, Energy_list, Title, energy_freq, energy_values, params2, matsd)
            except IndexError:
                sys.exit("Index Error")
    plt.show(block = True)
    os.chdir("..")
    os.chdir("..")
                
def run_program():
    inputs = sys.argv[1:]
    print("Run this program from the same directory point as you ran the file extraction (i.e. do not go into the Extracted_Files or Sorted folders for this. Use the 'help' command for input requirements")
    try:
        if inputs[0] == "help":
            helpme()
        else:
            if len(inputs) == 1:
                new_inputs = check_inputs_file(inputs)
                new_i = []
                new_i.append(new_inputs)
                new_inputs = new_i
            else:
                new_inputs = check_inputs_no_file(inputs)
                
            single_file_analysis(new_inputs)
    except IndexError:
        sys.exit("Please ensure you have entered the required number of inputs (1 or 5). Use 'help' command for more details ")

run_program()
    
