# -*- coding: utf-8 -*-
"""
FULL CODE FOR DATA EXTRACTION AND ANALYSIS - PYTHON VERSION
@author: Kathryn
"""
#Import libraries
import os, numpy as np, matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

#Get to directory
def change_directory(New_Path):
    bvariable = 0
    while bvariable == 0:
        Old_Path = os.getcwd()
        try:
            os.chdir(New_Path)
            bvariable = 1
        except OSError:
            print ("Please enter valid path")
            New_Path = input("Path: ")
    return Old_Path    
        
def correct_directory():            
    avariable = 0
    while avariable ==0:
        print("Current working directory is", os.getcwd())
        answer = input("Are you in the correct working directory?(Y/N): ")
        if  answer == "Y":
            Old_Path = os.getcwd()
            avariable = 1
            break
        elif answer == "N":
            New_Path = input("Path: ")
            Old_Path = change_directory(New_Path)
            avariable = 1
        else:
            avariable = 0
    return(Old_Path)

#Create Subfolder for extracted then sorted folders
def create_folder():
    CWD = os.getcwd()
    if not os.path.exists(f"{CWD}/Extracted_Files/Sorted"):
        os.makedirs(f"{CWD}/Extracted_Files/Sorted")

#Extract OUTput files and convert to list only (.out file)
def extract_file_out():
    create_folder()
    File_Directory= os.getcwd()
    for file in os.listdir(File_Directory):        
        splitfile = os.path.splitext(file)
        if splitfile[1]==".out":
              readingfile = open(file, "r")
              i = 0
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
                      i = i+1
                      print(i)
              readingfile.close()
              
#INPut to OUTput if INP files then to list only (.out file):
def run_duo(): 
    File_Directory= os.getcwd()
    for file in os.listdir(File_Directory):
        splitfile = os.path.splitext(file)
        if splitfile[1]==".inp":
            outputfile = splitfile[0]+".out"
            cmd = f"duo_dos2.exe < {file} > {outputfile}"
            os.system(cmd)
            
def extract_file_inp():
    #No option for generation here  
    run_duo()
    extract_file_out()

#Extract other file to out file
def extract_file_other(file_type):
    create_folder()
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
            #-----------------------------------------------------
            #From here to line check this makes sense with the file extraction!!!
            props = splitfile[0].split("_")               #this splits the filename by "_" to get molecule name and L value
            molecule = props[0]                           #this is the molecule name
            L_extra = props[1]                            #this is the L value
            L = L_extra.strip("L")
            for line in reading_file:
                line_parts = line.split(" ")
                #mol, state, pol, j, v, omega
                new_file_name_list = [molecule]
                line_parts[-1] = line_parts[-1].strip("||")
                line_parts[-1] = line_parts[-1].strip("\n")
                new_file_name_list.append(line_parts[-1]) #this is the state
                new_file_name_list.append(line_parts[-2]) #this is the polarity
                new_file_name_list.append(line_parts[0])  #this is the j value
                new_file_name_list.append(line_parts[4])  #this is the v value
                new_file_name_list.append(line_parts[-3]) #this is the omega value
                new_file_name = "&".join(new_file_name_list)
                line_to_write = []
                line_to_write.append(L)                   #this is the L value
                line_to_write.append(line_parts[2])       #this is the energy value
                line_to_write.append("\n")
                line_write = " ".join(line_to_write)
            #------------------------------------------
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
            print(i) #number of files you have dealt with

#GET DATA AND MODELS
def remove_anomalies(Energy_list, L_list): #this takes results that different from the main spread by a value of more than 100 and considers them anomalous.
    #Energy_ave = sum(Energy_list)/len(Energy_list)
    Energy_ave = np.median(Energy_list)
    Energy_list_new = []
    L_list_new = []
    for i in range(len(Energy_list)): #checking for strange deviations from the normal
        if abs(Energy_list[i]-Energy_ave)>100: #if the difference is too significant discard 
            Energy_ave = Energy_ave
        else:
            Energy_list_new.append(Energy_list[i])
            L_list_new.append(L_list[i])
    return(Energy_list_new, L_list_new)
    
def get_x_y(file): #Gets a list of Ls and energies for the state required
    L_list = []
    Energy_list = []
    with open(file, "r") as reading_file:
        for line in reading_file:
            splitline = line.split(" ")
            L_list.append(float(splitline[0]))
            Energy_list.append(float(splitline[1]))
            Title = os.path.splitext(file)[0].replace("&", " ")
    return(L_list, Energy_list, Title)

def lorentzian(energy, constant, energy_mean, energy_gamma): #Lorentzian function
    l_fit =  constant*(energy_gamma/((energy-energy_mean)**2 +(0.5*energy_gamma)**2))
    return(l_fit)

def get_loren_fit(Energy_list): #Generate histogram of energy and fit lorentzian function to it
    energy_freq, energy_bins = np.histogram(Energy_list, int(len(Energy_list)/6)+1)
    if sum(energy_freq) !=  0:
        energy_freq = energy_freq/sum(energy_freq)
    energy_values = []
    for i in range(len(energy_bins)-1):
        energy_val = (energy_bins[i]+energy_bins[i+1])/2
        energy_values.append(energy_val)
    energy_mean = sum(energy_values)/len(energy_values)
    if energy_freq.max() != 0:
        energy_gamma = 2/(np.pi*energy_freq.max())
    else:
        energy_gamma = 10000
    max_e = energy_freq.max()
    energy_freq = list(energy_freq)
    try:
        params2, matrix = curve_fit(lorentzian, energy_values, energy_freq, p0=[max_e,energy_mean,energy_gamma]) 
        matsd = np.sqrt(np.diag(matrix))  
    except RuntimeError: #If the fitting breaks it will not generate parameters
        params2 = [0,0,0]
        matsd = [0,0,0]
    except TypeError:
        params2 = [0,0,0]
        matsd =  [0,0,0]
    except ValueError:
        params2 = [0,0,0]
        matsd =  [0,0,0]
    #check for if the fit is bad (standard deviation of fits from point >1)    
    if matsd[2] > 10:
        params2 = [0,0,0]
        matsd = [0,0,0]
    return energy_freq, energy_values, params2, matsd

def get_loren_fit_small(Energy_list): # as above for very thin lines with little broadening -> i.e. v = 0
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
            weighting.append(0.01)
        else:
            weighting.append(1)
    try:
        params2, matrix = curve_fit(lorentzian, energy_values, energy_freq, p0=[max_e,energy_mean,energy_gamma], sigma=weighting, absolute_sigma=False) 
        matsd = np.sqrt(np.diag(matrix))
    except RuntimeError:
        params2 = [0,0,0]
        matsd = [0,0,0]
    except TypeError:
        params2 = [0,0,0]
        matsd = [0,0,0]
    except ValueError:
        params2 = [0,0,0]
        matsd =  [0,0,0]
    if matsd[2] > 10:
        params2 = [0,0,0]
        matsd = [0,0,0]
    return energy_freq, energy_values, params2, matsd

def get_line_broadening_l(params2,matsd): #get gamma (HWHM) and standard deviation
    HWHM = params2[2]
    HWHM = abs(HWHM)
    SD = matsd[2]
    return HWHM, SD

def plot_graph(L_list, Energy_list, Title, energy_freq, energy_values, params2, matsd, plots): #here you plot all graphs unlike other code because python handles it better than command line
            if plots == "Y":
                Fig1 = plt.figure() #L against energy
                plt.scatter(L_list, Energy_list)
                plt.xlabel("L")
                plt.ylabel("Energy")
                plt.title(Title)
                Fig2 =  plt.figure() #Line broadening
                plt.xlabel("Energy ($\mathregular{cm^{-1}}$)")
                plt.ylabel("Normalised Count")
                plt.axis([min(energy_values)-(max(energy_values)-min(energy_values))/50, max(energy_values)+(max(energy_values)-min(energy_values))/50, None, None])
                plt.scatter(energy_values, energy_freq)
            
            if [params2[0], params2[1], params2[2]]== [0,0,0]:
                resid = 0
            else:
                loren_profile = lorentzian(energy_values, *params2)
                if plots == "Y":
                    plt.plot(energy_values, loren_profile, "g")
            HWHM_l, SD = get_line_broadening_l(params2, matsd)
            Max_e = params2[1]
            Max_e_SD = matsd[1]
            os.chdir("..")
            os.chdir("..")
            with open("0Linelist.txt", "a") as f:
                if SD > 1.0:
                    f.write(f"{Title}: {HWHM_l:.6f} ± {SD:.6f} at {Max_e:.6f} ± {Max_e_SD:.6f} Query_Large_SD\n")
                elif HWHM_l >0 and SD/HWHM_l > 1:
                    f.write(f"{Title}: {HWHM_l:.6f} ± {SD:.6f} at {Max_e:.6f} ± {Max_e_SD:.6f} Query_Large_Rel_SD\n")
                else:
                    f.write(f"{Title}: {HWHM_l:.6f} ± {SD:.6f} at {Max_e:.6f} ± {Max_e_SD:.6f}\n")
            os.chdir(("./Extracted_Files/Sorted"))
            
#Run Program
def extraction():
    correct_directory()
    file_type = input("What is the input file extension? (with the dot please) ")
    if file_type == ".inp":
        extract_file_inp()
    elif file_type == ".out":
        extract_file_out()
    else:
        extract_file_other(file_type)
    get_keys_sorted()
    

def analysis():
    os.chdir("./Extracted_Files/Sorted")
    File_Directory= os.getcwd()
    state = input("What state would you like to analyse? ")
    plots = input("Would you like plots (Y/N)? ")
    while plots != "Y" and plots !="N":
        plots = input("Would you like plots (Y/N)? ")
    os.chdir("..")
    os.chdir("..")
    with open("0Linelist.txt", "w") as f:
        f.write("Mol State Pol J v Omega: Gamma ± StnDev at Peak_Energy ± StnDev Gamma_Warning\n")
    os.chdir("./Extracted_Files/Sorted")
    for file in os.listdir(File_Directory):
        splitfile = os.path.splitext(file)
        name_file = splitfile[0]
        name_file_split = name_file.split("&")
        try:
            if name_file_split[1] == state:
                L_list, Energy_list, Title = get_x_y(file)
                (Energy_list, L_list) = remove_anomalies(Energy_list, L_list)
                if int(name_file_split[4])== 0:
                    energy_freq, energy_values, params2, matsd = get_loren_fit_small(Energy_list)
                else:
                    energy_freq, energy_values, params2, matsd = get_loren_fit(Energy_list)
                plot_graph(L_list, Energy_list, Title, energy_freq, energy_values, params2, matsd, plots)
        except IndexError:
            state = state
            
    os.chdir("..")
    os.chdir("..")

def run_program():
    vers = 0
    while vers == 0:
        Full = input("Full Program including extraction? (Y/N) ") 
        if Full == "Y":
            vers = 1
            extraction()
            analysis()
        elif Full == "N":
            vers = 1
            correct_directory() #this should not go down to extracted files or the sorted, but the file with the original data in
            analysis()
        else:
            vers = 0

run_program()




