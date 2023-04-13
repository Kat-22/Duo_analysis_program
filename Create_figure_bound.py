# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 10:16:41 2023
Figure creation
@author: Kathryn
"""

#Import libraries
import os, matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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

#Generate linelist for this figure
def create_linelist():
    change_directory("./Extracted_Files")
    lista = []
    for file in os.listdir(os.getcwd()):
        splitfile = os.path.splitext(file)
        if splitfile[1] == ".out":
            section = []
            with open(f"./{file}", "r") as f:
                for line in f:
                    section.append(line)
                    if len(section) > 5:
                        section.pop(0)
                    if len(section) < 5:
                        section = section
                    else:
#This should be the specific search string for the bound state you 
#want to see the interactions with, followed by the J value you are interested
#in followed by a space
#------------------------------------------------------------------------------
                        if section[2].find("2 4 0 0.5 0.5 0.5 + ||A2Sigma+") != -1: 
                            if section[2][0:3] =="1.5":
#------------------------------------------------------------------------------
                                upperupperbound = section[0].strip("\n")
                                upperbound = section[1].strip("\n")
                                actualbound = section[2].strip("\n")
                                lowerbound = section[3].strip("\n")
                                lowerlowerbound = section[4]
                                lista.append([upperupperbound, upperbound, actualbound, lowerbound, lowerlowerbound])
    change_directory("..")
    with open("./Linelist_figure.txt", "w") as g:
        for entry in lista:
             new_line = "\t".join(entry)
             g.write(new_line)
#GET DATA AND MODELS
def get_bound():
#Choose the summary file to analyse (should be Linelist_figure.txt in the data folder)
#------------------------------------------------------------------------------
    with open("./Linelist_figure.txt", "r") as h:
#------------------------------------------------------------------------------
        unbound = []
        unboundish = []
        for line in h:
            linesplit = line.split("\t")
            unboundish.append(linesplit[0])
            unboundish.append(linesplit[1])
            unboundish.append(linesplit[3])
            unboundish.append(linesplit[4].strip("\n"))
        for line in unboundish:
            linesp = line.split(" ")
            new_line = ["SH"]
            state = linesp[-1].strip("\n")
            new_line.append(state.strip("||"))
            new_line.append(linesp[-2]) #polarity
            new_line.append(linesp[0])  #J
            new_line.append(linesp[4])  #v
            new_line.append(linesp[-3]) #omega
            nline = "&".join(new_line)
            counter = 0
            for entry in unbound:
                if entry == nline:
                    counter = 1
                    break
            if counter == 0:
                unbound.append(nline)
    return(unbound)

def get_x_y(file):
    L_list = []
    Energy_list = []
    with open(file, "r") as reading_file:
        for line in reading_file:
            splitline = line.split(" ")
            L_list.append(float(splitline[0]))
            Energy_list.append(float(splitline[1]))
            Title_sp = os.path.splitext(file)[0]
            Title = Title_sp.split("&")[1]
    return(L_list, Energy_list, Title)

def plot_graph(L_list, Energy_list, Title):
#If running this is python, these x limits and y limits can be chosen to 
#display the interactions with the bound state in more detail
#------------------------------------------------------------------------------
            #plt.xlim(10,10.25)
            #plt.ylim(36000,39000)
#------------------------------------------------------------------------------
#These Titles should be chosen for the bound and repulsive state names in the 
#extracted file
#------------------------------------------------------------------------------
            if Title == "A2Sigma+":
                plt.scatter(L_list, Energy_list, marker = ".", color = "r")
            elif Title == "b4Pi":
                plt.scatter(L_list, Energy_list, marker = ".", color = "b")
            elif Title == "1^2_Sigma-":
                plt.scatter(L_list, Energy_list, marker = ".", color = "g")
            elif Title == "a4Sigma-":
                plt.scatter(L_list, Energy_list, marker = ".", color = "k")
#------------------------------------------------------------------------------
            plt.xlabel("L")
            plt.ylabel("Energy")
            plt.title("Interactions between bound and replusive states")
            
            
#Run Program  
def analysis():
    correct_directory() #This should be the data folder
    line_list_creation = "0"
    while line_list_creation != "Y" and line_list_creation != "N":
        line_list_creation = input("Do you need to create the Linelist_figure.txt file? ")
    if line_list_creation == "Y":
        create_linelist()
    file_name_list = get_bound()
    for state_file in file_name_list:
#This should be the bound state to remove duplicates as repulsive states are
#the interactions to be displayed
#------------------------------------------------------------------------------
        if state_file.split("&")[1] == "A2Sigma+":
            file_name_list.pop(file_name_list.index(state_file))
#This should be the bound state you are interested in - this will ensure this
#state is above all the other repulsive states
#------------------------------------------------------------------------------
    file_name_list.append("SH&A2Sigma+&+&1.5&4&0.5")
#------------------------------------------------------------------------------
    change_directory("./Extracted_Files/Sorted")
    for file_name in file_name_list:
        for file in os.listdir(os.getcwd()):
            splitfile = os.path.splitext(file)
            name_file = splitfile[0]
            if name_file == file_name:
                L_list, Energy_list, Title = get_x_y(f"./{file}")
                plot_graph(L_list, Energy_list, Title)
#This creates the legend for the graph so ensure it matches the colour (color) and
#title (label) above
#------------------------------------------------------------------------------
    states_legend = [Line2D([0], [0], color="r", lw=2, label = "A2Sigma+"),
                     Line2D([0], [0], color="b", lw=2, label = "b4Pi" ),
                     Line2D([0], [0], color="g", lw=2, label = "1^2Sigma-"),
                     Line2D([0], [0], color="k", lw=2, label = "a4sigma-")]
#------------------------------------------------------------------------------
    plt.legend(handles = states_legend)
    plt.show(block = True)
    os.chdir("..")
    os.chdir("..")

analysis()






