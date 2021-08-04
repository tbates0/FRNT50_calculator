# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 17:37:00 2021

@author: Tim

This script takes an excel file with % neutralization values at different serum dilutions
and calculates the 50% focus reduction neutralization titer (FRNT50) of the combined replicates.
See the example data file for more information on input formatting

"""
#last tested in python 3.8.10 via anaconda
import numpy as np #last tested with numpy 1.20.2
from scipy.optimize import curve_fit #last tested with scipy 1.6.2
import matplotlib.pyplot as plt #last tested with matplotlib 3.3.4
import pandas as pd #last tested with pandas 1.2.5
import os, pathlib #os and pathlib are built in packages
    
'''

Useful functions

'''
#3 parameter logistic model for fitting and curve generation
def sigmoid(x, x0, k, z):
    '''
    Sigmoidal curve function
    1- at the beginning is such that the highest values are at low dilutions
    x = concentration input
    x0 = EC50 (log10)
    k = slope
    z = curve height
    if using concentration us the following formula:
        y = 1 / (1 + np.exp(-k*(x-x0)))
    '''
    y = z / (1 + np.exp(-k*(x-x0)))
    return y

#Generate the fit values for each of the samples, takes a single dataframe of samples as input
def fit(df):
    tempdf = []#temporary data holder
    df_fit = pd.DataFrame()#data frame to pass as output
    for i in df.axes[1]:#iterates over samples
        if i == "dilution" or i == "dilution-log":
            continue #attempting to fit these will throw errors
        if i not in things_to_plot:
            continue #skip things that we aren't plotting
        tempdf = df[["dilution-log",i]].dropna()
        x_fit = tempdf["dilution-log"]
        y_fit = tempdf[i]
        try:#this call can hpoave issues due to input data, its useful to see why it failed
            #bounds = [min][max] = [EC50(log10),slope,max height]
            #the syntax of the bounds are described above, they should be set manually for each type of data
            #note that EC50 is given as the log of EC50 while slope and height are normal
            popt, pcov = curve_fit(sigmoid,x_fit,y_fit,bounds=[[.01,-10,99],[5,-2,101]])
        except Exception as e:
            print (e)#tells you why it failed
        df_fit[i] = [popt,
                     pcov,#I dont currently use the covariance, but maybe you will
                     tempdf[["dilution-log"]].max(),
                     tempdf[["dilution-log"]].min()
                     ]
        print(i,popt)#keeps track of progress
    return(df_fit)#pass the dataframe of fit parameters out

#Generate the fit curve for plotting for each of the samples
def curve(df_fit):
    df_curve = pd.DataFrame()
    for i in df_fit.axes[1]:
        #x_max = df_fit.get(i)[2][0] #manually set maximum
        #x_min = df_fit.get(i)[3][0] #and minimum
        x_line = np.linspace(1, 4, num = 100) #makes 100 evenly spaced x values between 1 and 4
        y = df_fit.get(i)[0] #returns popt (curve fit parameters) for the sample wiht name i
        #The next line takes the x values and fit parameters to calculate y values for the fit curve
        #it uses the sigmoid function we defined above, the same one we used for curve fitting
        y_line = sigmoid(x_line, *y) #the * lets you pass the popt array directly to our sigmoid function
        df_curve[i] = [x_line, y_line]
    return(df_curve)#returns a dataframe with x,y values for the fit curves
    
'''

Data file and labeling

'''
#locate the csv file containing the data to be fit and plotted
#This is for windows systems (copy paste from explorer, including forward slashes)
data_folder = pathlib.Path(r"C:\Users\Tim\Desktop".replace('\\','/'))
#This is for Mac/Linux that use backslashes natively
#data_folder = pathlib.Path("/Users/Tim/Desktop")
infile = "Example data.xlsx"

#these are the variants being used. Make sure the sheet names match this list
#see the example spreadsheet for sheet and data layout
inputs = [
    'var1',
    'var2',
    'var3',
    ]

#what do you want to save the plot as?
#PDF saves as a unflattened (editable, vector) pdf file
outfile = 'combined-plots.pdf'
save = True #set to false to supress saving plots

#Plot labels
title = "Neutralization of SARS-CoV-2 by patient sera"
x_label = 'Log$_{10}$ of dilution factor'
y_label = "Percent neutralization"

'''

Data processing

'''

#initialize some data-frames
#Each variant gets its own dataframe and they are grouped into a dictionary
frames = {} #frame carries raw data
fits = {} #fits carries the curve fit parameters
curves = {} #curves carries x,y coordinates of the curves to plot
for i in inputs:
    frames[i] = pd.DataFrame()
    fits[i] = pd.DataFrame()
    curves[i] = pd.DataFrame()

#read the data
os.chdir(data_folder)
for i in inputs:
    frames[i] = dfWT = pd.read_excel(infile, sheet_name=i, na_values='VALUE!')

#This tells the graphing step which samples to plot
#by default it generates this list from the WT sheet, ignoring the dilution and catching a few other common issues
things_to_plot = []
for i in frames[inputs[0]].axes[1]:
    if i != 'dilution' and i[0:7] != 'Unnamed' and i[0] != "0" and i[0] != '_':
        things_to_plot.append(i)
things_to_plot.sort()

#This loop calculates the log of x values (concentration or dilution factors)
for i in inputs:
    frames[i]["dilution-log"]=np.log10(frames[i]["dilution"])

for i in inputs: #Runs the two functions consecutively
    fits[i]=fit(frames[i]) #fits logistic function to data
    curves[i]=curve(fits[i]) #makes x,y coordingates of function to plot

'''

Plot generation

'''
i=0#start sample iterable
tempdf=[] #used in plotting loop to temporarily store the averages and error bar sizes
#colors need to be defined like this to have dot and line colors be consistent within variants
color = ['#1f77b4', '#ff7f0e', '#2ca02c']
#adding different markers makes plots color blind friendly
symbol = ['o','s','v']
#initialize the plot
#Have to manually set the dimensions based on how many things you are plotting
fig,axes = plt.subplots(nrows=2,ncols=2,sharex='all',sharey='all',figsize=(14,16))

#This loop plots each sample, with all variants on the same plot
for row in axes:#iterates over rows in the plot
    for ax in row:#iterates through plots in the selected row
        try:
            thing = things_to_plot[i]#pick a sample based on the number of times through the loop
            #you could make this a for loop, but the level of indent is getting a bit silly
            for j,val in enumerate(curves):#plot the curve of the sample
                ax.plot(
                    curves[val][thing][0],#x values
                    curves[val][thing][1],#y values
                    c=color[j],
                    #label = "_hidden" #supresses the curve from making its own legend entry
                    )
            for j,val in enumerate(frames):#plot the averages with standard error bars
                tempdf = frames[val][["dilution-log",thing,"dilution"]].dropna() #pick the data for the sample
                #make the average and standard error
                means = tempdf.groupby("dilution").mean() #calculates means
                means["sem"] = tempdf.groupby("dilution").sem()[thing] #calculates standard error
                ax.errorbar(
                    means["dilution-log"],#x values
                    means[thing],#y values (average of all points at any given x value)
                    yerr = means["sem"], #error bars
                    fmt = symbol[j] #marker shapes
                    )
                '''#use this one if any points have no data (error bars plots can't run with null points)
                ax.scatter(
                    tempdf["dilution-log"],
                    tempdf[thing],
                    marker=symbol[j],
                    s = 20,
                    c=color[j]
                    )'''
            ax.set_title(thing)#puts a title of the sample name on each subplot
            #its easier to just manually define x and y ticks
            ticks = [1,2,3,4]
            ax.set_xticks(ticks)
            ax.set_xticklabels(ticks)
            ax.set_ylim(-5,105)
            ax.set_yticks([0,20,40,60,80,100])
            ax.set_yticklabels(['0','20','40','60','80','100'])
            i+=1 #moves to the next sample in the next subplot
        except Exception as e:
            print ("error")
            print(e)
            break

#put axis labels on to the plot as a whole rather than individual subplots
#axis and title names are defined at the top of the script
fig.text(0.5, 0.09, x_label, va='center', ha='center', fontsize=36)
fig.text(0.08, 0.5, y_label, va='center', ha='center', rotation='vertical', fontsize=36)
plt.figtext(.5,.93,title, fontsize=40, ha='center')
#These give a quick key for visualizing which color is which variant
fig.text(0.1,0.04,inputs[0],c=color[0],fontsize=36)
fig.text(0.3,0.04,inputs[1],c=color[1],fontsize=36)
fig.text(0.5,0.04,inputs[2],c=color[2],fontsize=36)

#This makes a list of EC50s in the console to be copied to excel
for i in inputs:#iterate over variants
    print(i)
    for j in things_to_plot:#iterates over samples
        print(10**(fits[i].get(j)[0][0]))#pops out and antilog transforms EC50 for sample j
print('Sample IDs')
for i in things_to_plot:#gives the order in which the above lists were output
    print(i)
    

#Saves as a layered pdf that can be edited in Illustrator
if save == True:
    fig.savefig(outfile, bbox_inches='tight')
