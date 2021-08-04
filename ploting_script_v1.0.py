# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 17:37:00 2021

@author: Tim

This script takes an excel file with % neutralization values at different serum dilutions
and calculates the FRNT50 of the combined replicates.
See 

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
infile = "Source data.xlsx"

#these are the variants being used. Make sure the sheet names match this list
#see the example spreadsheet for sheet and data layout
inputs = [
    'WT',
    'UK',
    'SA',
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
for i in frames['WT'].axes[1]:
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
fig,axes = plt.subplots(nrows=9,ncols=6,sharex='all',sharey='all',figsize=(28,32))

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
fig.text(0.1,0.04,'WT',c=color[0],fontsize=36)
fig.text(0.3,0.04,'UK',c=color[1],fontsize=36)
fig.text(0.5,0.04,'SA',c=color[2],fontsize=36)

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


"""depreciated functions (use at your own risk)

def normalized_OD(raw, plate_bg=0.21, scale_cofactor=1, curve_max=8.2):
    '''
    Normalization function 
    Subtracts minimum value in array, assuming that represents plate background. 
        If signal hasn't decreased to plate background, call function with assigned value for plate_bg
    Then, normalizes top value to 1, assuming max value is plateau of EC50 curve. 
        If not, I ADVISE AGAINST PLOTTING NORMALIZED. Instead, input raw array into plotELISA.
        Sometimes, if you have enough other points to know where top curve is, use 0<scale_cofactor<1 to set scale, 
        estimating where points should hit along curve.
    '''
    if plate_bg == []:
        plate_bg = min(raw)
    raw_bg = [x-plate_bg for x in raw]
    #factor = scale_cofactor/max(raw_bg)
    factor = 1/(curve_max-plate_bg)
    norm = [x*factor for x in raw_bg]
    return norm

def molarity_dilutions(startnM, dilutions, factor):
    '''
    Generates an array of log concentrations (nM) matching absorbance array dilutions, used as x-values for scatterplot.
    Inputs:
        startnM - starting VHH conc (nM)
        dilutions - number of dilutions (start is 0, so if 12 wells, input 12)
        factor - dilution factor 
    '''
    concentrations = []
    for i in range(dilutions):
        concentrations.append(startnM/(factor**i))
    log_conc = [np.log10(x) for x in concentrations]
    return log_conc

def calculate_EC50(data, VHH_name='', xlim_low=-1.5, xlim_high=1):
    '''
    Uses input pairs of molarity and absorbance value LISTS in order to incorporate disparate experimental reps - See example if confused
        (not in form of [100nm, OD=5], but [[list of mol for 01-18 exp, list of abs values for 01-18 exp],[01-23 mol, 01-23 abs]])
    Optional inputs:
        VHH_name - str of ID for VHH to use in file-name, title, and to plot with other VHHs
        xlim_low and xlim_high - number values
            bounds are in nM, because fit function breaks if there is no start conc is below ~10 (num value), 
            so it is easier to represent in nM
    '''
    scatterdata = []
    for rep in data:
        pairs = list(np.transpose(rep))
        for pair in pairs:    
            scatterdata.append(list(pair))  #Creates one big list of pairs of mol and abs to be scatterplotted
    

    dilutions = list(np.transpose(scatterdata)[0]) #separates x data into single list
    absorbances = list(np.transpose(scatterdata)[1]) #separates y data into single list

    popt, pcov = curve_fit(sigmoid, dilutions, absorbances)

    x = np.linspace(xlim_low,xlim_high, num = 100)
    y = sigmoid(x, *popt)
    
    print('EC50 of ' + VHH_name + ' = ' + str(10**(popt[0]))[:5] + " ug/mL")

    # create plot
    fig, ax = plt.subplots()

    plt.scatter(dilutions,absorbances, color= 'xkcd:cornflower blue')
        
    plt.plot(x,y, c = 'xkcd:dirty orange')
    #ax.set_xscale('log')
                 
    ax.plot([xlim_low, xlim_high], [0.5,0.5], "k--")    #Adds line at 0.5 to mark EC50
    ax.set_xticklabels(['.032','.100','.316','1.00','3.16','10.0'])
    plt.xlim((xlim_low, xlim_high))                #Set x-limits of plot (lowbound, highbound)  NOTE: in M, not nM
    plt.ylim((-0.1,1.2))                        #Set y-limits of plot (lowbound, highbound) 
    plt.xlabel('Antibody Concentration ug/mL')           #Xlabel (string)
    plt.ylabel('Normalized OD450 Absorbance')   #Ylabel (string)
    plt.title(VHH_name+' Best Fit')                  #Plot title, input
    plt.text(-1.4,1, "EC50 = " + str(10**(popt[0]))[:5] + " ug/mL", )
        
    plt.tight_layout()
    plt.show()
#    fig.savefig(VHH_name+' Best Fit.pdf')
    fig.savefig(VHH_name+' Best Fit.png')
    
#This function plots the curves and points
x_axis_size = []
ec50list = []
color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#000000']
for idx, val in enumerate(things_to_plot):
    #Plot the fit curve
    if curve == True:
        plt.plot(
            df_curve[val][0], 
            df_curve[val][1], 
            label = "_hidden", #does not show up in legend
            c = color[idx],
            )
    if points == True:
        tempdf = df[["dilution-log",val,"dilution"]].dropna()
        plt.scatter(
            tempdf[["dilution-log"]], 
            tempdf[[val]],
            label = val,
            c = color[idx]
            )
        x_axis_size = np.append(x_axis_size, tempdf.get('diluiton').unique())
        x_axis_size = np.unique(x_axis_size)
    if errorbar == True:
        tempdf = df[["dilution-log",val,"dilution"]].dropna()
        #make the average and standard error
        means = tempdf.groupby("dilution").mean()
        means["sem"] = tempdf.groupby("dilution").sem()[val]
        plt.errorbar(
            means["dilution-log"],
            means[val],
            yerr = means["sem"],
            label = val,
            fmt = "o",
            c = color[idx],
            capsize = 3
            )
        x_axis_size = means.index     
    if EC_50 == True:
        ec50str = 'EC$_{50}$ of ' + val + ' = ' + str(10**(df_fit.get(val)[0][0]))[:5]
        print(ec50str)
        ec50list.append(ec50str)

#setting the x axis to include the highest and lowest values in the plotted data
ticknums = np.sort(x_axis_size)
ticklabs=[]
for i in ticknums:
    a ='%.5f'%i
    b = a[:4].rstrip('0')
    if a[0] == '0':
        if a[2] == '0':
            b = a[:5].rstrip('0')
        if a[3] == '0':
            b = a[:5].rstrip('0')
        if a[4] == '0':
            b = a[:6].rstrip('0')
    c = b.rstrip('.')
    ticklabs.append(c)


#Labeling parameters
#plt.gca().invert_xaxis()
#plots on a log axis with un-logged labels
fig.xticks(ticks = [np.log10(i) for i in ticknums], labels = ticklabs)
fig.xlabel(x_label)
fig.ylabel(y_label)
ec50s='\n'.join(ec50list)
#plt.ylim(-0.05,1.05)
"""
