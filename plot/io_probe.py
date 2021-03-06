# Written 30/09/17 by dh4gan
# Module for handling data output from lotka_volterra_probes (C++)

import numpy as np
import glob
import re
import matplotlib.pyplot as plt
from scipy.signal import blackman, periodogram

itime = 0
iprey = 1
ipred = 2
ipreydot = 3
ipreddot = 4
ipreyout = 5
ipredout = 6

preycolor = '#669900'
predcolor = '#ff3300'


xlabel = 'x (kpc)'
ylabel = 'y (kpc)'

timelabel = 'Time (Myr)'
numlabel = 'Number (Thousands)'

preylabel = 'Prey (Thousands)'
predlabel = 'Predators (Thousands)'

preyleg = 'Prey'
predleg = 'Predators'

axesize=20

def sort_nicely(l):
    """
    Sort the given list in a 'Natural' order
    (Ned Batchelder's Compact Human Python Sort)
    """
    
    convert = lambda text:int(text) if text.isdigit() else text
    
    alphanum_key = lambda key:[ convert(c) for c in re.split('([0-9]+)',key)]
    
    l.sort(key=alphanum_key)
    
    
def find_sorted_local_input_fileset(stringmatch):
    '''Given a matching string (e.g. '*.txt') the function will create a sorted list
    of all matches '''
    
    filechoices = glob.glob(stringmatch)
    
    # Number of matches
    nmatch = len(filechoices)
    
    # Sort using Natural sort (see function at top)
    sort_nicely(filechoices)
    for i in range (nmatch):
        print '(',i+1,'): ', filechoices[i]
    
    
    print 'Detected ', nmatch, ' potential inputfiles in this directory'    
    
    
    return filechoices


def find_sorted_local_input_file(stringmatch):
    '''Given a matching string (e.g. '*.txt') the function allows user to select a single file'''
    
    filechoices = glob.glob(stringmatch)
    
    # Number of matches
    nmatch = len(filechoices)
    
    # Sort using Natural sort (see function at top)
    sort_nicely(filechoices)

    print 'Detected ', nmatch, ' potential inputfiles in this directory'
    print 'Here are the options: '

    for i in range (nmatch):
        print '(',i+1,'): ', filechoices[i]
        
    if(nmatch>0): print 'If none of these files suit:'
    print '(',nmatch+1,'):  Manually enter a filename'
    
    userselect = input('Make a selection: ')
    
    if userselect==nmatch+1:
        filename = raw_input('Manually enter filename: ')
    else:
        filename = filechoices[userselect-1]
        
    print 'File ',filename, ' selected for read-in'
    return filename 



def read_logfile(filename):
    '''Reads log file for a given star'''
    return np.genfromtxt(filename)

def read_logfiles():
    '''Reads all star files (matching to *.log)'''
    
    logfiles = find_sorted_local_input_fileset("*.log")
    
    alldata = []
    
    for filename in logfiles:
        print "Reading file ",filename
        alldata.append(read_logfile(filename))
        
    return alldata
    
    

def read_graph_file(filename):
    '''Reads the graph file, returning vertices, edges etc'''
    
    print 'Reading graph file ',filename

    data = np.genfromtxt(filename, skiprows = 1)

    vertexID = data[:,0]
    vertexPositions = data[:,1:4]
    edges = data[:,4:]

    nvertices = vertexID.shape[0]
    nedges = np.zeros(nvertices)
    
    # List Vertex Data
    for i in range(nvertices):
        nedges[i] = int(np.count_nonzero(edges[i,:]))
        print "Vertex ",i, "Edge Count",  int(nedges[i]), " position: ", vertexPositions[i,:]

    print "Total Number of edges: ", int(np.sum(nedges))
    print "Total Edge Length: ", np.sum(edges)
    
    return vertexID, vertexPositions, edges


def plot_graph(graphfile,outputfile=None,markerscale=0.1,vertexID=None,vertexPositions=None,edges=None):
    
    if(vertexID==None):
        vertexID, vertexPositions,edges = read_graph_file(graphfile)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.set_xlabel(xlabel,fontsize=axesize)
    ax.set_ylabel(ylabel,fontsize=axesize)
    
    ax.set_xlim((-1.5*np.amax(vertexPositions[:,0]),1.5*np.amax(vertexPositions[:,0])))
    ax.set_ylim((-1.5*np.amax(vertexPositions[:,1]),1.5*np.amax(vertexPositions[:,1])))
    
    # Plot edges first
        
        
    for i in range(len(edges)):
        
        for j in range(len(vertexID)):
            
            # If edge exists, draw line connecting vertices
            if(edges[i,j]>0.0):
        
                xdata = [vertexPositions[i,0],vertexPositions[j,0]]
                ydata = [vertexPositions[i,1],vertexPositions[j,1]]
                ax.plot(xdata,ydata, color='black',zorder=-1)
    
    
    ax.scatter(vertexPositions[:,0],vertexPositions[:,1],color='black')        
    
    if(outputfile==None):
        plt.show()
    else:
        fig.savefig(outputfile)
        
    plt.close()


def plot_graph_population(index,graphfile,outputfile=None, data=None,  markerscale=0.1,vertexID=None,vertexPositions=None,edges=None):
    '''Plots the graph at time index `index', with each star's population represented as a pie chart'''
    
    if(vertexID==None):
        vertexID, vertexPositions,edges = read_graph_file(graphfile)
    
    if(data==None):
        data = read_logfiles()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.set_xlabel(xlabel,fontsize=axesize)
    ax.set_ylabel(ylabel,fontsize=axesize)
    
    ax.set_xlim((-1.5*np.amax(vertexPositions[:,0]),1.5*np.amax(vertexPositions[:,0])))
    ax.set_ylim((-1.5*np.amax(vertexPositions[:,1]),1.5*np.amax(vertexPositions[:,1])))
    
    # Plot edges first
        
        
    for i in range(len(edges)):
        
        for j in range(len(vertexID)):
            
            # If edge exists, draw line connecting vertices
            if(edges[i,j]>0.0):
        
                xdata = [vertexPositions[i,0],vertexPositions[j,0]]
                ydata = [vertexPositions[i,1],vertexPositions[j,1]]
                ax.plot(xdata,ydata, color='black',zorder=-1)
    
    
    # Plot vertices (with circles representing predator/prey populations)
    
    for i in range(len(vertexID)):
        sizes = [data[i][index,iprey],data[i][index,ipred]]
        
        total = np.sum(sizes)
        
        if(total>1.0e-30):
            pred2prey = sizes[1]/total
        
            circle1 = plt.Circle(vertexPositions[i,0:2],radius = markerscale*total, color=preycolor,zorder=1)
            circle2 = plt.Circle(vertexPositions[i,0:2],radius = markerscale*pred2prey*total, color=predcolor,zorder=1)
            ax.add_artist(circle1)
            ax.add_artist(circle2)
        else:
            circle1 = plt.Circle(vertexPositions[i,0:2],radius = markerscale, color='gray',zorder=1)
            ax.add_artist(circle1)
            
    
    if(outputfile==None):
        plt.show()
    else:
        fig.savefig(outputfile)
        
    plt.close()
    

def plot_graph_population_movie(graphfile,interval=1,markerscale=0.1):
    '''Makes a movie of the system from multiple snapshots'''
    
    vertexID, vertexPositions,edges = read_graph_file(graphfile)
    data = read_logfiles()
     
    nsnap = len(data[0][:,0])
    
    nzeroes = int(np.log10(nsnap))+2
     
    print "There are ",nsnap, " snapshots"
    for i in range(0,nsnap,interval):
    
        outputfile = "snapshot."+str(i+1).zfill(nzeroes)+".png"
        
        print "Generating file ",outputfile
        plot_graph_population(i,graphfile,outputfile=outputfile, data=data,  markerscale=0.1,vertexID=vertexID,vertexPositions=vertexPositions,edges=edges)
    
def plot_population_star():

    filename = find_sorted_local_input_file('*.log')

    data=read_logfile(filename)
        
    time = data[:,itime]
    nprey = data[:,iprey]
    npred = data[:,ipred]
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel(timelabel,fontsize=axesize)
    ax1.set_ylabel(numlabel,fontsize=axesize)
    ax1.plot(time,nprey,label=preyleg,color=preycolor)
    ax1.plot(time,npred,label=predleg,color=predcolor)
    ax1.legend()

    ax1.set_yscale('log')
    plt.show()
    outputfile = 'pp_'+filename+'.png'

    fig.savefig(outputfile)
    
def plot_periodogram_star():
    
    filename = find_sorted_local_input_file('*.log')

    data=read_logfile(filename)
        
    time = data[:,itime]
    nprey = data[:,iprey]
    npred = data[:,ipred]
    
    print 'Computing periodogram'

    npoints = len(time)
    # Compute the periodogram of y along the x axis    

    dt = time[-1]/npoints

    # First, smooth and apply a window function

    w = blackman(npoints)

    
    freqprey, periodprey = periodogram(w*nprey, fs = 1.0/dt) 
    freqpred, periodpred = periodogram(w*npred, fs = 1.0/dt)

    # Calculate periods from these frequencies
    periods = np.divide(1.0,freqprey)

    
    prey_outputfile = 'prey_period'+filename+'.png'
    pred_outputfile = 'pred_period'+filename+'.png'

    # Make figures

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    
    ax1.set_xlabel("Period (years)", fontsize = 16)
    ax1.set_ylabel("Periodogram for Prey", fontsize = 16)
    
    print "Largest prey peak: ",periods[np.argmax(periodprey[periods<100.0])]

    ax1.plot(periods,periodprey)

    fig1.savefig(prey_outputfile, format='png')
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    
    ax2.set_xlabel("Period (years)", fontsize = 16)
    ax2.set_ylabel("Periodogram for Predators", fontsize = 16)

    ax2.plot(periods,periodpred)
    
    print "Largest predator peak: ",periods[np.argmax(periodpred[periods<100.0])]  

    fig2.savefig(pred_outputfile, format='png')
    
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.set_xlabel("Period (years)", fontsize = 16)
    ax3.set_ylabel("Periodogram Ratio", fontsize = 16)

    ax3.plot(periods,periodpred*periodprey)

    fig3.savefig('periodratio_'+filename+'.png', format='png')
    
    plt.show()
            
    
    
    
    
def plot_all_populations(data=None):
    if(data==None):
        data = read_logfiles()
    
    preyfig = plt.figure()
    ax1 = preyfig.add_subplot(111)
    
    predfig = plt.figure()
    ax2 = predfig.add_subplot(111)
    
    ax1.set_xlabel(timelabel,fontsize=axesize)
    ax1.set_ylabel(numlabel,fontsize=axesize)
    ax1.set_yscale('log')
        
    ax2.set_xlabel(timelabel,fontsize=axesize)
    ax2.set_ylabel(numlabel,fontsize=axesize)
    ax2.set_yscale('log')
    
    for i in range(len(data)):
        
        time = data[i][:,itime]
        nprey = data[i][:,iprey]
        npred = data[i][:,ipred]
        
        ax1.plot(time,nprey,color=preycolor)    
        ax2.plot(time,npred,color=predcolor)
        
    preyfig.savefig("allprey_vs_t.png")
    predfig.savefig("allpred_vs_t.png")
     
    
def calculate_total_populations(data=None):
    
    if(data==None):
        data = read_logfiles()
    
    for i in range(len(data)):
        if(i==0):
            time = data[i][:,itime]
            nprey = data[i][:,iprey]
            npred = data[i][:,ipred]
        else:
            nprey = nprey + data[i][:,iprey]
            npred = npred + data[i][:,ipred]

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(time,nprey,label=preyleg, color=preycolor)
    ax1.plot(time,npred,label=predleg,color=predcolor)
    ax1.set_xlabel(timelabel,fontsize=axesize)
    ax1.set_ylabel(numlabel,fontsize=axesize)
    ax1.set_yscale('log')
    ax1.legend()
    
    fig1.savefig("n_vs_t.png")
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(nprey,npred)
    ax2.set_xlabel(preylabel)
    ax2.set_ylabel(predlabel)
    fig2.savefig("prey_vs_predators.png")
    
    return nprey, npred


def calculate_total_periodogram(data=None):
    
    if(data==None):
        data = read_logfiles()
    
    for i in range(len(data)):
        if(i==0):
            time = data[i][:,itime]
            nprey = data[i][:,iprey]
            npred = data[i][:,ipred]
        else:
            nprey = nprey + data[i][:,iprey]
            npred = npred + data[i][:,ipred]


    npoints = len(time)
    # Compute the periodogram of y along the x axis    

    dt = time[-1]/npoints

    # First, smooth and apply a window function

    w = blackman(npoints)
    
    freqprey, periodprey = periodogram(w*nprey, fs = 1.0/dt) 
    freqpred, periodpred = periodogram(w*npred, fs = 1.0/dt)

    # Calculate periods from these frequencies
    periods = np.divide(1.0,freqprey)

    prey_outputfile = 'prey_period.png'
    pred_outputfile = 'pred_period.png'

    # Make figures

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    
    ax1.set_xlabel("Period (years)", fontsize = 16)
    ax1.set_ylabel("Periodogram for Prey", fontsize = 16)
    
    print "Largest prey peak: ",periods[np.argmax(periodprey[periods<100.0])]

    ax1.plot(periods,periodprey)

    fig1.savefig(prey_outputfile, format='png')
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    
    ax2.set_xlabel("Period (years)", fontsize = 16)
    ax2.set_ylabel("Periodogram for Predators", fontsize = 16)

    ax2.plot(periods,periodpred)
    
    print "Largest predator peak: ",periods[np.argmax(periodpred[periods<100.0])]  

    fig2.savefig(pred_outputfile, format='png')
    
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.set_yscale('log')
    ax3.set_xscale('log')
    ax3.set_xlabel("Period (years)", fontsize = 16)
    ax3.set_ylabel("Periodogram Ratio", fontsize = 16)

    ax3.plot(periods,periodpred*periodprey)

    fig3.savefig('periodratio.png', format='png')
    
    plt.show()
            

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(time,nprey,label=preyleg, color=preycolor)
    ax1.plot(time,npred,label=predleg,color=predcolor)
    ax1.set_xlabel(timelabel,fontsize=axesize)
    ax1.set_ylabel(numlabel,fontsize=axesize)
    ax1.set_yscale('log')
    ax1.legend()
    
    fig1.savefig("n_vs_t.png")
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(nprey,npred)
    ax2.set_xlabel(preylabel)
    ax2.set_ylabel(predlabel)
    fig2.savefig("prey_vs_predators.png")
    
    return nprey, npred
