import io_probe as io
import filefinder as ff
# Written 23/10/17 by dh4gan

# Plots the system at a specific snapshot

graphfile = ff.find_local_input_files("*.graph")

itime = input("Which snapshot is to be plotted? ")

outputfile = graphfile+".snap."+str(itime)

io.plot_graph_population(itime, graphfile)