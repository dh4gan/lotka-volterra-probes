import io_probe as io
import filefinder as ff
# Written 23/10/17 by dh4gan

# Plots the system at a specific snapshot

graphfile = ff.find_local_input_files("*.graph")

io.plot_graph(graphfile)