import io_probe as io
import filefinder as ff

# Written 23/10/17 by dh4gan

# Makes snapshots of the graph evolution for
# creating animations

graphfile = ff.find_local_input_files("*.graph")

#interval = input("What is the snapshot interval?")
interval = 1

io.plot_graph_population_movie(graphfile,interval=interval)