import io_probe as io
import matplotlib.pyplot as plt


finished = False

while(not(finished)):
    fig = io.plot_periodogram_star()
    savechoice = raw_input('Try another file? (y/n) ')

    if(savechoice=='y' or savechoice=='Y'):
        finished = False
    else:
        finished = True




