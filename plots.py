from tkinter.tix import ButtonBox
from typing import Iterable
import matplotlib.pyplot as plt
import numpy as np
from pyparsing import trace_parse_action
# from constants import constants



def temps_tube_fluid(temps: dict, row_idx: Iterable[int], fluid: str) -> plt.figure:
    """Plots the temperature of a specific fluid along a tube

    Args:
        temps (dict): temperature data
        row_idx (Iterable): list of indeces of row to plot
        fluid (str): name of the fluid to plot

    Returns:
        plt.figure: figure
    """    

    # create segments vector
    # x = range(constants.n_segments) # need to fix constants first
    x = range(100)

    dict_key = 't_'+fluid

    fig = plt.figure()
    ax = fig.add_subplot(111) 
    lines = []
    rows = []
    for i in range(len(row_idx)):
        line, = ax.plot(x,temps[dict_key][:,i])
        tube_label = 'row '+str(i)
        lines.append(line)
        rows.append(tube_label)

    ax.set_title(f'Temperature of {fluid} along tube')
    ax.set_xlabel('Segment')
    ax.set_ylabel('Temperature [k]')
    ax.legend(lines, rows)

    plt.savefig(f't_{fluid}_tube.png') 
    return fig


def temps_tube_both(temps: dict, row_idx: int) -> plt.figure:
    """Plot co2 and air along one tube

    Args:
        temps (dict): temperature data
        tube_idx (int): index of row to plot

    Returns:
        plt.figure: figure
    """      

    x = range(100)

    fig = plt.figure()
    ax = fig.add_subplot(111)   
    air, = ax.plot(x,temps['t_air'][:,row_idx])
    co2, = ax.plot(x,temps['t_co2'][:,row_idx])
    
    ax.set_title('Temperatures along tube of both fluids')
    ax.set_xlabel('Segment')
    ax.set_ylabel('Temperature [K]')

    ax.legend([air,co2],['air','co2'])

    plt.savefig('tube_co2_air.png')

    return fig



def bundle_temps(temps: dict) -> plt.figure:
    """Plots temperature distribution across the entire bundle

    Args:
        temps (dict): temperature data

    Returns:
        plt.figure: figure
    """    

    n_rows = 4
    n_segments = 100

    temps_array = np.empty((n_segments,0))
    for i in range(n_rows):
        temps_array = np.c_[temps_array, temps["t_air"][:,i]]
        temps_array = np.c_[temps_array, temps["t_co2"][:,i]]
        
    fig = plt.figure()
    ax = fig.add_subplot(111)   
    im = ax.imshow(temps_array, cmap='hot', interpolation='nearest')
    ax.set_title('Heat exchanger temperature distribution')
    ax.set_xlabel('Tube number')
    ax.set_ylabel('Segment')
    fig.colorbar(im)

    plt.savefig('heat_map.png')

    return fig



if __name__ == "__main__":
    ## generate test data (TODO: change to real simulator data)
    test = np.linspace(0,1,100)
    t_co2_1 = np.sin(test)
    t_air_1 = np.sin(test)+0.1

    t_co2_2 = t_co2_1+0.1
    t_air_2 = t_air_1+0.1

    t_co2_3 = t_co2_2+0.1
    t_air_3 = t_air_2+0.1

    t_co2_4 = t_co2_3+0.1
    t_air_4 = t_air_3+0.1

    temps = {}
    temps["t_co2"] = np.vstack([t_co2_1,t_co2_2,t_co2_3,t_co2_4]).T
    temps["t_air"] = np.vstack([t_air_1,t_air_2,t_air_3,t_air_4]).T


    ##### 

    # plot co2 temp in one tube

    fig1 = temps_tube_fluid(temps=temps, row_idx=[0], fluid='co2')

    # plot air temp in all tubes

    fig2 = temps_tube_fluid(temps=temps, row_idx=[0,1,2,3], fluid='air')

    # plot co2 and air across a tube

    fig3 = temps_tube_both(temps=temps, row_idx=0)

    # plot all temps in a bundle

    fig4 = bundle_temps(temps=temps)






