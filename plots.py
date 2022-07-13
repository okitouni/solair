from tkinter.tix import ButtonBox
from typing import Iterable
import matplotlib.pyplot as plt
import numpy as np
from pyparsing import trace_parse_action
# from constants import constants



def temps_tube_fluid(temps: dict, tube_idx: Iterable, fluid: str):
    """Plots the temperature of a specific fluid along a tube

    Args:
        temps (dict): temperature data
        tube_idx (Iterable): list of indeces of tubes to plot
        fluid (str): name of the fluid to plot
    """    

    # create segments vector
    # x = range(constants.n_segments)
    x = range(100)

    dict_key = 't_'+fluid

    plt.figure()
    plt.plot(x,temps[dict_key][:,tube_idx])
    plt.title(f'Temperature of {fluid} along tube')
    plt.xlabel('Segment')
    plt.ylabel('Temperature [k]')
    
    plt.show()
    plt.savefig(f't_{fluid}_tube.png') 


def temps_tube_both(temps: dict, tube_idx: int):
    """Plot co2 and air along one tube

    Args:
        temps (dict): temperature data
        tube_idx (int): index of tube to plot
    """    

    x = range(100)

    plt.figure()
    air = plt.plot(x,temps['t_air'][:,tube_idx])
    co2 = plt.plot(x,temps['t_co2'][:,tube_idx])
    plt.title('Temperatures along tube of both fluids')
    plt.xlabel('Segment')
    plt.ylabel('Temperature [k]')
    # plt.legend([air,co2],['Air','CO2'])
    plt.show()
    plt.savefig('tube_co2_air.png')



def bundle_temps(temps: dict):
    """Plots temperature distribution across the entire bundle

    Args:
        temps (dict): temperature data
    """    

    n_tubes = 4
    n_segments = 100

    temps_array = np.empty((n_segments,0))
    for i in range(n_tubes):
        temps_array = np.c_[temps_array, temps["t_air"][:,i]]
        temps_array = np.c_[temps_array, temps["t_co2"][:,i]]
        
    plt.imshow(temps_array, cmap='hot', interpolation='nearest')
    plt.title('Heat exchanger temperature distribution')
    plt.xlabel('Tube number')
    plt.ylabel('Segment')
    plt.colorbar()

    plt.savefig('heat_map.png')



## generate test data
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

temps_tube_fluid(temps=temps, tube_idx=0, fluid='co2')

# plot air temp in all tubes

temps_tube_fluid(temps=temps, tube_idx=range(4), fluid='air')

# plot co2 and air across a tube

temps_tube_both(temps=temps, tube_idx=0)

# plot all temps in a bundle

bundle_temps(temps=temps)






