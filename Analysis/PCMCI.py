#!/usr/bin/env python
# coding: utf-8

# In[16]:    


# Modules and data
import pandas as pd
import os
import pyarrow.feather as feather
from matplotlib import pyplot as plt
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr
import numpy as np

data_all = pd.read_csv("chinaflu_norm.csv")
data_all.rename(columns={'fluP': 'Flu',
                         'o3': '$\mathregular{O_3}$',
                         'temp': 'T',
                         'pm2.5': '$\mathregular{PM_{2.5}}$',
                         'ah': 'AH'}, inplace=True)


# Function for preparing the values of the variables specified by 'var_names'
# in the province specified by 'state' as numpy array of shape (T, N), where
# T is the number of time steps and N the number of variables
def get_data(data_all, state, var_names):
    
    # Select the state
    if state == "Overall":
        data_out = data_all
    else:
        data_out = data_all.loc[data_all[r'state'] == state]

    # Select the columns
    data_out = data_out[var_names]

    # Turn into numpy array
    data_out = data_out.values

    # Return
    return data_out

def get_selected_links(var_names, tau_min, tau_max):

    # Get index of the temperature variable, if it exists
    if 'T' in var_names:
        temp_idx = np.argwhere(np.array(var_names) == 'T')[0, 0]
    else:
        temp_idx = None

    # Get index of the humidity variable, if it exists
    if 'AH' in var_names:
        humid_idx = np.argwhere(np.array(var_names) == 'AH')[0, 0]
    else:
        humid_idx = None

    # Build dictionary
    selected_links = {}
    
    for idx, var in enumerate(var_names):

        if var == 'Flu':
            # Flu may be influenced by all variables at lags
            selected_links[idx] = [(other_idx, -tau) for other_idx, other_var in enumerate(var_names)
                                   for tau in range(tau_min, tau_max + 1) ]
        elif var == 'AH':
            # Humiditiy may be influenced by itself at all lags
            selected_links[idx] = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]

            # Humidity may also be influenced by temperature
            selected_links[idx] = [(temp_idx, -tau) for tau in range(tau_min, tau_max + 1)]

        elif var == 'T':
            # Temperature may be influenced by itself at all lags
            selected_links[idx] = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]

            # Temperature may also be influenced by humidity
            selected_links[idx] = [(humid_idx, -tau) for tau in range(tau_min, tau_max + 1)]

        else:
            # All other variables, here this are O3 and PM2.5, may be influenced 
            # by all variables other than Flu 
            selected_links[idx] = [(other_idx, -tau) for other_idx, other_var in enumerate(var_names)
                                   for tau in range(tau_min, tau_max + 1) if
                                   (other_var != 'Flu')]

    # Return
    return selected_links

def apply_pcmci(data_all,
                state,
                var_names,
                tau_min,
                tau_max,
                pc_alpha,
                verbosity):
    
    # Get the data and mask
    data = get_data(data_all=data_all,
                                   state=state,
                                   var_names=var_names)

    # Prepare the DataFrame object
    dataframe = pp.DataFrame(data,
                             var_names=var_names,
                             missing_flag=999.)

    # Prepare the independence test and PCMCI object
    parcorr = ParCorr()
    pcmci = PCMCI(dataframe=dataframe,
                  cond_ind_test=parcorr,
                  verbosity=verbosity)

    # Get the selected_links arguement
    selected_links = get_selected_links(var_names,
                                        tau_min,
                                        tau_max)

    # Run PCMCI with these parameters
    results = pcmci.run_pcmciplus(tau_min=tau_min,
                                  tau_max=tau_max,
                                  pc_alpha=pc_alpha,
                                  selected_links=selected_links)
    # Plot
    tp.plot_graph(
        arrow_linewidth=8.0,
        figsize=(10*0.5, 10*0.5),
        vmin_edges=-0.5,
        vmax_edges=0.5,
        node_label_size=16,
        node_size=0.3,
        link_label_fontsize=13,
        val_matrix=results['val_matrix'],
        graph=results['graph'],
        var_names=var_names,
        link_colorbar_label='cross-MCI (edges)',
        node_colorbar_label='auto-MCI (nodes)',
        label_fontsize=15,
        network_lower_bound=0.2,
        show_colorbar=1
    );
    
    plt.show()

    return results


# In[17]:


# Province-wise analysis
states = ["Anhui", "Beijing", "Chongqing", "Fujian", "Gansu", "Guangdong", "Guangxi", "Guizhou", 
          "Hainan", "Hebei", "Heilongjiang", "Henan", "Hubei", "Hunan", "Inner Mongolia", "Jiangsu", 
          "Jiangxi", "Jilin", "Liaoning", "Ningxia", "Qinghai", "Shaanxi", "Shandong", "Shanghai", 
          "Shanxi", "Sichuan", "Tianjin", "Xinjiang", "Yunnan", "Zhejiang"]

pc_alpha = 0.05 

for state in states:
    results = apply_pcmci(data_all=data_all,
                      state=state,
                      var_names=['Flu', '$\mathregular{O_3}$', 'T', 'AH','$\mathregular{PM_{2.5}}$'],
                      tau_min=0,
                      tau_max=2,
                      pc_alpha=pc_alpha,
                      verbosity=0
                      )

    


# In[18]:


# Nationwide analysis
pc_alpha = 0.001

results = apply_pcmci(data_all=data_all,
                      state="Overall",
                      var_names=['Flu', '$\mathregular{O_3}$', 'T', 'AH','$\mathregular{PM_{2.5}}$'],
                      tau_min=0,
                      tau_max=2,
                      pc_alpha=pc_alpha,
                      verbosity=0
                      )
print(pc_alpha)


# In[ ]:




