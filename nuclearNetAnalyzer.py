# -*- coding: utf-8 -*-
"""
@author: Dmitry Shishkin
"""
#%% Imports:

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import nuclear_aux as nuc_aux


#%%

'''
Add "show_net_species_info = .true." to inlist to print out a list, or use
    below to work with raw data from the mesa network files.
Will then plot the network on a grid.
'''

def analyzeNetwork(net, mesa=True, netdir=None):
    """Main analyzing function.
    
    Reads a nuclear network file and returns a dict
    
    Parameters
    ----------
    net : str
        .net file name to be read in.
    mesa : bool, optional
        Whether to assume MESA tree as the source or in dir. 
        Default is False, i.e, look for the $MESA_DIR enviromental variable and
        navigate into the standard network local. If True, prvoide...
    netdir : str, optional
        Directory where network file is located.
        
    Returns
    -------
    data_dict : dict
        Dictionary containing the elements (keys) and isotopes (values).
    """
    
    if mesa:
        mesa_dir = os.environ.get('MESA_DIR')
        file_path = os.path.join(mesa_dir,'data','net_data','nets')
    else:
        if netdir == None:
            file_path = os.getcwd()
        else:
            file_path = netdir
    
    data_dict = {}
    data_dict = nuc_aux.read_file(data_dict,file_path,net)
    if "" in data_dict.keys():
        del data_dict[""]
    print(data_dict)
    
    return data_dict | {"network":net}

# Call function to get the network data_dict 
# In the example: mesa_80.net, do not go to MESA folder - instead look for in cwd.
# data_dict = analyzeNetwork('example/mesa_80', mesa=False, netdir=os.getcwd())

#########################################################################
#%% Processing for the plotting part - building on the extract data_dict
#########################################################################

def getDataStructs(data_dict, periodicTable, fullIsotopeList, nrows=200, ncols=200, minProt=0, maxProt=200, minNeut=0, maxNeut=200):
    """Function to get the data structures for plotting.
    
    Parameters
    ----------
    periodicTable : pd.DataFrame, optional
        DataFrame containing the periodic table information. Default is read from
        'sources/periodicTableNames.csv'.
    fullIsotopeList : pd.DataFrame, optional
        DataFrame containing the full list of isotopes. Default is read from
        'sources/nndc_nudat_data_export.csv'.
    nrows : int, optional
        Number of rows in the grid. Default is 200.
    ncols : int, optional
        Number of columns in the grid. Default is 200.
        
    Returns
    -------
    data : np.array
        Array containing the isotopes data as a grid.
    elementName : list
        List of element names corresponding to the isotopes.
    eleLabelPos : list
        List of positions for the element labels.
    """

    elementName=[]
    eleLabelPos=[]

    # Isotopes are displayed in dark blue
    Iso_id=[]
    # Numbers are extracted from a csv : "periodicTableNames.csv"
    # Source: Self
    # periodicTable = pd.read_csv('sources/periodicTableNames.csv') 
    for element in list(data_dict.keys())[:-1:1]:
        y = np.array(periodicTable.AtomicNumber[periodicTable.Symbol==element.capitalize()])#prot
        x = np.array(data_dict[element])-y                                                  #neut
        if element == 'neut':
            y = np.array(0); x = np.array(1);
        elif element == 'prot':
            y = np.array(1); x = np.array(0);
        if element=='neut' or element=='prot':
            elementName.append(element[0])
        else:
            elementName.append(element.capitalize())
        eleLabelPos.append([y,np.min(x)-1])
        minProt = min(y,minProt); maxProt = max(y,maxProt);                                 #for limits
        minNeut = np.min([np.min(x),minNeut]); maxNeut = np.max([np.max(x),maxNeut]);       #for limits
        Iso_id.append(x+y*ncols)

    Iso_id = np.array(list(pd.core.common.flatten(Iso_id)))
    Iso_val = Iso_id*0+1
    data = np.zeros(nrows*ncols)
    data[Iso_id] = Iso_val
    data = np.ma.array(data.reshape((nrows, ncols)), mask=data==0)

    Back_id=[]
    # Values are extracted from a csv : "nndc_nudat_data_export.csv"
    # Source: National Nuclear Data Center, information extracted from the NuDat database, https://www.nndc.bnl.gov/nudat/
    # fullIsotopeList = pd.read_csv('sources/nndc_nudat_data_export.csv') 
    for i,pro in enumerate(fullIsotopeList.z):
        Back_id.append(fullIsotopeList.n[i]+pro*ncols)
        
    Back_id = np.array(list(pd.core.common.flatten(Back_id)))
    Back_val = Back_id*0+0.1
    dataB = np.zeros(nrows*ncols)
    dataB[Back_id] = Back_val
    dataB = np.ma.array(dataB.reshape((nrows, ncols)), mask=dataB==0)

    data[data.mask] = dataB[data.mask]

    return data.data, elementName, eleLabelPos


