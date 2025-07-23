# mesa_net_reader

A reader (and simple plotter) of formatted nuclear network files from [MESA](https://github.com/MESAHub/mesa).

### Requires:
numpy, matplotlib, pandas

### Basic usage:
```
import nuclearNetAnalyzer as nuc_an  
import nuclear_aux as nuc_aux  
data_dict = nuc_an.analyzeNetwork('example/mesa_80', mesa=False, netdir=os.getcwd())  
        # data_dict : dict ; Dictionary containing the elements (keys) and isotopes (values).  
data, elementName, eleLabelPos = nuc_an.getDataStructs(data_dict,  
                                                        periodicTable=pd.read_csv('sources/periodicTableNames.csv'),  
                                                        fullIsotopeList=pd.read_csv('sources/nndc_nudat_data_export.csv'))  
        # data : np.array ; Array containing the isotopes data as a grid.  
        # elementName : list ; List of element names corresponding to the isotopes.  
        # eleLabelPos : list ; List of positions for the element labels.  
plt.imshow(data, origin="lower", vmin=0)  
for j,iso in enumerate(np.array(elementName)):  
    plt.text(eleLabelPos[j][1],eleLabelPos[j][0],iso,ha='center',va='center',fontsize=6)  
```
#### [[[See example in example/mesa_80_readplot.ipynb]]]

### Structure:  
- mesa_net_reader/  
    - nuclear_aux.py                  # parsing functions  
    - nuclearNetAnalyzer.py           # main call, data init  
    - sources/                        # data sources and csv files  
        - nndc_nudat_data_export.csv  # isotope database from [NuDat](https://www.nndc.bnl.gov/nudat/)  
        - periodicTableNames.csv      # name to number convert  
    - example/                        # read-and-plot example for mesa_80 network  
        - network_read_plot.ipynb       
        - mesa_80.net                 # network example from https://github.com/MESAHub/mesa/tree/main/data/net_data/nets/mesa.80.net  
