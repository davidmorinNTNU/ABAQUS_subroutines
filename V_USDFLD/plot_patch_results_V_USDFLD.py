import numpy as np
import matplotlib.pyplot as plt
#
def read_csv_file(filename):
    """
    Read the csv file from the ABAQUS post-processing script
    """
    # Open text file for reading
    fp = open(filename+'.csv','r')
    # Scroll through the lines in the file
    for line in fp:
        if line[0:1] == '*':
           # Extract the variable names from the first line of the file
           keys = [x.strip() for x in line.replace('*','').split(',')]
           # Create a dictionnary with a list for each variable/key
           dico = {}
           for key in keys:
               dico[key] = []
        else:
           # Store each value in the corresponding dictionnary key
           for key,value in zip(keys,line.split(',')):
               dico[key].append(float(value))
    # Scroll through the lines in the file
    fp.close()
    # Change the lists to arrays
    for key in list(dico.keys()):
        dico[key] = np.array(dico[key])
    return dico
##########################################################################################
# START OF THE SCRIPT
##########################################################################################
# Define files to be read
filenames = ['example_VUSDFLD_V2_data','example_USDFLD_V2_data']
# Read csv files
results = {}
for filename in filenames:
    results[filename] = read_csv_file(filename)
# Plot results
plt.rcParams["figure.figsize"] = (5, 4)
for key in list(results.keys()):
    plt.plot(results[key]['disp'],results[key]['force']/1000.0,label=key.split('_')[1]+' force')
plt.xlabel('Displacement (mm)')
plt.ylabel('Force (kN)')
plt.xlim([0.0,1.0])
plt.ylim([0.0,0.8])
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('USDFLD_VUSDFLD_V2.pdf')
plt.show()
##########################################################################################
# END OF THE SCRIPT
##########################################################################################
exit()