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
filenames = ['example_VUHARD_V2_data']
# Read csv files
results = {}
for filename in filenames:
    results[filename] = read_csv_file(filename)
# Plot results
plt.rcParams["figure.figsize"] = (8, 4)
fig,axs = plt.subplots(1,2)
k = 0
for key in list(results.keys()):
    axs[0].plot(results[key]['eps_rate'],results[key]['sig_rate'],label=key.split('_')[1]+' Thermo-Visco-plastic')
    axs[0].plot(results[key]['eps_ref'],results[key]['sig_ref'],label=key.split('_')[1]+' Visco-Plastic')
    axs[0].set_xlabel('True strain (-)')
    axs[0].set_ylabel('True stress (MPa)')
    axs[0].set_xlim([0.0,0.4])
    axs[0].set_ylim([0.0,1500.0])
    axs[0].legend()
    axs[0].grid()
#
    axs[1].plot(results[key]['eps_rate'],results[key]['tmp_rate']-273.0,label=key.split('_')[1]+' Thermo-Visco-plastic')
    axs[1].plot(results[key]['eps_ref'],results[key]['tmp_ref']-273.0,label=key.split('_')[1]+' Visco-Plastic')
    axs[1].set_xlabel('True strain (-)')
    axs[1].set_ylabel('Temperature (Degree C)')
    axs[1].set_xlim([0.0,0.4])
    axs[1].set_ylim([0.0,200.0])
    axs[1].legend()
    axs[1].grid()
    k += 1
plt.tight_layout()
plt.savefig('VUHARD_V2.pdf')
plt.show()
##########################################################################################
# END OF THE SCRIPT
##########################################################################################
exit()