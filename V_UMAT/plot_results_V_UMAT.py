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
filenames = ['example_UMAT_data','example_VUMAT_data']
# Read csv files
results = {}
for filename in filenames:
    results[filename] = read_csv_file(filename)
# Plot results
plt.rcParams["figure.figsize"] = (8, 4)
fig, axs = plt.subplots(1,2)
ax1 = axs[0]
ax2 = axs[1]
for key in list(results.keys()):
    ax1.plot(results[key]['disp'],results[key]['force']/1000.0,label=key.split('_')[1]+' force')
    ax2.plot(results[key]['disp'],results[key]['damage'],label=key.split('_')[1]+' damage')
ax1.set_xlabel('Displacement (mm)')
ax1.set_ylabel('Force (kN)')
ax2.set_xlabel('Displacement (mm)')
ax2.set_ylabel('Damage (-)')
ax1.set_xlim([0.0,2.0])
ax1.set_ylim([0.0,10.0])
ax2.set_xlim([0.0,2.0])
ax2.set_ylim([0.0,1.0])
ax1.legend()
ax2.legend()
ax1.grid()
ax2.grid()
plt.tight_layout()
plt.savefig('UMAT_VUMAT.pdf')
plt.show()
##########################################################################################
# END OF THE SCRIPT
##########################################################################################
exit()