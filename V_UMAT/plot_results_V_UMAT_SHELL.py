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
filenames = ['example_VUMAT_SHELL_data']
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
    ax1.plot(results[key]['true_strain'],results[key]['S11'],label=r'$\sigma_{11}$')
    ax1.plot(results[key]['true_strain'],results[key]['S22'],label=r'$\sigma_{22}$')
    ax1.plot(results[key]['true_strain'],results[key]['S12'],label=r'$\sigma_{12}$')
    ax2.plot(results[key]['true_strain'],results[key]['thick'])
ax1.set_xlabel('True longitudinal strain (-)')
ax1.set_ylabel('Stress (MPa)')
ax2.set_xlabel('True longitudinal strain (-)')
ax2.set_ylabel('Thickness (mm)')
ax1.set_xlim([0.0,0.4])
ax1.set_ylim([0.0,400.0])
ax2.set_xlim([0.0,0.4])
ax2.set_ylim([0.0,0.101])
ax1.legend()
ax1.grid()
ax2.grid()
plt.tight_layout()
plt.savefig('VUMAT_SHELL.pdf')
plt.show()
##########################################################################################
# END OF THE SCRIPT
##########################################################################################
exit()