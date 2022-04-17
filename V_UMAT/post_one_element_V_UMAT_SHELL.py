from odbAccess import *
import sys
import numpy as np
#
def write_results(filename,data,keys):
    header = '*'+','.join(keys)+'\n'
    nvars,ntime  = np.shape(data)[1],np.shape(data)[0]
    format_to_write = '{},'*(nvars-1)+'{}\n'
    fp = open(filename+'_data.csv','w')
    fp.write(header)
    for i in range(0,ntime):
        fp.write(format_to_write.format(*[x for x in data[i,:]]))
    fp.close()
    return
#
def post_one_element(filename):
    #-------------------------------------------------------------------------------
    # Open odb file
    #-------------------------------------------------------------------------------
    odb  = openOdb(path=filename+'.odb')
    #-------------------------------------------------------------------------------
    # Load the step
    #-------------------------------------------------------------------------------
    stepname = 'LOADING'
    step = odb.steps[stepname]
    #-------------------------------------------------------------------------------
    # Load the assembly
    #-------------------------------------------------------------------------------
    myAssem = odb.rootAssembly
    #-------------------------------------------------------------------------------
    # Load and define the instance
    #-------------------------------------------------------------------------------
    myInst   = myAssem.instances
    instance = myAssem.instances.keys()[0]
    myInst   = myInst[instance]
    #-------------------------------------------------------------------------------
    # Extract data from the element
    #-------------------------------------------------------------------------------
    print(step.historyRegions)
    for i,element in enumerate(myInst.elementSets['ELEMENTS'].elements):
        ELEMENTHist = step.historyRegions['Element '+instance+'.%s'%element.label+' Int Point 1 Section Point 1']
        LE22   = np.array(ELEMENTHist.historyOutputs['LE22'].data)[:,1]
        S11    = np.array(ELEMENTHist.historyOutputs['S11'].data)[:,1]
        S22    = np.array(ELEMENTHist.historyOutputs['S22'].data)[:,1]
        S12    = np.array(ELEMENTHist.historyOutputs['S12'].data)[:,1]
        time   = np.array(ELEMENTHist.historyOutputs['LE22'].data)[:,0]
        ELEMENTHist = step.historyRegions['Element '+instance+'.%s'%element.label+' Int Point 1']
        thick  = np.array(ELEMENTHist.historyOutputs['STH'].data)[:,1]
    #-------------------------------------------------------------------------------
    # Close odb file
    #-------------------------------------------------------------------------------
    odb.close()
    data = np.transpose(np.array([time,LE22,S11,S22,S12,thick]))
    keys = ['time','true_strain','S11','S22','S12','thick']
    write_results(filename,data,keys)
    return
####################################################################################
####################################################################################
# START OF SCRIPT
####################################################################################
####################################################################################
filename = sys.argv[-1]
post_one_element(filename)
exit()
####################################################################################
####################################################################################
# END OF SCRIPT
####################################################################################
####################################################################################