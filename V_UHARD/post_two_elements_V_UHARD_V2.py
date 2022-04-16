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
def post_two_elements(filename):
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
    instance = 'TENS_1-1'
    myInst   = myAssem.instances[instance]
    #-------------------------------------------------------------------------------
    # Extract data from the element
    #-------------------------------------------------------------------------------
    for i,element in enumerate(myInst.elementSets['ELEMENTS'].elements):
        ELEMENTHist = step.historyRegions['Element '+instance+'.%s'%element.label+' Int Point 1']
        eps_rate = np.array(ELEMENTHist.historyOutputs['LE22'].data)[:,1]
        sig_rate = np.array(ELEMENTHist.historyOutputs['S22'].data)[:,1]
        tmp_rate = np.abs(np.array(ELEMENTHist.historyOutputs['TEMP'].data)[:,1])
        time = np.array(ELEMENTHist.historyOutputs['S22'].data)[:,0]
    #-------------------------------------------------------------------------------
    # Load and define the instance
    #-------------------------------------------------------------------------------
    instance = 'TENS_2-1'
    myInst   = myAssem.instances[instance]
    #-------------------------------------------------------------------------------
    # Extract data from the element
    #-------------------------------------------------------------------------------
    for i,element in enumerate(myInst.elementSets['ELEMENTS'].elements):
        ELEMENTHist = step.historyRegions['Element '+instance+'.%s'%element.label+' Int Point 1']
        eps_ref = np.abs(np.array(ELEMENTHist.historyOutputs['LE22'].data)[:,1])
        sig_ref = np.abs(np.array(ELEMENTHist.historyOutputs['S22'].data)[:,1])
        tmp_ref = np.abs(np.array(ELEMENTHist.historyOutputs['TEMP'].data)[:,1])
    #-------------------------------------------------------------------------------
    # Close odb file
    #-------------------------------------------------------------------------------
    odb.close()
    data = np.transpose(np.array([time,eps_rate,sig_rate,tmp_rate,eps_ref,sig_ref,tmp_ref]))
    keys = ['time','eps_rate','sig_rate','tmp_rate','eps_ref','sig_ref','tmp_ref']
    write_results(filename,data,keys)
    return
####################################################################################
####################################################################################
# START OF SCRIPT
####################################################################################
####################################################################################
filename = sys.argv[-1]
post_two_elements(filename)
exit()
####################################################################################
####################################################################################
# END OF SCRIPT
####################################################################################
####################################################################################