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
    # Extract data from reference point
    #-------------------------------------------------------------------------------
    for i,node in enumerate(myAssem.nodeSets['REFPT_TOP'].nodes[0]):
        REFPTHist = step.historyRegions['Node ASSEMBLY.%s'%node.label]
        force     = np.array(REFPTHist.historyOutputs['RF2'].data)[:,1]
        disp      = np.array(REFPTHist.historyOutputs['U2'].data)[:,1]
        time      = np.array(REFPTHist.historyOutputs['U2'].data)[:,0] 
    #-------------------------------------------------------------------------------
    # Extract data from the element
    #-------------------------------------------------------------------------------
    for i,element in enumerate(myInst.elementSets['ELEMENTS'].elements):
        ELEMENTHist = step.historyRegions['Element '+instance+'.%s'%element.label+' Int Point 1']
        damage = np.array(ELEMENTHist.historyOutputs['SDV_D'].data)[:,1]
    #-------------------------------------------------------------------------------
    # Close odb file
    #-------------------------------------------------------------------------------
    odb.close()
    data = np.transpose(np.array([time,disp,force,damage]))
    keys = ['time','disp','force','damage']
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