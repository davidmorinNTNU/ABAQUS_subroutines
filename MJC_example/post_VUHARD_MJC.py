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
def post_round_specimen(filename):
#-------------------------------------------------------------------------------
#   Open odb file
#-------------------------------------------------------------------------------
    odb  = openOdb(path=filename+'.odb')
#-------------------------------------------------------------------------------
#   Load the step
#-------------------------------------------------------------------------------
    stepname = 'LOADING'
    step = odb.steps[stepname]
#-------------------------------------------------------------------------------
#   Load Assembly and Instances
#-------------------------------------------------------------------------------
    myAssem    = odb.rootAssembly
    myInst     = myAssem.instances
    part       = 'R2'
    instance   = part+'-1'
    myInst     = myInst[instance]
#-----------------------------------------------------------------------
#   Extract displacement from diameter
#-----------------------------------------------------------------------
    for i,node in enumerate(myInst.nodeSets['DIAMETER'].nodes):
        radius = node.coordinates[0]
        diamHist = step.historyRegions['Node {}-1.{}'.format(part,node.label)]
        diam_red =-2.0*np.array(diamHist.historyOutputs['U1'].data)[:,1]
#-----------------------------------------------------------------------
#   Extract force from reference point
#-----------------------------------------------------------------------
    for i,node in enumerate(myAssem.nodeSets[part].nodes):
        forceHist = step.historyRegions['Node ASSEMBLY.%s'%node[0].label]
        force     = np.array(forceHist.historyOutputs['RF2'].data)[:,1]/1000.0
        time      = np.array(forceHist.historyOutputs['RF2'].data)[:,0]
    odb.close()
#-----------------------------------------------------------------------
#   Extract force from reference point
#-----------------------------------------------------------------------
    data = np.transpose(np.array([time,diam_red,force]))
    keys = ['time','diam_red','force']
    write_results(filename,data,keys)
    return
####################################################################################
####################################################################################
# START OF SCRIPT
####################################################################################
####################################################################################
filename = sys.argv[-1]
post_round_specimen(filename)
exit()
####################################################################################
####################################################################################
# END OF SCRIPT
####################################################################################
####################################################################################