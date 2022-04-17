from abaqus import *
from abaqusConstants import *
viewport = session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=169.068740844727, height=128.548233032227)
viewport.makeCurrent()
viewport.maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()

def post_pictures(filename):
    o1 = session.openOdb(name=filename+'.odb')
    viewport.setValues(displayedObject=o1)
    viewport.viewportAnnotationOptions.setValues(triad=OFF, title=OFF, state=OFF, annotations=OFF, compass=OFF)
    viewport.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    viewport.odbDisplay.setPrimaryVariable(variableLabel='TEMP', outputPosition=INTEGRATION_POINT, )
    session.viewports[session.currentViewportName].odbDisplay.setFrame(step='LOADING', frame=81)
    viewport.view.fitView()
    viewport.odbDisplay.commonOptions.setValues(visibleEdges=FREE)
    session.printOptions.setValues(reduceColors=False)
    session.printToFile(fileName=filename+'_TEMP_R2', format=PNG, canvasObjects=(viewport, ))
    viewport.odbDisplay.setPrimaryVariable(variableLabel='SDV_D', outputPosition=INTEGRATION_POINT, )
    session.printToFile(fileName=filename+'_D_R2', format=PNG, canvasObjects=(viewport, ))
    viewport.view.setValues(nearPlane=40.5815, farPlane=61.547, width=6.5314, height=2.96063, viewOffsetX=0.140359, viewOffsetY=-6.26133)
    session.printToFile(fileName=filename+'_D_ZOOM_R2', format=PNG, canvasObjects=(viewport, ))
    viewport.odbDisplay.setPrimaryVariable(variableLabel='TEMP', outputPosition=INTEGRATION_POINT, )
    session.printToFile(fileName=filename+'_TEMP_ZOOM_R2', format=PNG, canvasObjects=(viewport, ))
    session.odbs[filename+'.odb'].close()
    return

filename='example_VUHARD_MJC'
post_pictures(filename)