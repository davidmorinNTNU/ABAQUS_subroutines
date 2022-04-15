import sys
from abaqus import *
from abaqusConstants import *
viewport = session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=169.068740844727, height=128.711563110352)
viewport.makeCurrent()
viewport.maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
#
def post_pictures(filename,frame):
    o1 = session.openOdb(name=filename+'.odb')
    viewport.setValues(displayedObject=o1)
    viewport.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    viewport.view.setValues(session.views['Front'])
    viewport.viewportAnnotationOptions.setValues(triad=OFF, title=OFF, state=OFF, annotations=OFF, compass=OFF)
    viewport.odbDisplay.setPrimaryVariable(variableLabel='SDV_D', outputPosition=INTEGRATION_POINT, )
    viewport.odbDisplay.setFrame(step='LOADING', frame=frame)
    viewport.odbDisplay.contourOptions.setValues(contourType=QUILT, maxValue=1.0, minValue=0.0)
    session.printOptions.setValues(reduceColors=False)
    session.printToFile(fileName=filename, format=PNG, canvasObjects=(viewport, ))
    session.odbs[filename+'.odb'].close()
    return
#
filenames = ['example_VUSDFLD_V2','example_USDFLD_V2']
frames = [50,695]
for filename,frame in zip(filenames,frames):
    post_pictures(filename,frame)
