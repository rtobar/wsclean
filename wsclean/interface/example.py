import sys
import numpy

from pywsclean import *

if len(sys.argv)<2:
	print 'Syntax: example.py <ms>\n'
else:
	parameters = ImagingParameters()
	parameters.msPath = sys.argv[1]
	parameters.imageWidth = 512
	parameters.imageHeight = 512
	parameters.pixelScaleX = '1amin'
	parameters.pixelScaleY = '1amin'
	parameters.extraParameters = '-weight natural'

# Test the operator
	o = Operator(parameters)
	
	data,weights = o.read()

	image = numpy.zeros(parameters.imageWidth*parameters.imageHeight)
	
	o.backward(image, data)
	
	o.forward(data, image)
		
	data = numpy.ones(o.data_size(), dtype=numpy.complex128)
	
	o.backward(image, data)
	
	o.write("name.fits", image)

# Test the full cleaning command
	wsc=WSClean()
	wsc.width=2048
	wsc.height=2048
	wsc.scale='5asec'
	wsc.datacolumn='DATA'
	wsc.niter = 1000

	wsc.set_uniform_weighting()
	# Alternatively: wsc.set_natural_weighting()

	# This makes wsclean-dirty.fits (et al) from the DATA column
	wsc.image([sys.argv[1]], 'wsclean')
	
	# This predicts from wsclean-model.fits into the MODEL_DATA column
	# (note that the column is always called MODEL_DATA, the datacolumn parameter
	#  only sets the column used in the 'image' command)
	wsc.predict([sys.argv[1]], 'wsclean')

