import sys
import numpy

from pywsclean import *

if len(sys.argv)<2:
	print 'Syntax: example.py <ms>\n'
else:
	parameters = ImagingParameters()
	parameters.msPath = sys.argv[1]
	parameters.imageWidth = 1024
	parameters.imageHeight = 1024
	parameters.pixelScaleX = '1amin'
	parameters.pixelScaleY = '1amin'
	parameters.extraParameters = '-weight natural'

# Test the operator
	with Operator(parameters) as o:
	
		data,weights = o.read()
		
		image = numpy.zeros(parameters.imageWidth*parameters.imageHeight)
		
		o.backward(image, data)
		
		o.forward(data, image)
			
		data = numpy.ones(o.data_size(), dtype=numpy.complex128)
		
		o.backward(image, data)
		
		o.write("name.fits", image)
		
		try:
			o.backward(data, image) # data and image are swapped
		except Exception as e:
			print '- Ok, specifying wrong input raised:', e

	try:
		o.backward(image, data) # Outside with block, should raise
	except Exception as e:
		print '- Ok, using operator outside with block raised:', e

# Test if a second use of operator works fine
	with Operator(parameters) as o:
		data,weights = o.read() # Reread data
		try:
			o.backward(image, numpy.ones(o.data_size(), dtype=numpy.complex64)) # Use wrong type
		except Exception as e:
			print '- Ok, wrong type in backward raised:', e
		try:
			o.forward(numpy.ones(o.data_size(), dtype=numpy.complex64), image) # Use wrong type
		except Exception as e:
			print '- Ok, wrong type in forward raised:', e

# Test the full cleaning command
	wsc=WSClean()
	wsc.width=1536
	wsc.height=1536
	wsc.scale='5asec'
	wsc.datacolumn='DATA'
	wsc.niter = 250

	wsc.set_uniform_weighting()
	# Alternatively: wsc.set_natural_weighting()

	# This makes wsclean-dirty.fits (et al) from the DATA column
	wsc.image([sys.argv[1]], 'wsclean')
	
	# This predicts from wsclean-model.fits into the MODEL_DATA column
	# (note that the column is always called MODEL_DATA, the datacolumn parameter
	#  only sets the column used in the 'image' command)
	wsc.predict([sys.argv[1]], 'wsclean')


