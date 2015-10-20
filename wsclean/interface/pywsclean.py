"""
A Python interface to WSClean
This wrapper can be used to call the (C++) WSClean imager.
"""

import _wsclean
import numpy

class ImagingParameters(object):
	"""Parameters for imaging"""
	msPath = ''
	imageWidth = 512
	imageHeight = 512
	pixelScaleX = '1asec'
	pixelScaleY = '1asec'
	extraParameters = ''

class ImagingData(object):
	"""Information about the imaging operation"""
	dataSize = 0
	
class Operator(object):
	"""Class that wraps WSClean as an operator, so that it is easy
	to get an image from data 'in memory' (and the inverse). Currently, the
	operator will write that data to the MODEL_DATA Measurement Set before imaging."""
	_userdata = None;
	_parameters = None;
	_imagingdata = None;
	
	def __init__(self, parameters):
		"""Constructor: initialize WSClean"""
		self._parameters = parameters
		self._userdata,self._imagingdata = _wsclean.initialize(parameters)
		return
	
	def __del__(self):
		"""Destructor: release WSClean resources"""
		if self._userdata != None:
			print 'Releasing resources for WSClean...'
			_wsclean.deinitialize(self._userdata)

	def data_size(self):
		"""Get the number of visibilities"""
		return self._imagingdata.dataSize

	def read(self):
		"""Read the visibilities and return as a (data,weight) tuple. """
		print 'Reading '+str(self.data_size())+' samples...'
		data = numpy.ascontiguousarray(numpy.zeros(self._imagingdata.dataSize, dtype=numpy.complex128))
		weights = numpy.ascontiguousarray(numpy.zeros(self._imagingdata.dataSize, dtype=numpy.float64))
		_wsclean.read(self._userdata, data, weights)
		return data,weights

	def write(self, filename, data):
		"""Write a FITS image with the correct keywords etc."""
		dataCont = numpy.ascontiguousarray(data)
		_wsclean.write(self._userdata, filename, data)

	def forward(self, dataOut, dataIn):
		"""Perform the forward operation. This is 'prediction': convert
		an image into visibilities. dataOut should be an complex double array
		that will be filled with visibilities, dataIn should be an array
		of doubles, representing the image for the operator input."""
		dataOutCont = numpy.ascontiguousarray(dataOut)
		dataInCont = numpy.ascontiguousarray(dataIn)
		_wsclean.operator_A(self._userdata, dataOut, dataIn)
		
	def backward(self, dataOut, dataIn):
		"""Perform the backward operation. This is the 'imaging' step:
		convert visibilities into an image. dataOut should be an array
		of doubles, which will be filled with the image, dataOut should be an array
		of complex doubles, representing the visibilities for the operator input."""
		dataOutCont = numpy.ascontiguousarray(dataOut)
		dataInCont = numpy.ascontiguousarray(dataIn)
		_wsclean.operator_At(self._userdata, dataOut, dataIn)


class WSClean(object):
	"""The Python interface to WSClean
	"""
	
	datacolumn=''
	"""The column used for imaging; empty means CORRECTED_DATA if it exists, otherwise
	use DATA."""
	
	width=1024
	"""Image width"""
	
	height=1024
	"""Image height"""
	
	scale='1asec'
	"""Pixel scale of image. Units can e.g. be deg, amin, asec, masec. There should not
	be a space between the number and its unit."""
	
	niter=0
	"""Number of clean or moresane iterations"""
	
	gain=-1
	"""Gain per minor iteration. -1 means use WSClean's default."""
	
	mgain=-1
	"""Gain per major iteration. -1 means use WSClean's default."""
	
	__weightpar=''
	
	def __init__(self):
		return;
	
	def image(self, msnames, nameprefix):
		"""Run WSClean to make an image on the specified list of measurement sets"""
		plist=self.__get_parameterlist(nameprefix)
		import os
		msnamelist=' '.join(msnames)
		cmd='wsclean '+str(plist)+' '+msnamelist
		print cmd
		os.system(cmd)
		return;
	
	def predict(self, msnames, nameprefix):
		"""Run WSClean to predict"""
		plist=self.__get_parameterlist(nameprefix)
		import os
		msnamelist=' '.join(msnames)
		cmd='wsclean -predict '+str(plist)+' '+msnamelist
		print cmd
		os.system(cmd)
		return;
	
	def set_uniform_weighting(self):
		"""Enable uniform weighting"""
		self.__weightpar='-weight uniform'
		
	def set_natural_weighting(self):
		"""Enable natural weighting"""
		self.__weightpar='-weight natural'
		
	def set_briggs_weighting(self, robustness):
		"""Enable Briggs' weighting with a given robustness"""
		self.__weightpar='-weight briggs '+str(robustness)
	
	def __get_parameterlist(self, prefixname):
		plist = '-size '+str(self.width)+' '+str(self.height)+' -scale '+str(self.scale);
		if self.datacolumn!='':
			plist += ' -datacolumn '+self.datacolumn;
		if self.__weightpar!='':
			plist += ' '+self.__weightpar;
		if self.niter!=0:
			plist += ' -niter '+str(self.niter);
		if self.gain!=-1:
			plist += ' -gain '+str(self.gain);
		if self.mgain!=-1:
			plist += ' -mgain '+str(self.mgain);
		if prefixname!='':
			plist += ' -name '+prefixname;
		return plist;
