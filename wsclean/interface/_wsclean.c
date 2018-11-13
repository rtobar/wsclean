#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>

#include <numpy/arrayobject.h>

#include "wscleaninterface.h"

static PyObject* getAttribute(PyObject* o, const char* attr)
{
	PyObject* attrObj = PyObject_GetAttrString(o, attr);
	if(attrObj == NULL)
	{
		PyErr_SetString(PyExc_RuntimeError, "missing attribute");
	}
	return attrObj;
}

static const char* getStringAttribute(PyObject* o, const char* attr)
{
	PyObject* attrObj = getAttribute(o, attr);
	if(attrObj == NULL)
		return NULL;
	else
		return PyString_AsString(attrObj);
}

static int getIntAttribute(PyObject* o, const char* attr, int* isValid)
{
	PyObject* attrObj = getAttribute(o, attr);
	if(attrObj == NULL)
	{
		*isValid = 0;
		return 0;
	}
	else {
		if(PyInt_Check(attrObj))
		{
			*isValid = 1;
			return PyInt_AsLong(attrObj);
		}
		else {
			*isValid = 0;
			PyErr_SetString(PyExc_RuntimeError, "not an int");
			return 0;
		}
	}
}

static void* getUserData(PyObject* capsule)
{
	if(!PyCapsule_CheckExact(capsule))
	{
		PyErr_SetString(PyExc_RuntimeError, "Invalid parameter specified; expecting a userdata structure");
		return NULL;
	}
	return PyCapsule_GetPointer(capsule, "pywsclean.userdata");
}

PyObject* setException(const char* c)
{
	PyErr_SetString(PyExc_RuntimeError, c);
	return NULL;
}

static PyObject* wsclean_initialize_func(PyObject* self, PyObject* args)
{
	PyObject *parameters;
  int ok = PyArg_ParseTuple(args, "O", &parameters);
	if(!ok)
		return NULL;
	
	//PyObject_Print(parameters, stdout, 0); printf("\n");
	//PyObject_Print(data, stdout, 0); printf("\n");
	
	imaging_parameters p;
	p.msPath = getStringAttribute(parameters, "msPath");
	if(!p.msPath) return NULL;
	int valid;
	p.imageWidth = getIntAttribute(parameters, "imageWidth", &valid);
	if(!valid) return NULL;
	p.imageHeight = getIntAttribute(parameters, "imageHeight", &valid);
	if(!valid) return NULL;
	p.doNormalize = getIntAttribute(parameters, "doNormalize", &valid);
	if(!valid) return NULL;
	const char* scaleStr = getStringAttribute(parameters, "pixelScaleX");
	if(!scaleStr) return NULL;
	p.pixelScaleX = wsclean_parse_angle(scaleStr);
	scaleStr = getStringAttribute(parameters, "pixelScaleY");
	if(!scaleStr) return NULL;
	p.pixelScaleY = wsclean_parse_angle(scaleStr);
	p.extraParameters = getStringAttribute(parameters, "extraParameters");
	if(!p.extraParameters) return NULL;

	void *userData;
	imaging_data d;
	wsclean_initialize(&userData, &p, &d);
	
	PyObject* capsule = PyCapsule_New(userData, "pywsclean.userdata", NULL);
	
	PyObject* wscleanModuleString = PyString_FromString((char*)"pywsclean");
	if(!wscleanModuleString)
		return setException("failed creating wsclean string");
	PyObject* wscleanModule = PyImport_Import(wscleanModuleString);
	if(!wscleanModule)
		return setException("failed to import module pywsclean");
	PyObject* imagingDataFunc = PyObject_GetAttrString(wscleanModule, "ImagingData");
	if(!imagingDataFunc)
		return setException("failed to get attribute pywsclean.ImagingData");
	
	PyObject* imagingDataObj = PyType_GenericNew(imagingDataFunc, Py_BuildValue(""), Py_BuildValue(""));
	if(!imagingDataObj)
		return NULL;
	PyObject_SetAttrString(imagingDataObj, "dataSize", PyInt_FromLong(d.dataSize));
	
	// No return value:
	return Py_BuildValue("OO", capsule, imagingDataObj);
}

static PyObject* wsclean_deinitialize_func(PyObject* self, PyObject* args)
{
	PyObject *pyUserData;
  int ok = PyArg_ParseTuple(args, "O", &pyUserData);
	if(!ok)
		return NULL;
	void* userData = getUserData(pyUserData);
	if(!userData)
		return NULL;
	wsclean_deinitialize(userData);
	return Py_BuildValue("");
}

static PyObject* wsclean_read_func(PyObject* self, PyObject* args)
{
	PyObject *pyUserData, *pyData, *pyWeights;
  int ok = PyArg_ParseTuple(args, "OOO", &pyUserData, &pyData, &pyWeights);
	if(!ok)
		return NULL;
	void* userData = getUserData(pyUserData);
	if(!userData)
		return NULL;
	
	DCOMPLEX *data = (DCOMPLEX *) PyArray_DATA(pyData);
	if(!data)
		return NULL;
	double *weights = (double *) PyArray_DATA(pyWeights);
	if(!weights)
		return NULL;
	
	wsclean_read(userData, data, weights);
	
	return Py_BuildValue("");
}

static PyObject* wsclean_write_func(PyObject* self, PyObject* args)
{
	PyObject *pyUserData, *pyFilename, *pyImage;
  int ok = PyArg_ParseTuple(args, "OOO", &pyUserData, &pyFilename, &pyImage);
	if(!ok)
		return NULL;
	
	void* userData = getUserData(pyUserData);
	if(!userData)
		return NULL;
	
	const char* filename = PyString_AsString(pyFilename);
	if(!filename)
		return NULL;
	
	double *image = (double *) PyArray_DATA(pyImage);
	if(!image)
		return NULL;
	
	wsclean_write(userData, filename, image);
	
	return Py_BuildValue("");
}

static PyObject* wsclean_operator_A_func(PyObject* self, PyObject* args)
{
	PyObject *pyUserData, *pyDest, *pySrc;
  int ok = PyArg_ParseTuple(args, "OOO", &pyUserData, &pyDest, &pySrc);
	if(!ok)
		return NULL;
	void* userData = getUserData(pyUserData);
	if(!userData)
		return NULL;
	
	DCOMPLEX *dest = (DCOMPLEX *) PyArray_DATA(pyDest);
	if(!dest)
		return NULL;
	
	double *src = (double *) PyArray_DATA(pySrc);
	if(!src)
		return NULL;
	
	wsclean_operator_A(userData, dest, src);
	
	return Py_BuildValue("");
}


static PyObject* wsclean_operator_At_func(PyObject* self, PyObject* args)
{
	PyObject *pyUserData, *pyDest, *pySrc;
  int ok = PyArg_ParseTuple(args, "OOO", &pyUserData, &pyDest, &pySrc);
	if(!ok)
		return NULL;
	void* userData = getUserData(pyUserData);
	if(!userData)
		return NULL;
	
	double *dest = (double *) PyArray_DATA(pyDest);
	if(!dest)
		return NULL;
	
	DCOMPLEX *src = (DCOMPLEX *) PyArray_DATA(pySrc);
	if(!src)
		return NULL;
	
	wsclean_operator_At(userData, dest, src);
	
	return Py_BuildValue("");
}

static PyMethodDef WSCleanMethods[] = {
 { "initialize", wsclean_initialize_func, METH_VARARGS, "Initialize the WSClean interface" },
 { "deinitialize", wsclean_deinitialize_func, METH_VARARGS, "Deinitialize the WSClean interface" },
 { "read", wsclean_read_func, METH_VARARGS, "Read the initial data (visibilities) from disk" },
 { "write", wsclean_write_func, METH_VARARGS, "Write an image result to disk" },
 { "operator_A", wsclean_operator_A_func, METH_VARARGS, "Apply the forward operator (predict visibilities)" },
 { "operator_At", wsclean_operator_At_func, METH_VARARGS, "Apply the reverse operator (make a dirty image)" },
 { NULL, NULL, 0, NULL }
};

DL_EXPORT(void) init_wsclean(void)
{
  Py_InitModule("_wsclean", WSCleanMethods);
}

