#include "wscleaninterface.h"

#include "math.h"

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		printf("Specify measurement set: interfaceexample <ms>\n");
	}
	else {
		imaging_parameters dinfo;
		imaging_data format;

		void* userdata;
		dinfo.msPath=argv[1];
		dinfo.imageWidth = 512;
		dinfo.imageHeight = 512;
		dinfo.pixelScaleX = 2.0 * M_PI/(180.0*60.0); // 2 amin
		dinfo.pixelScaleY = 2.0 * M_PI/(180.0*60.0);
		dinfo.extraParameters="-weight natural";
		
		wsclean_initialize(&userdata, &dinfo, &format);
		
		complex double* mydata = (complex double*) malloc(format.dataSize * sizeof(complex double));
		complex double* emptydata = (complex double*) malloc(format.dataSize * sizeof(complex double));
		double* myweights = (double*) malloc(format.dataSize * sizeof(double));
		double* myimage = (double*) malloc(dinfo.imageWidth*dinfo.imageHeight * sizeof(double));
		
		wsclean_operator_At(userdata, myimage, emptydata);
		wsclean_operator_A(userdata, mydata, myimage);
		
		wsclean_read(userdata, mydata, myweights);
		wsclean_operator_At(userdata, myimage, mydata);
		wsclean_operator_A(userdata, mydata, myimage);
		
		size_t i;
		for(i=0; i!=format.dataSize; ++i)
			mydata[i] = 1.0;
		
		wsclean_operator_At(userdata, myimage, mydata);
		wsclean_operator_A(userdata, mydata, myimage);
		
		wsclean_write(userdata, "wsclean-interface-test.fits", myimage);
		
		free(myimage);
		free(mydata);
		
		wsclean_deinitialize(userdata);
	}
}
