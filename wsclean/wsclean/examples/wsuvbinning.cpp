#include "../binneduvoutput.h"

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		std::cout
			<< "Syntax: wsuvbinningexample <uvcoverage prefix> <dirty imaging prefix>\n"
			<< "Advice for uvcoverage image, is to make it with parameters:\n"
			<< "\t-no-small-inversion -padding 1 -grid-mode nn -nwlayers 1 -make-psf-only\n"
			<< "For the dirty image:\n"
			<< "\t-no-small-inversion -weight natural -make-psf\n";
	}
	else {
		BinnedUVOutput::Make(argv[1], argv[2]);
	}
}

