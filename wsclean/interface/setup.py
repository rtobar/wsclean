from distutils.core import setup, Extension

# define the extension module
_wsclean = Extension('_wsclean', sources=['_wsclean.c'], libraries=['../build/wsclean-shared','casa_ms','casa_tables','casa_casa','casa_measures','casa_fits','fftw3','fftw3_threads','boost_date_time','boost_system','boost_thread','boost_filesystem','boost_program_options','cfitsio','gsl','stationresponse'])

# run the setup
setup(ext_modules=[_wsclean])
