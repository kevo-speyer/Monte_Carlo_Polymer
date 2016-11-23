F2PY INSCTRUCTIONS:
==== =============
1) Create dynamic library from a fortran file:

    $ f2py -c fort_lib.f90 --fcompiler=intelem -m fort_lib

or if intel fortran compiler is not installed or does not work:    

    $ f2py -c fort_lib.f90 -m fort_lib
  
This will create a new file fort_lib.so wich is callable in any python script

2) Inside the python script import routines from the library:

    >>> from test_f2py import f_routine as corr_f
    >>> y = corr_f(x)
