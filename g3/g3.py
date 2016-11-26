#!/bin/python
def main(): # main_script
    """
    This routine reads 'drop_fluct.mide', and writes out the time correlation function
    """
    import numpy as np
    from fort_lib import g_3 as fort_g3    
    import sys
   
    #Check if arguments are passed corectly to script
    if len(sys.argv) != 3:
        print "ERROR: This script must be called so:"
        print "python g3.py file_in col_nr"
        print "col_nr is number of data columns. First column should be time"
        exit()
           
    #Set out file
    fl_out = 'g3_vs_t.dat'

    #Read data
    input_file = str(sys.argv[1])
    n_cols = int(sys.argv[2])
    time = read_time(input_file,1) 
    x = read_data(input_file, n_cols)
	     
    #Do heavy calculation
    g3_data = fort_g3(np.transpose(x),x.shape[0],n_cols) # Call fortran routine # As fast as Muhammad Ali
    #autocorr = corr_np(x) # Call nupy calculation    

    #Save data
    save_data(fl_out, g3_data, time)    

def corr_np(x):
    """
    Numpy implementation of autocorrelation.
    Source: http://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation
    """
    import numpy as np
    yunbiased = x-np.mean(x)
    ynorm = np.sum(yunbiased**2)
    acor = np.correlate(yunbiased, yunbiased, "same")/ynorm
    acor = acor[len(acor)/2:]
    return acor

def save_data(file_name,data_in,fst_col=0):
    """Save data (2nd column) to  file passed as string in first argument. 
    OPTIONAL: If a 3rd argument is passed (must be a numpy array of the same length as data), this will be printed as
    first column. If no 3rd argument is passed, the first column will be 0, 1, 2, ..., n """
    import numpy as np
    n = data_in.shape[0]
    g = open(file_name, 'w+')
    
    if type(fst_col) is int:
        fst_col = range(n)

    for i in range(n):
        s = str(str(fst_col[i])+' '+ str(data_in[i]) + '\n') 
        g.write(s)

    g.close()   


def read_data(file_name,col_nr=1):
    """Read file passed as argument and return it as a numpy array"""
    import numpy as np
    import gzip
     
    col_nr = col_nr - 1

    data = []
    
    if file_name[len(file_name)-3:len(file_name)] == '.gz':   
        with gzip.open(file_name,'r') as f:
            for line in f:
                data.append( np.asarray( line.split()[1:col_nr+2] ) )
    else:
        with open(file_name,'r') as f:
            for line in f:
                data.append( np.asarray( line.split()[1:col_nr+2] ) )

    x = np.asarray(data)
    
    return x

def read_time(file_name,col_nr=1):
    """Read file passed as argument and return it as a numpy array"""
    import numpy as np
    import gzip
     
    col_nr = col_nr - 1

    data = []
    
    if file_name[len(file_name)-3:len(file_name)] == '.gz':   
        with gzip.open(file_name,'r') as f:
            for line in f:
                data.append( float( line.split()[col_nr] ) )
    else:
        with open(file_name,'r') as f:
            for line in f:
                data.append( float( line.split()[col_nr] ) )

    x = np.asarray(data)

    return x


########### END SUPPORT ROUTINES ##########

if __name__ == "__main__":
    main()

