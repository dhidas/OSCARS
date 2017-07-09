import os.path

def read_file_list(ifile, idir=None):
    """
    read a list of parameters and filenames from a file.
    
    The file format should be two columns separated by whitespace.
    In the first column is the parameter value as a floating point number
    and in the second column the file name.  Whitespace in the name is
    unforgivably wrong, but technically allowed.
    
    Parameters
    ----------
    ifile : str
        Full path to file containing the list of parameters and filenames
        
    idir : str
        Path to directory where the files are contained if not in the 'ifile' directory
        
    Returns
    -------
    file_list : list
        List of [parameter, filename]s
    """
    
    # Directory where files are
    mydir = os.path.dirname(ifile)
    if idir is not None:
        mydir = idir
    
    # List we will return
    mylist=[]
    with open(ifile, 'r') as fi:
        for l in fi:
            ls = l.split()
            pv = float(ls[0])
            fn = os.path.join(mydir, ' '.join(ls[1:]))
            
            mylist.append([pv, fn])
            
    return mylist





def read_file_list2(ifile, gap=None, phase_mode=None, phase=None, idir=None):
    """
    read a list of parameters and filenames from a file.
    
    The file format should be four columns separated by whitespace.
    The columns are gap value, phase mode, phase value, file name.
    Whitespace in the phase mode is not allowed.

    One can specify the following:
        gap + phase_mode
        phase + phase_mode
    
    Parameters
    ----------
    ifile : str
        Full path to file containing the list of parameters and filenames

    gap : float
        Gap value of interest (must be exact match, not interpolated)

    phase_mode : str
        Phase mode of interest

    phase : float
        Phase value
        
    idir : str
        Path to directory where the files are contained if not in the 'ifile' directory
        
    Returns
    -------
    file_list : list
        List of [gap, phase_mode, phase, filename]s
    """
    
    # Directory where files are
    mydir = os.path.dirname(ifile)
    if idir is not None:
        mydir = idir
    
    print(mydir)
    # List we will return
    mylist=[]
    with open(ifile, 'r') as fi:
        for l in fi:
            ls = l.split()
            if len(ls) < 4:
                continue
                
            g  = float(ls[0])
            pm = ls[1]
            ph = float(ls[2])
            fn = os.path.join(mydir, ' '.join(ls[3:]))
     
            if phase_mode is None and phase is None and gap == None:
                mylist.append([g, pm, ph, fn])
            elif gap is not None and phase_mode is not None:
                if g == gap and phase_mode == pm:
                    mylist.append([ph, fn])
            elif phase is not None and phase_mode is not None:
                if phase == ph and phase_mode == pm:
                    mylist.append([g, fn])
            
    return mylist
