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





def read_file_list2(ifile, idir=None):
    """
    read a list of parameters and filenames from a file.
    
    The file format should be four columns separated by whitespace.
    The columns are gap value, phase mode, phase value, file name.
    Whitespace in the phase mode is not allowed.  
    
    Parameters
    ----------
    ifile : str
        Full path to file containing the list of parameters and filenames
        
    idir : str
        Path to directory where the files are contained if not in the 'ifile' directory
        
    Returns
    -------
    file_list : list
        List of [gap, phase mode, phase, filename]s
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
                
            gap = float(ls[0])
            pm = ls[1]
            phase = float(ls[2])
            fn = os.path.join(mydir, ' '.join(ls[3:]))
            
            mylist.append([gap, pm, phase, fn])
            
    return mylist
