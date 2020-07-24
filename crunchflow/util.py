import os
import re

def correct_exponent(filename, folder='.', verbose='med'):
    """CrunchFlow has trouble outputting triple-digit exponents and omits 
    the 'E'. For example, '2.5582E-180' prints as '2.5582-180'. Correct 
    all such instances in the provided file.
    
    params:
        filename [str]: name of the file to be processed
        folder [str]: folder containing the file, either relative or 
            absolute path
        verbose [str]: Print each correction as it's performed ('high'), 
            print total number of corrections ('med'), or print nothing 
            ('low'). Default: 'med'
        
    """
    
    crunch_file = os.path.join(folder, filename)
    n_repl = 0 # Count number of replacements

    # Pre-compile regex to save time
    neg_search = re.compile('([0-9][0-9])\-([0-9][0-9][0-9])')
    pos_search = re.compile('([0-9][0-9])\+([0-9][0-9][0-9])')

    # Read in crunch input file
    with open(crunch_file, 'r') as fin:
        cf = fin.readlines()

    with open(crunch_file, 'w') as fout:
        for line in cf:
            # before any changes, store line for printing 
            tmp = line

            if re.search(neg_search, line):  
                line = re.sub(neg_search, r'\1E-\2', line)
                n_repl += 1

            if re.search(pos_search, line):  
                line = re.sub(pos_search, r'\1E+\2', line)
                n_repl += 1
                    
            if verbose == 'high':
                print("Original: \n\t " + tmp)
                print("New: \n\t " + line)

                
            fout.write(line)
            
    if verbose=='med':
        print("Made {} replacements in {}".format(n_repl, crunch_file))


def crunch_input_block(line, block):
    """While reading a CrunchFlow input file, return the input block to 
    which the current line belongs.
    
    params:
        line [str]: current line in the input file
        block [str]: block to which the previous line belongs. 
            Initialize as empty str
        
    returns:
        block [str]: one of the default CrunchFlow blocks (e.g., 
            "RUNTIME") or an empty string if between blocks. If 
            within a geochemical condition block, returns the name 
            of the geochemical condition.
    
    Example usage:
    
        >>> block = ""
        >>> with open(input_file, "r") as f:
        >>>     for line in f:
        >>>         block = crunch_input_block(line, block)
        >>>         
        >>>         if block == "BOUNDARY_CONDITIONS":
        >>>             print(line)
    
    """
    
    # Strip whitespace
    line = line.strip()
    
    # List of valid CrunchFlow input blocks
    cr_fields = ['TITLE', 'MINERALS', 'AQUEOUS_KINETICS', 'RUNTIME', 
                 'PEST', 'PRIMARY_SPECIES', 'SECONDARY_SPECIES', 
                 'SURFACE_COMPLEXATION', 'GASES', 'TRANSPORT', 
                 'BOUNDARY_CONDITIONS', 'TEMPERATURE', 'ISOTOPES', 
                 'ION_EXCHANGE', 'POROSITY', 'DISCRETIZATION', 
                 'INITIAL_CONDITIONS', 'FLOW', 'OUTPUT', 'EROSION']

    # If it's all uppercase, we are at the beginning or end of a block
    if line.isupper():
        if line in cr_fields:
            block = line
            
        elif line == 'END':
            block = ''
    
    # Geochemical conditions block
    elif 'Condition' in line:
        if line.split()[0] == 'Condition':
            block = line.split()[1]
    
    return block
