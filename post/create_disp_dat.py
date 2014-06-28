#!/bin/env python
"""
create_disp_dat.py

Create disp.dat file from nodout file.

This is replacing StuctPost, which relied on LS-PREPOST, to extract data from
d3plot* files.  (LS-PREPOST no longer works gracefully on the cluster w/o
GTK/video support.)  Instead of working with d3plot files, this approach now
utilizes ASCII nodout files.  Also replaced the Matlab scritps, so this should
run self-contained w/ less dependencies.

EXAMPLE
=======
create_disp_dat.py

=======
Copyright 2014 Mark L. Palmeri (mlp6@duke.edu)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

__author__ = "Mark Palmeri"
__email__ = "mlp6@duke.edu"
__license__ = "Apache v2.0"


def main():
    import sys

    # lets read in some command-line arguments
    args = parse_cli()

    # open dispout for binary writing
    dispout = open(args.dispout, 'wb')

    generate_header(args.nodedyn, args.dynadeck)

    write_header(header, dispout)

    # open nodout file
    # REFACTOR: try/except or with
    if args.nodout.endswith('gz'):
        import gzip
        print("Extracting gzip-compressed data . . .\n")
        nodout = gzip.open(args.nodout, 'r')
    else:
        print("Extracting data . . .\n")
        nodout = open(args.nodout, 'r')

    timestep_count = 0
    written_count = 0
    
    write_to_cli('Processing Time Step: ')    
    
    for line in nodout:
        if 'nodal' in line:
            timestep_read = True
            write_to_cli('%i ' % timestep_count)
            timestep_count = timestep_count + 1
            data = []
            node_count = 0
            continue
        if timestep_read is True:
            raw_data = line.split()
            corrected_raw_data = correct_Enot(raw_data)
            data.append(map(float, corrected_raw_data))
            node_count = node_count + 1
            if node_count == header['numnodes']:
                timestep_read = False
                process_timestep_data(data, dispout)
                written_count = written_count + 1

    assert (written_count == header['numtimesteps']), 'Mismatch in number of timesteps'

    dispout.close()
    nodout.close()


def parse_cli():
    """
    parse command-line interface arguments
    """
    import sys
    
    if sys.version_info[:2] < (2, 7):
        sys.exit("ERROR: Requires Python >= 2.7")
    
    import argparse

    parser = argparse.ArgumentParser(description="Generate disp.dat "
                                     "data from an ls-dyna nodout file.",
                                     formatter_class=
                                     argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--nodout",
                        help="ASCII file containing nodout data",
                        default="nodout.gz")
    parser.add_argument("--dispout", help="name of the binary displacement "
                        "output file", default="disp.dat")
    parser.add_argument("--nodedyn", help="nodes.dyn input file", 
                        default="nodes.dyn")    
    parser.add_argument("--dynadeck", help="dynadeck.dyn master input deck",
                        default="dynadeck.dyn")
                        
    args = parser.parse_args()

    return args


def generate_header(nodedyn, dynadeck):
    """
    generate headers & write to disp.dat
    
    INPUTS: dispout ('disp.dat') [binary filename to write to]
            nodedyn ('nodes.dyn') [used to determine # nodes & timesteps]
            
    OUTPUTS: header w/ # nodes, dims and timesteps written to disp.dat
             header (dict: numnodes, numdims, numtimesteps)
    """
    header = {}
    header['numnodes'] = count_nodes(nodedyn)
    header['numdims'] = 4  # node ID, x-val, y-val, z-val
    header['numtimesteps'] = calc_timesteps(dynadeck)

    return header


def calc_timesteps(dynadeck):
    """
    calculate # of expected time steps in nodout based on *CONTROL_TERMINATION
    and *DATABASE_NODOUT entries (first non-comment entry on line after
    keyword) in dynadeck.dyn
    
    INPUT: dynadeck ('dynadeck.dyn')
    
    OUTPUT: numtimesteps (int)
    """
    term_time = float(read_keyword_value('\*CONTROL_TERMINATION', 0, dynadeck))
    dt = float(read_keyword_value('\*DATABASE_NODOUT', 0, dynadeck))
    
    numtimesteps = int(term_time/dt)
    return numtimesteps


def read_keyword_value(keyword, column, dynadeck):
    """
    read first entry for a keyword, skipping any comment lines if present
    
    INPUTS: keyword ('\*CONTROL_TERMINATION')
                !! make sure that * in keyword is escaped as \* !!
            column (int) [comma-delimited column to read for value, 0 = 1st col]
            dynadeck ('dynadeck.dyn')
            
    OUTPUT: keyword_value
    """
    import re
    r = re.compile(keyword)    
    readNextNonCommentLine = False
    with open(dynadeck, 'r') as d:
        for line in d:
            if readNextNonCommentLine is True:
                # skip next line if it is a comment
                if line[0] is '$':
                    continue
                else:
                    keyword_value = float(line.split(',')[column])
                    return keyword_value
            if r.match(line):
                readNextNonCommentLine = True
    

def count_nodes(nodefile):
    """
    count # nodes from nodes.dyn
    """
    import fem_mesh
    import numpy as n
    header_comment_skips = fem_mesh.count_header_comment_skips(nodefile)
    nodeIDcoords = n.loadtxt(nodefile,
                             delimiter=',',
                             skiprows=header_comment_skips,
                             comments='*',
                             dtype=[('id', 'i4'), ('x', 'f4'),
                                    ('y', 'f4'), ('z', 'f4')])
    numNodes = len(nodeIDcoords)
    return numNodes
    
    
def write_header(header, outfile):
    '''
    write binary header information to reformat things on read downstream
    'header' is a dictionary containing the necessary information
    '''
    import struct
    outfile.write(struct.pack('fff', header['numnodes'],
                                     header['numdims'], 
                                     header['numtimesteps']))


def process_timestep_data(data, outfile):
    '''
    operate on each time step data row
    '''
    import struct
    # write all node IDs, then x-val, then y-val, then z-val
    [outfile.write(struct.pack('f', data[j][i]))
        for i in [0, 1, 2, 3]
        for j in range(len(data))]

def correct_Enot(raw_data):
    '''
    ls-dyna seems to drop the 'E' when the negative exponent is three digits,
    so check for those in the line data and change those to 'E-100' so that
    we can convert to floats
    '''
    import re
    for i in range(len(raw_data)):
        raw_data[i] = re.sub(r'(?<!E)\-[1-9][0-9][0-9]', 'E-100', raw_data[i])
    return raw_data

def write_to_cli(s):
    """
    write string (s) to CLI
    """
    import sys
    sys.stdout.write(s)
    sys.stdout.flush()

if __name__ == "__main__":
    main()
