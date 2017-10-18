"""Parses a canopy clusters file and a corresponding gene catalogue.
Outputs a directory with each of the bins as fasta file and/or a "query" amino
acid fasta file with all the genes present in any cluster.

All sequences are stored in memory, you so might want to give it a couple of GB
to work with.
"""


__author__ = 'Jakob Nybo Nissen, DTU Bioinformatics'


import argparse
import os
import sys
from time import time
from collections import defaultdict


genetic_code = {
    ('A', 'A', 'A'): 'K', ('A', 'A', 'G'): 'K', ('A', 'A', 'T'): 'N', ('A', 'A', 'C'): 'N', 
    ('A', 'G', 'A'): 'R', ('A', 'G', 'G'): 'R', ('A', 'G', 'T'): 'S', ('A', 'G', 'C'): 'S', 
    ('A', 'T', 'A'): 'I', ('A', 'T', 'G'): 'M', ('A', 'T', 'T'): 'I', ('A', 'T', 'C'): 'I', 
    ('A', 'C', 'A'): 'T', ('A', 'C', 'G'): 'T', ('A', 'C', 'T'): 'T', ('A', 'C', 'C'): 'T', 
    ('G', 'A', 'A'): 'E', ('G', 'A', 'G'): 'E', ('G', 'A', 'T'): 'D', ('G', 'A', 'C'): 'D', 
    ('G', 'G', 'A'): 'G', ('G', 'G', 'G'): 'G', ('G', 'G', 'T'): 'G', ('G', 'G', 'C'): 'G', 
    ('G', 'T', 'A'): 'V', ('G', 'T', 'G'): 'V', ('G', 'T', 'T'): 'V', ('G', 'T', 'C'): 'V', 
    ('G', 'C', 'A'): 'A', ('G', 'C', 'G'): 'A', ('G', 'C', 'T'): 'A', ('G', 'C', 'C'): 'A', 
    ('T', 'A', 'A'):  '', ('T', 'A', 'G'):  '', ('T', 'A', 'T'): 'Y', ('T', 'A', 'C'): 'Y', 
    ('T', 'G', 'A'):  '', ('T', 'G', 'G'): 'W', ('T', 'G', 'T'): 'C', ('T', 'G', 'C'): 'C', 
    ('T', 'T', 'A'): 'L', ('T', 'T', 'G'): 'L', ('T', 'T', 'T'): 'F', ('T', 'T', 'C'): 'F', 
    ('T', 'C', 'A'): 'S', ('T', 'C', 'G'): 'S', ('T', 'C', 'T'): 'S', ('T', 'C', 'C'): 'S', 
    ('C', 'A', 'A'): 'Q', ('C', 'A', 'G'): 'Q', ('C', 'A', 'T'): 'H', ('C', 'A', 'C'): 'H', 
    ('C', 'G', 'A'): 'R', ('C', 'G', 'G'): 'R', ('C', 'G', 'T'): 'R', ('C', 'G', 'C'): 'R', 
    ('C', 'T', 'A'): 'L', ('C', 'T', 'G'): 'L', ('C', 'T', 'T'): 'L', ('C', 'T', 'C'): 'L', 
    ('C', 'C', 'A'): 'P', ('C', 'C', 'G'): 'P', ('C', 'C', 'T'): 'P', ('C', 'C', 'C'): 'P', 
    }


def mkdir(name, indent=False):
    """Creates a new directory in a threadsafe way."""
    
    try:
        os.mkdir(name)
    except FileExistsError:
        if os.path.isfile(name):
            raise
        print('\t'*indent + 'Directory {} already exists, skipping creation.'.format(name))
    else:
        print('\t'*indent + 'Creating directory "{}".'.format(name))


def timed(function):
    """Decorator adding a timer to a function,
    and prints the time elapsed in the terminal. Just eye candy."""
    
    def inner(*args, **kwargs):
        begin = time()
        result = function(*args, **kwargs)
        print('\tDone in {:,.2f} seconds'.format(time() - begin))
        return result
    
    return inner


def iterfasta(filehandle):
    "Iterate over a fasta file, yielding (first_word_inheader, seq) tuples"
    
    buffer = list()
    
    header = next(filehandle).strip()
    if not header.startswith('>'):
        raise ValueError('First line is not a header')
    header = header.split()[0][1:]
    
    for line in map(str.strip, filehandle):
        if line.startswith('>'): 
            yield header, ''.join(buffer)
            buffer.clear()
            header = line.split()[0][1:]
            
        else:
            buffer.append(line)

    yield header, ''.join(buffer)


@timed
def translatedict(seqdict, code=genetic_code):
    print('Translating sequence dictionary.')
    
    if args.progress:
        status = '\t{{}}/{} k sequences translated.'.format(len(seqdict)//1000)
        for n, (name, sequence) in enumerate(seqdict.items()):
            if n % 1000 == 0:
                print(status.format(n//1000), end='\r')
            seqdict[name] = ''.join([code[codon] for codon in zip(*[iter(sequence)]*3)])
        
        print(' '*70, end='\r') # clear output line
    
    else:
        for name, sequence in seqdict.items():
            seqdict[name] = ''.join([code[codon] for codon in zip(*[iter(sequence)]*3)])


@timed
def init_dicts(mingenes, clusterspath):
    """Reads a canopy cluster file and initializes two dicts:
    bindict: a binname:[genes] dict
    seqdict: a name:seq dict
    """
    
    print('Parsing clusters.')

    seqdict = dict()
    bindict = defaultdict(list)
    
    with open(clusterspath) as file:
        currentbin, gene = next(file).split()
        currentbinbuffer = [(currentbin, gene)]
        
        for bin, gene in map(str.split, file):
            if bin != currentbin:
                
                if len(currentbinbuffer) >= mingenes:
                    for bufbin, bufgene in currentbinbuffer:
                        seqdict[bufgene] = None
                        bindict[bufbin].append(bufgene)
                
                currentbinbuffer.clear()
                currentbin = bin
            
            currentbinbuffer.append((bin, gene))
        
        for bufbin, bufgene in currentbinbuffer:
            seqdict[bufgene] = None
            bindict[bufbin].append(bufgene)

    return bindict, seqdict


@timed
def fill_seqdict(seqdict, cataloguein, as_aa=False, code=genetic_code):
    "Fills the seqdict from a gene catalogue."
    
    progress = args.progress
    print('Filling in genes from gene catalogue.')
    
    # Parsing gene catalogue
    with open(cataloguein) as inputfile:
        for n, (name, sequence) in enumerate(iterfasta(inputfile)):
            if progress and n % 1000 == 0:
                print('\tProcessed {}k genes.'.format(n//1000), end='\r')
                
            if name in seqdict:
                if as_aa:
                    seqdict[name] = ''.join([code[codon] for codon in zip(*[iter(sequence)]*3)])
                else:
                    seqdict[name] = sequence
                
    if progress:
        print(' '*50, end='\r') # clear the line
    
    return None


@timed
def write_bins(bindict, seqdict, bindir):
    "Given a bindict and a filled seqdict, writes nucleotide gene bins"
    
    progress = args.progress
    print('Creating bins.')
    mkdir(bindir, indent=True)
        
    # Check which files exist:
    presentfiles = [i.name for i in os.scandir(bindir)]
    missingfiles = [name+'.fna' for name in bindict if name+'.fna' not in presentfiles]
    
    if progress:
        status = '\tCreated {{}}/{} MGSs'.format(len(missingfiles))
    
    if not missingfiles:
        print('\tAll files already exists. Skipping creation.')
        return
    elif len(missingfiles) != len(bindict):
        nexist = len(bindict) - len(missingfiles)
        print('\t{} files already exist. Creating rest.'.format(nexist))
    
    buffer = list()
    
    for processed, filename in enumerate(missingfiles):
        for gene in bindict[filename[:-4]]:
            buffer.append('>{}\n{}'.format(gene, seqdict[gene]))
        
        fullname = os.path.join(bindir, filename)
        with open(fullname, 'w') as file:
            print('\n'.join(buffer), file=file)
        
        buffer.clear()
        if progress:
            print(status.format(processed), end='\r')
    
    if progress:
        print(' '*50, end='\r') # clear line
    return None


@timed
def write_query(seqdict, queryout):
    """Given a filled seqdict, writes the query."""
    
    print('Creating amino acid query file.')
    
    buffer = list()
    
    with open(queryout, 'w') as query:
        for gene, sequence in seqdict.items():
            buffer.append('>{}\n{}'.format(gene, sequence))
            
            if len(buffer) == 10000:
                print('\n'.join(buffer), file=query)
                buffer.clear()
        
        print('\n'.join(buffer), file=query)
        
    return None


def main(args, code=genetic_code):
    """Execute the entire workflow of this script."""
    
    begin_time = time()
    
    # In all cases, init the dicts
    bindict, seqdict = init_dicts(args.mingenes, args.clusters)
    
    # If bins are desired, construct directory, then load as NT, then convert
    if args.bindir is not None:
        fill_seqdict(seqdict, args.cataloguein, as_aa=False)
        write_bins(bindict, seqdict, args.bindir)
    
    # If query is desired
    if args.queryout is not None:
        # If bins have been created, seqdict has been created, must be translated
        if args.bindir is not None:
            translatedict(seqdict)
        
        # Else, fill it directly with amino acids
        else:
            fill_seqdict(seqdict, args.cataloguein, as_aa=True)

        write_query(seqdict, args.queryout)
        
    print('Complete. Total elapsed time:', round(time() - begin_time, 2), 'seconds.')
    return None


# Arg parsing, checking necessary files etc
if __name__ == '__main__':
    usage = 'python parsecanopy.py clusters catalogue [-q query] [-b bindir] [options]'
    
    # Parse input
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=usage)
    
    # Required, positional arguments
    parser.add_argument('clusters', help='Canopy clusters path',
                        metavar='clusters')
    parser.add_argument('cataloguein', help='path to gene catalogue',
                        metavar='catalogue')
    
    outputgroup = parser.add_argument_group('output')
    outputgroup.add_argument('-q', dest='queryout', help='path to write query to',
                            metavar='query')
    outputgroup.add_argument('-b', dest='bindir', help='directory to write bins to',
                            metavar='bindir')
    
    parser.add_argument('-m', dest='mingenes', type=int, default=1, metavar='mingenes',
                        help='minimum size of genes per cluster [1]')
    parser.add_argument('--progress', action='store_true',
        help='Print progress continually to stdout [False]')
    
    # Invoke help if called with no argument
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    
    # Must state either make bins or query
    if not (args.queryout or args.bindir):
        raise ValueError('Must make either bins or query file.')
    
    # Input files must actually exist
    for path in (args.clusters, args.cataloguein):
        if not os.path.exists(path):
            raise FileNotFoundError(path)
            
    if args.queryout is not None and os.path.exists(args.queryout):
        raise FileExistsError(args.queryout)
        
    if args.bindir is not None and os.path.isfile(args.bindir):
        raise FileExistsError('{} is a regular existing file.'.format(args.queryout))
        
    main(args)

