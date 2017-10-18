"""Parses a directory of bins, either of contigs or DNA genes. If former, creates
a directory of gene bins. In either case, then parses the gene bins, and
optionally writes two files:
1) a "query" amino acid fasta file with all the genes present in any cluster.
2) A "bin table" - a table of the content of each bin for quick parsing later.

Contigs with same first word in the header are considered identical. Beware.
"""


__author__ = 'Jakob Nybo Nissen, DTU Bioinformatics'


import os
import sys
import argparse
import subprocess
import multiprocessing
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
    "Iterate over a fasta file, yielding (first word in header, seq) tuples"
    
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


def callshell(command):
    "Calls the shell. Must be a top-level function to be parallelizable."
    
    return subprocess.run(command, shell=True, check=True)


@timed
def makebins(contigdir, bindir, extension):
    "Runs prodigal on each of the bins in parallel."
    
    mkdir(args.bindir)
    mkdir(args.logdir)
    
    # Get all files to process:
    infiles = [i.name for i in os.scandir(contigdir) if i.is_file()]
    filenames = [f.rpartition('.')[0] for f in infiles]
    
    infiles = [os.path.join(contigdir, f) for f in infiles if f.endswith(extension)]
    outfiles = [os.path.join(bindir, n + '.fna') for n in filenames]
    
    exist = [os.path.exists(f) for f in outfiles]
    infiles = [f for f, exists in zip(infiles, exist) if not exists]
    outfiles = [f for f, exists in zip(outfiles, exist) if not exists]
    filenames = [f for f, exists in zip(filenames, exist) if not exists]

    if not outfiles:
        print('All aa files already exists, skipping creation...')
        return
    elif len(outfiles) != len(exist):
        print('{}/{} aa files already exists, creating rest.'.format(
              len(exist) - len(outfiles), len(exist)))
    else:
        print('Creating aa files with prodigal ({} cores).'.format(args.cores))
    
    # Create a list of the commands to be executed in shell
    template = '{exec} -p meta -i {infile} -d {outpath} > {stdout} 2> {stderr}'
    commands = list()
    
    for filename, infile, outfile in zip(filenames, infiles, outfiles):
        stdout = '{}/prodigal.{}.stdout'.format(args.logdir, filename)
        stderr = '{}/prodigal.{}.stderr'.format(args.logdir, filename)
        commands.append(template.format(exec=args.prodigalpath,
                                        infile=infile, outpath=outfile,
                                        stdout=stdout, stderr=stderr))
    
    # Execute prodigal in parallel
    if args.progress:
        processes_done = 0
        
        def callback(result, totalps=len(outfiles)):
            "Generator yielding processed"
            nonlocal processes_done
            processes_done += 1
            end = '\n' if processes_done == totalps else '\r'
            print('\tBins processed: {}/{}'.format(processes_done, totalps), end=end)
            return None
    else:
        def callback(result):
            pass

    results = list()
    with multiprocessing.Pool(processes=args.cores) as pool:
        for command in commands:
            results.append(pool.apply_async(callshell, (command,),
                                            callback=callback, error_callback=callback))
        
        pool.close()
        pool.join()
    
    for result in results:
        if not result.successful():
            print('One or more prodigal processed failed.')    
            result.get() # raise its error


@timed
def makequery(bindir, queryout, extension, code=genetic_code):
    "Creates the query file from the gene bins and returns a dict of bins."
    
    print('Parsing bins and creating query file.')
    progress = args.progress
    
    bindict = defaultdict(list)
    seengenes = set()
    buffer = list()
    
    # Get the files with the right extension
    infiles = [os.path.join(bindir, i.name) for i in os.scandir(bindir)
               if i.is_file() and i.name.endswith('.' + extension)]
    print('\tFound {} files with extension {}.'.format(len(infiles), extension))
    
    with open(queryout, 'w') as queryfile:
        # Iterate over the input file, construct bindir and queryfile
        for n, infile in enumerate(infiles):
            basename = os.path.basename(infile).rpartition('.')[0]

            with open(infile) as file:
                for name, sequence in iterfasta(file):
                    bindict[basename].append(name)

                    if name not in seengenes:
                        seengenes.add(name)
                        translated = ''.join([code[codon] for codon in zip(*[iter(sequence)]*3)])
                        buffer.append('>{}\n{}'.format(name, translated))

                        if len(buffer) == 10000:
                            print('\n'.join(buffer), file=queryfile)
                            buffer.clear()

            if progress:
                print('\t{}/{} bins processed.'.format(n, len(infiles)), end='\r')

        print('\n'.join(buffer), file=queryfile)
        
    if progress:
        print(' '*70, end='\r')
    
    return bindict


@timed
def makebintable(bindict, bintableout):
    "Create the bin table from the dictionary returned by makequery."
    
    print('Creating clusters file.')
    
    buffer = list()
    
    for bin, genes in bindict.items():
        separator = '\n' + bin + '\t'
        buffer.append(bin + '\t')
        buffer.append(separator.join(genes))
        buffer.append('\n')
    
    with open(bintableout, 'w') as clustersfile:
        print(''.join(buffer), file=clustersfile, end='')


def main(args):
    "Executes the entire workflow of this script."
    
    begin_time = time()
    # First call prodigal
    if args.contigdir is not None:
        makebins(args.contigdir, args.bindir, args.extension)
    
    # If query is chosen, create that. Extension is .fna when bins are created
    # by makebins, else user-set extension.
    if args.queryout is not None:
        if args.contigdir is not None:
            extension = 'fna'
        else:
            extension = args.extension
            
        bindict = makequery(args.bindir, args.queryout, extension)
        makebintable(bindict, args.bintableout)
    print('Complete. Total elapsed time:', round(time() - begin_time, 2), 'seconds.')


# Command line parsing
if __name__ == '__main__':
    usage = "python parsebins.py bindir [contigdir] [-q query -b bintable] [options]"
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=usage)
    
    cpus = os.cpu_count()
    
    # Positional arguments
    parser.add_argument('bindir', help='directory to find/write gene bins')
    parser.add_argument('contigdir', nargs='?', help='path to bins of contigs (optional)')
    
    # For creation of bin table and query
    bintablegroup = parser.add_argument_group('Query and bintable')
    bintablegroup.add_argument('-q', dest='queryout', metavar='query',
                               help='path to write query to')
    bintablegroup.add_argument('-b', dest='bintableout', metavar='bintable',
                               help='path to write list of bins')
    
    # Optional arguments
    parser.add_argument('-e', dest='extension', default="fna", metavar='extension',
                        help='extension of input bins [fna]')
    parser.add_argument('-p', dest='prodigalpath', metavar='prodigal',
                        default='/services/tools/prodigal/2.6.3/prodigal',
                        help='path to prodigal executable [preset]')
    parser.add_argument('-c', dest='cores', type=int, default=cpus, metavar='cores',
        help='no. of parallel prodigal jobs [{}]'.format(cpus))
    parser.add_argument('--progress', action='store_true',
        help='print progress continually to stdout [False]')

    # If no arguments, print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    
    # Check number of cores
    if args.cores < 1:
        raise argsparse.ArgumentTypeError('Zero or negative cores provided. Exiting')
        
    # You must state either contigdir or queryout
    if args.contigdir is None and args.queryout is None:
        raise ValueError('Input is gene bins, no output specified.')
        
    # If you state query or bintableout, you must also state the other
    if bool(args.queryout) ^ bool(args.bintableout):
        raise ValueError('query and bintable must be stated together.') 
        
    # Prodigal executable must exist in filesystem if called upon
    if args.contigdir is not None and not os.path.isfile(args.prodigalpath):
        raise FileNotFoundError(args.prodigalpath)
    
    # Input directory must refer to an actual directory
    if args.contigdir is None:
        if not os.path.isdir(args.bindir):
            raise FileNotFoundError('bindir: {}'.format(args.bindir))
    else:
        if not os.path.isdir(args.contigdir):
            raise FileNotFoundError('contigdir: {}'.format(args.contigdir))
    
    # Output files must not exist
    if args.queryout is not None and os.path.exists(args.queryout):
        raise FileExistsError(args.queryout)
        
    if args.bintableout is not None and os.path.exists(args.bintableout):
        raise FileExistsError(args.bintableout)
        
    # Make sure the log directory is not a file
    args.logdir = os.path.join(args.bindir, 'prodigallog')
    if args.contigdir is not None and os.path.isfile(args.logdir):
        raise FileExistsError('File exists and is a regular file: {}'.format(args.logdir))

    main(args)

