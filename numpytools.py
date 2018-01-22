"""Bioinformatics toolbox

This is a loose collection of tools that require the numpy package.

Author: Jakob Nybo Nissen, DTU Bioinformatics
"""


import numpy as _np
import collections as _collections
from misctools_c import tetnucfq_safe, tetnucfq_unsafe


def byte_iterfasta(filehandle):
    "Yields writeable Numpy char arrays from a binary opened fasta file."

    # Skip to first header
    for probeline in filehandle:
        if probeline.startswith(b'>'):
            break
    else: # nobreak
        raise ValueError('No headers in this file.')

    header = probeline.strip(b'>\n')
    buffer = list()

    # Iterate over lines
    for line in map(bytes.rstrip, filehandle):
        if line.startswith(b'>'):
            array = _np.frombuffer(b''.join(buffer), dtype=_np.uint8)
            array.flags['WRITEABLE'] = True
            yield header, array
            buffer.clear()
            header = line[1:]

        else:
            buffer.append(line)

    array = _np.frombuffer(b''.join(buffer), dtype=_np.uint8)
    array.flags['WRITEABLE'] = True
    yield header, array


AssemblyStats = _collections.namedtuple('AssemblyStats', ['size', 'n50', 'ncontigs', 'largest', 'smallest', 'sizes'])


def assemblystats(path, xmax=10000, step=100):
    """Expects the output from following command as input:
    For SPAdes: grep '^>' contigs.fasta | cut -d_ -f4 | sort -nr | uniq -c
    For MEGAHIT: grep '^>' final.contigs.fa | cut -d= -f4 | sort -nr | uniq -c
    
    Returns n_contigs, largest_contig, N50, n_nucleotides_above_each_threshold,
    total nucleotides."""
    
    countlengths = list()
    with open(path) as file:
        for line in file:
            fields = line.split()
            countlengths.append((int(fields[0]), int(fields[1])))
    
    lengths = _np.zeros(xmax//step + 1, dtype=_np.int)
    assemblysize = sum(count*length for count, length in countlengths)
    ncontigs = sum(count for count, length in countlengths)
    largestcontig = countlengths[0][1]
    smallestcontig = countlengths[-1][1]
    N50 = None
    
    # Length distribution
    for count, length in countlengths:
        lengths[:(length//step)+1] += length * count
        if lengths[0] >= assemblysize and N50 is None:
            N50 = length
            
    return AssemblyStats(assemblysize, N50, ncontigs, largestcontig, smallestcontig, lengths)




