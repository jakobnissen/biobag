
"""Bioinformatics toolbox

This is a loose collection of tools for tasks that the author have had to do
too many times, e.g. parsing and translating a Fasta file.

Author: Jakob Nybo Nissen, DTU Bioinformatics
"""



import collections as _collections
import gzip as _gzip
from misctools_c import reverse_complement_kmer
from misctools_c import kmercounts as _kmercounts
from misctools_c import threemerfreq as _threemerfreq
from misctools_c import fourmerfreq as _fourmerfreq
from misctools_c import freq_432mers as _freq_432mers



class Reader:
    "Use this instead of `open` for files which may be gzipped or not."
    
    def __init__(self, filename, readmode='r'):
        if readmode not in ('r', 'rb'):
            raise ValueError("the reader cannot write, set mode to 'r' or 'rb'")
        
        self.filename = filename
        self.readmode = readmode
    
    def __enter__(self):
        with open(self.filename, 'rb') as f:
            signature = f.peek(2)[:2]
        
        # Gzipped files begin with the two bytes 1F8B
        if tuple(signature) == (31, 139):
            if self.readmode == 'r':
                self.filehandle = _gzip.open(self.filename, 'rt')
                
            else:
                self.filehandle = _gzip.open(self.filename, self.readmode)
                
        else:
            self.filehandle = open(self.filename, self.readmode)
            
        return self.filehandle
    
    def __exit__(self, type, value, traceback):
        self.filehandle.close()



def streamprint(iterator, filehandle, bufferlength=10000, sep='\n'):
    """Given an iterator of lines and a filehandle, prints the content of the
    iterator to the file in a memory-efficient way.
    
    Return the number of iterator-given strings processed."""
    
    buffer = list()
    operations = 0
    
    for thing in iterator:      
        buffer.append(str(thing))
        
        if len(buffer) == bufferlength:
            print(sep.join(buffer), file=filehandle, end=sep)
            operations += bufferlength
            buffer.clear()
    
    print(sep.join(buffer), file=filehandle, end=sep)
    return operations + len(buffer)



def forkprint(iterator, *filenames, bufferlength=100, sep='\n'):
    """Given an iterator returning (index, line) and N filenames, prints
    the line to the index'th filename memory-efficiently.
    
    Returns a tuple of the number of printing operations.
    """
    
    operations = [0]*len(filenames)
    
    buffers = list()
    for fh in filenames:
        buffers.append(list())
    
    for index, thing in iterator:
        buffers[index].append(str(thing))
        
        if len(buffers[index]) == bufferlength:
            with open(filenames[index], 'a') as file:
                print(sep.join(buffers[index]), file=file)
                
            operations[index] += bufferlength
            buffers[index].clear()
            
    for index, (buffer, filename) in enumerate(zip(buffers, filenames)):
        if not buffer:
            continue
            
        with open(filename, 'a') as file:
            print(sep.join(buffer), file=file)
            
        operations[index] += len(buffer)
        buffer.clear()
        
    return tuple(operations)



def significant(n, digits=3):
    """"Returns a number rounded to some (default 3) significant digits.
    
    significant(5882) == '5880'
    significant(0.40032, 4) == '0.400'"""
    
    if digits < 1:
        raise ValueError('Digits must be >= 1')
    
    integer_digits = len(str(int(n))) - (n < 0)
    
    if integer_digits >= digits:
        return str(int(round(n, digits - integer_digits)))
    
    else:
        return ('{' + '0:.{}f'.format(digits - integer_digits) + '}').format(n)



class FastaEntry:
    """One single FASTA entry. The header is immutable, and FastaEntries are
    grouped in dicts and sets by their header. Entries are compared equal by
    their sequence
    
    >>> a, b, c = FastaEntry('>one', 'TAG'), FastaEntry('one', 'TAG'), FastaEntry('two', 'CTA')
    >>> a.header == c.header, a.header == b.header # ">" in header removed on instantiation
    False, True
    >>> a == b == c.reversecomplement() # Entries compared by sequence
    True
    >>> a is b, set((a, b)) # distinct objects, hashes same
    False, {<Fasta Entry one>}
    >>> a.translate().sequence, b.translate(endatstop=True)
    '*', None
    """
    
    __slots__ = ['header', 'sequence']
    
    genetic_code = {
    (65, 65, 65): 75, (65, 65, 71): 75, (65, 65, 84): 78, (65, 65, 67): 78,
    (65, 71, 65): 82, (65, 71, 71): 82, (65, 71, 84): 83, (65, 71, 67): 83,
    (65, 84, 65): 73, (65, 84, 71): 77, (65, 84, 84): 73, (65, 84, 67): 73,
    (65, 67, 65): 84, (65, 67, 71): 84, (65, 67, 84): 84, (65, 67, 67): 84,
    (71, 65, 65): 69, (71, 65, 71): 69, (71, 65, 84): 68, (71, 65, 67): 68,
    (71, 71, 65): 71, (71, 71, 71): 71, (71, 71, 84): 71, (71, 71, 67): 71,
    (71, 84, 65): 86, (71, 84, 71): 86, (71, 84, 84): 86, (71, 84, 67): 86,
    (71, 67, 65): 65, (71, 67, 71): 65, (71, 67, 84): 65, (71, 67, 67): 65,
    (84, 65, 65): 42, (84, 65, 71): 42, (84, 65, 84): 89, (84, 65, 67): 89,
    (84, 71, 65): 42, (84, 71, 71): 87, (84, 71, 84): 67, (84, 71, 67): 67,
    (84, 84, 65): 76, (84, 84, 71): 76, (84, 84, 84): 70, (84, 84, 67): 70,
    (84, 67, 65): 83, (84, 67, 71): 83, (84, 67, 84): 83, (84, 67, 67): 83,
    (67, 65, 65): 81, (67, 65, 71): 81, (67, 65, 84): 72, (67, 65, 67): 72,
    (67, 71, 65): 82, (67, 71, 71): 82, (67, 71, 84): 82, (67, 71, 67): 82,
    (67, 84, 65): 76, (67, 84, 71): 76, (67, 84, 84): 76, (67, 84, 67): 76,
    (67, 67, 65): 80, (67, 67, 71): 80, (67, 67, 84): 80, (67, 67, 67): 80,
    }

    
    dna_alphabet = b'ACGT'
    iupacdna_alphabet = b'ACGTMRWSYKVHDBN'
    aa_alphabet = b'ACDEFGHIKLMNPQRSTVWY'
    rna_alphabet = b'ACGU'
    complementtable = bytes.maketrans(b'ACGTMRWSYKVHDBN', b'TGCAKYWSRMBDHVN')
    
    def __init__(self, header, sequence):
        if header[0] in ('>', '#') or header[0].isspace():
            raise ValueError('Header cannot begin with #, > or whitespace')
        self.header = header
            
        if isinstance(sequence, bytearray):
            self.sequence = sequence
        
        elif isinstance(sequence, str):
            self.sequence = bytearray(sequence.encode())
            
        elif isinstance(sequence, bytes):
            self.sequence = bytearray(sequence)
                
        else:
            raise ValueError('sequence must be str, bytes or bytearray')
        
    def __len__(self):
        return len(self.sequence)
    
    def __str__(self):
        return '>{}\n{}'.format(self.header, self.sequence.decode())
    
    def format(self, width=60):
        sixtymers = range(0, len(self.sequence), width)
        spacedseq = '\n'.join([self.sequence[i: i+width].decode() for i in sixtymers])
        return '>{}\n{}'.format(self.header, spacedseq)
    
    # Two entries with same header cannot co-exist in same set/dict!
    def __hash__(self):
        return hash(self.header)
    
    def __contains__(self, other):
        if isinstance(other, str):
            return other.encode() in self.sequence
        
        elif isinstance(other, bytes) or isinstance(other, bytearray):
            return other in self.sequence
        
        else:
            raise TypeError('Can only compare to str, bytes or bytearray')
    
    # Entries are compared equal by their sequence.
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.sequence == other.sequence
        else:
            raise TypeError('Cannot compare to object of other class')
        
    def __getitem__(self, index):
        return self.sequence[index]
        
    def __repr__(self):
        return '<FastaEntry {}>'.format(self.header)
    
    def reversecomplemented(self):
        stripped = self.sequence.translate(None, delete=self.iupacdna_alphabet)
        if len(stripped) > 0:
            raise ValueError("Non-IUPAC DNA char found: '" + stripped[0] +"'")
        
        complemented = self.sequence[::-1].translate(self.complementtable)
        
        return FastaEntry(self.header, complemented)
    
    def check(self, alphabet):
        """This is not done at instantiation because it takes time."""
        
        if alphabet not in (FastaEntry.dna_alphabet, FastaEntry.rna_alphabet,
                            FastaEntry.aa_alphabet, FastaEntry.iupacdna_alphabet):
            
            raise ValueError(('Only accepts dna_alphabet, iupacdna_alphabet,'
                              'rna_alphabet or aa_alphabet of the FastaEntry class'))
            
        # Check if any characters survives removal of entire alphabet
        stripped = self.sequence.translate(None, delete=alphabet)
        if len(stripped) > 0:
            raise ValueError("Invalid character found: '" + chr(stripped[0]) + "'")
    
    def translated(self, endatstop=True):
        codons = zip(*[iter(self.sequence)] * 3)
        try:
            translated_bytes = [self.genetic_code.get(codon) for codon in codons]
            translated = bytearray(translated_bytes)
        except KeyError as exception:
            exception.args = (f'{exception.args[0]} is not a valid DNA codon',)
            raise
        
        if endatstop:
            stoppos = translated.find(42)
            if stoppos == 0:
                return None
            
            elif stoppos != -1:
                translated = translated[:stoppos]
                
        return FastaEntry(self.header, translated)
    
    def kmercounts(self, k):
        if k < 1 or k > 10:
            raise ValueError('k must be between 1 and 10 inclusive')
        return _kmercounts(self.sequence, k)
    
    def fourmer_freq(self):
        return _fourmerfreq(self.sequence)
    
    def threemer_freq(self):
        return _threemerfreq(self.sequence)
    
    def freq_432mers(self):
        return _freq_432mers(self.sequence)



def iterfasta(filehandle, alphabet=None, comment='#'):
    """A generator which yields FastaEntries from an open fasta file.
    
    Usage:
    >>> with open('myfile.fasta') as fastafile:
    ...     entries = iterfasta(fastafile)
    ...
    ...     for entry in entries:
    ...         [ DO STUFF ]
    """
    
    if alphabet is not None:
        if alphabet not in (FastaEntry.dna_alphabet, FastaEntry.rna_alphabet,
                            FastaEntry.aa_alphabet, FastaEntry.iupacdna_alphabet):
            
            raise ValueError(('Only accepts dna_alphabet, iupacdna_alphabet,'
                              'rna_alphabet or aa_alphabet of the FastaEntry class'))
    
    # Skip to first header
    for probeline in filehandle:
        stripped = probeline.lstrip()
        if stripped.startswith(comment):
            pass

        elif probeline[0] == '>':
            break

        else:
            raise TypeError('First non-comment line is not a Fasta header')

    else: # nobreak
        raise TypeError('Empty or outcommented file')
    
    header = probeline[1:-1]
    buffer = list()
    
    # Iterate over lines
    for line in map(str.rstrip, filehandle):
        
        # If line is header, yield the last sequence
        if line[0] == '>':
            yield FastaEntry(header, b''.join(buffer))
            buffer.clear()
            header = line[1:]
        
        # Else check the line and add it to current sequence
        else:
            byteline = line.encode()
            
            if alphabet is not None:
                stripped = byteline.translate(None, delete=alphabet)
                if len(stripped) > 0:
                    raise ValueError("Invalid character found: '" + chr(stripped[0]) +"'")

            buffer.append(byteline)
            
    yield FastaEntry(header, b''.join(buffer))



def simple_iterfasta(filehandle):
    """Yields (header, sequence) tuples from an open fasta file.
    
    Usage:
    >>> with open('myfile.fasta') as fastafile:
    ...     entries = iterfasta(fastafile)
    ...
    ...     for entry in entries:
    ...         [ DO STUFF ]
    """
    
    # Skip to first header
    for probeline in filehandle:
        if probeline.startswith('>'):
            break
    else: # nobreak
        raise ValueError('No headers in this file.')
    
    header = probeline.strip('>\n')
    buffer = list()
    
    # Iterate over lines
    for line in map(str.rstrip, filehandle):
        if line.startswith('>'): 
            yield header, ''.join(buffer)
            buffer.clear()
            header = line[1:]
            
        else:
            buffer.append(line)
            
    yield header, ''.join(buffer)



AssemblyStats = _collections.namedtuple('AssemblyStats', ['size', 'n50', 'ncontigs', 'largest', 'smallest', 'sizemax', 'sizestep', 'sizes'])



def assemblystats(fasta_path, xmax=10000, step=100):
    """Returns statistics about a fasta file from an assembly."""
    
    length_counter = _collections.Counter()
    
    with open(fasta_path) as filehandle:
        for header, sequence in simple_iterfasta(filehandle):
            length_counter[len(sequence)] += 1
    
    lengthcounts = sorted(length_counter.items(), reverse=True)
    
    # Initialize and calculate all other variables than "lengths" and "N50"
    assemblysize = sum(length * count for length, count in lengthcounts)
    ncontigs = sum(count for length, count in lengthcounts)
    largestcontig = lengthcounts[0][0]
    smallestcontig = lengthcounts[-1][0]
    N50 = None
    
    # Length distribution
    lengths = list()
    thresholds = reversed(range(0, xmax + 1, step))
    threshold = next(thresholds)
    currentsize = 0
    
    for length, count in lengthcounts:
        while length < threshold:
            lengths.append(currentsize)
            threshold = next(thresholds)
            
        currentsize += length * count
        
        if N50 is None and currentsize >= assemblysize/2:
            N50 = length
    
    lengths.append(currentsize)
    for threshold in thresholds:
        lengths.append(currentsize)
    
    lengths = lengths[::-1]
            
    return AssemblyStats(assemblysize, N50, ncontigs, largestcontig, smallestcontig, xmax, step, lengths)

