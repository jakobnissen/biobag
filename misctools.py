"""Bioinformatics toolbox

This is a loose collection of tools for tasks that the author have had to do
too many times, e.g. parsing and translating a Fasta file.

Author: Jakob Nybo Nissen, DTU Bioinformatics
"""


import collections as _collections


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
    ('A', 'A', 'A'): 'K', ('A', 'A', 'G'): 'K', ('A', 'A', 'T'): 'N', ('A', 'A', 'C'): 'N',
    ('A', 'G', 'A'): 'R', ('A', 'G', 'G'): 'R', ('A', 'G', 'T'): 'S', ('A', 'G', 'C'): 'S',
    ('A', 'T', 'A'): 'I', ('A', 'T', 'G'): 'M', ('A', 'T', 'T'): 'I', ('A', 'T', 'C'): 'I',
    ('A', 'C', 'A'): 'T', ('A', 'C', 'G'): 'T', ('A', 'C', 'T'): 'T', ('A', 'C', 'C'): 'T',
    ('G', 'A', 'A'): 'E', ('G', 'A', 'G'): 'E', ('G', 'A', 'T'): 'D', ('G', 'A', 'C'): 'D',
    ('G', 'G', 'A'): 'G', ('G', 'G', 'G'): 'G', ('G', 'G', 'T'): 'G', ('G', 'G', 'C'): 'G',
    ('G', 'T', 'A'): 'V', ('G', 'T', 'G'): 'V', ('G', 'T', 'T'): 'V', ('G', 'T', 'C'): 'V',
    ('G', 'C', 'A'): 'A', ('G', 'C', 'G'): 'A', ('G', 'C', 'T'): 'A', ('G', 'C', 'C'): 'A',
    ('T', 'A', 'A'): '*', ('T', 'A', 'G'): '*', ('T', 'A', 'T'): 'Y', ('T', 'A', 'C'): 'Y',
    ('T', 'G', 'A'): '*', ('T', 'G', 'G'): 'W', ('T', 'G', 'T'): 'C', ('T', 'G', 'C'): 'C',
    ('T', 'T', 'A'): 'L', ('T', 'T', 'G'): 'L', ('T', 'T', 'T'): 'F', ('T', 'T', 'C'): 'F',
    ('T', 'C', 'A'): 'S', ('T', 'C', 'G'): 'S', ('T', 'C', 'T'): 'S', ('T', 'C', 'C'): 'S',
    ('C', 'A', 'A'): 'Q', ('C', 'A', 'G'): 'Q', ('C', 'A', 'T'): 'H', ('C', 'A', 'C'): 'H',
    ('C', 'G', 'A'): 'R', ('C', 'G', 'G'): 'R', ('C', 'G', 'T'): 'R', ('C', 'G', 'C'): 'R',
    ('C', 'T', 'A'): 'L', ('C', 'T', 'G'): 'L', ('C', 'T', 'T'): 'L', ('C', 'T', 'C'): 'L',
    ('C', 'C', 'A'): 'P', ('C', 'C', 'G'): 'P', ('C', 'C', 'T'): 'P', ('C', 'C', 'C'): 'P',
    }
    
    def __init__(self, header, sequence):
        if not header or not sequence:
            raise ValueError('Header and sequence must be nonempty')
            
        if header.startswith('>'):
            object.__setattr__(self, 'header', header[1:])#self.header = header[1:]
        else:
            object.__setattr__(self, 'header', header)#self.header = header
        self.sequence = sequence
        
    def __repr__(self):
        return '<Fasta Entry {}>'.format(self.header)
    
    def __len__(self):
        return len(self.sequence)
    
    def __str__(self):
        return '>{}\n{}'.format(self.header, self.sequence)
    
    def format(self, width=60):
        """Returns the entry as a string in FASTA format with newlines every
        width-th sequence character."""
        
        sixtymers = range(0, len(self.sequence), width)
        spacedseq = '\n'.join([self.sequence[i: i+width] for i in sixtymers])
        return '>{}\n{}'.format(self.header, spacedseq)
    
    # Two entries with same header cannot co-exist in same set/dict!
    def __hash__(self):
        return hash(self.header)
    
    def __contains__(self, other):
        return other in self.sequence
    
    # Entries are compared equal by their sequence.
    def __eq__(self, other):
        try:
            return self.sequence == other.sequence
        except AttributeError:
            return self.sequence == other
        
    def __getitem__(self, index):
        return self.sequence[2]
    
    def __setattr__(self, attribute, value):
        if attribute == 'header':
            raise AttributeError('FastaEntries headers cannot be modified.')
        else:
            object.__setattr__(self, attribute, value)
    
    def reversecomplement(self, d=dict(zip('ACGTMRWSYKVHDBN', 'TGCAKYWSRMBDHVN'))):
        try:
            complemented = ''.join([d[n] for n in reversed(self.sequence)])
        except KeyError as exception:
            exception.args = (f'{exception.args[0]} is not a valid DNA base',)
            raise
            
        self.sequence = complemented
        return self
    
    def translate(self, endatstop=False):
        codons = zip(*[iter(self.sequence)] * 3)
        try:
            translated = ''.join(self.genetic_code.get(codon) for codon in codons)
        except KeyError as exception:
            exception.args = (f'{exception.args[0]} is not a valid DNA codon',)
            raise
        
        if endatstop:
            stoppos = translated.find('*')
            if stoppos == 0:
                return None
            
            elif stoppos != -1:
                translated = translated[:stoppos]
            
        self.sequence = translated
        return self


def iterfasta(filehandle, FastaEntry=FastaEntry):
    """A generator which yields FastaEntries from an open fasta file.
    
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
            yield FastaEntry(header, ''.join(buffer))
            buffer.clear()
            header = line[1:]
            
        else:
            buffer.append(line)
            
    yield FastaEntry(header, ''.join(buffer))


def simplefastaiter(filehandle):
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
            yield header, ''.join(buffer
            buffer.clear()
            header = line[1:]
            
        else:
            buffer.append(line)
            
    yield header, ''.join(buffer)


SamLineBase = _collections.namedtuple('SamLineBase', ['qname', 'flag', 'rname', 'pos',
                                         'mapq', 'cigar', 'rnext', 'pnext',
                                         'tlen', 'seq', 'qual', 'optional'])


class SamLine(SamLineBase):
    """A single alignment file in a SAM file. Instantiate with the `fromstring`
    method.
    
    >>> line = 'readname1\\t77\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\tAAAGC [...]'
    >>> samline = SamLine.fromstring(line)
    >>> samline.bits
    ['paired', 'unmapped', 'partner unmapped', 'forward']
    >>> samline.hasflag('paired') and samline.hasflag(12)
    True
    """
    
    flagdescriptions = {0x1: 'paired',
                        0x2: 'both aligned',
                        0x4: 'unmapped',
                        0x8: 'partner unmapped',
                        0x10: 'complemented',
                        0x20: 'partner complemented',
                        0x40: 'forward',
                        0x80: 'reverse',
                        0x100: 'secondary',
                        0x200: 'quality failed',
                        0x400: 'duplicate',
                        0x800: 'supplementary'}
    
    descriptionbits = {v: k for k, v in flagdescriptions.items()}
    
    def __repr__(self):
        return self.qname
    
    def __str__(self):
        return '\t'.join([str(field) for field in self])
    
    @property
    def bits(self):
        """A list of the applicable bitwise flags"""
        
        return [word for bit, word in SamLine.flagdescriptions.items()
                if self.flag & bit]
    
    def hasflag(self, flag):
        """Given a flag string, returns whether that flag applies.
        Alternatively, given an integer, returns whether all bits in the integer apply."""
        
        try:
            return flag & self.flag == flag
        except TypeError:
            try:
                return bool(SamLine.descriptionbits[flag] & self.flag)
            except KeyError as exception:
                exception.args = ((f'{exception.args[0]} is not a proper flag.'
                                   f'Choose among {", ".join(SamLine.descriptionbits)}.'),)
            raise
    
    @classmethod
    def fromstring(cls, string):
        """Creates a SamLine from a tab-separed string. This is the preferred
        way of instantiating a SamLine"""
        
        fields = string.split('\t')
        
        if len(fields) < 11:
            raise ValueError(f'Too few fields in SAM line {string}')
        elif len(fields) == 11:
            optional = ''
        else:
            optional = '\t'.join(fields[11:])
            
        for numericalindex in (1, 3, 4, 7, 8):
            fields[numericalindex] = int(fields[numericalindex])
        
        return cls(*fields[:11], optional)


class SamParser:
    """This class is for opening and parsing SAM files. You can iterate over
    the SamParser. It will yield strings when encountering a SAM header, and
    SamLine objects when encountering alignment lines
    
    >>> with open('myfile.sam') as samfile:
    ...     parser = SamParser(samfile)
    ...     parser.consumeheaders()
    ...     for samline in parser:
    ...         [ DO STUFF HERE ]
    """
    
    def __init__(self, filehandle):
        self.filehandle = filehandle
        self.state = None
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.state is not None:
            samline = self.state
            self.state = None
            return samline
        
        line = next(self.filehandle).strip()

        if line.startswith('@'):
            return line

        else:
            return SamLine.fromstring(line)
            
    def iterheaders(self):
        """Returns a generator with all the headers in the SAM file."""
        
        while True:
            if self.state is not None:
                return

            line = next(self)
        
            if isinstance(line, str):
                yield line
        
            else:
                self.state = line
                return
    
    def consumeheaders(self):
        "Consumes all headers and returns the number of headers consumed."
        return sum(1 for i in self.iterheaders())


#Right now, this class is not needed. Implement it later if I need it.
#A fastq entry can just be a string for my current needs.

class FastqEntry:
    """TODO
    include deletions in either qual or seq mirrors in the other
    """
    
    __slots__ = ['header', 'sequence', 'quality']
       
    def __init__(self, header, sequence, quality):
        if not header or not sequence or not quality:
            raise ValueError('Header, sequence and quality must be nonempty')
        
        header = header[1:] if header[0] == '@' else header
        self.header = header
        self.sequence = sequence
        self.quality = quality
        
    def check(self, phred=33):
        """Checks whether the format is OK. This is a function separate from
        __init__ in order to speed the latter up."""
        
        if len(self.sequence) != len(self.quality):
            return False
        
        elif any(ord(ch) < phred for ch in self.quality):
            return False
        
        else:
            return True
    
    @property
    def logprobs(self):
        #return [-(ord(ch)-33)/4.3429448190325175 for ch in self.quality]
        return [(33 - ord(ch)) / 4.3429448190325175 for ch in self.quality ]
        
    def __repr__(self):
        return '<Fastq Entry {}>'.format(self.header)
    
    def __len__(self):
        return len(self.sequence)
    
    def __str__(self):
        return '@{}\n{}\n+\n{}'.format(self.header, self.sequence, self.quality)
    
    # Two entries with same header cannot co-exist in same set/dict!
    def __hash__(self):
        return hash(self.header)
    
    def __contains__(self, other):
        return other in self.sequence
    
    # Entries are compared equal by their sequence and quality.
    def __eq__(self, other):
        try:
            return self.sequence == other.sequence and self.quality == other.quality
        except AttributeError as exception:
            exception.args = ('Must compare to Fastq-like object')
            raise
        
    def __getitem__(self, index):
        return self.sequence[2]


def iterfastq(filehandle, singleline=False, FastqEntry=FastqEntry):
    """A generator which yield FastqEntries from an open fastq file.
    Is about 50% faster if single line (2.6 µs vs 4 µs / entry).
    
    Usage:
    >>> with open('myfile.fastq') as fastqfile:
    ...     entries = iterfastq(fastqfile)
    ...
    ...     for entry in entries:
    ...         [ DO STUFF ]
    """

    if singleline:
        strippedlines = map(str.strip, filehandle)
        for header, sequence, plus, quality in zip(*[strippedlines]*4):
            yield FastqEntry(header, sequence, quality)
            
    else:
        seqbuffer = list()
        qualbuffer = list()
        readingquality = False
        
        header = next(filehandle).strip()
        if not header.startswith('@'):
            raise ValueError('First line is not a header')
        header = header[1:]
        
        for line in map(str.rstrip, filehandle):
            if line.startswith('@'):
                yield FastqEntry(header, ''.join(seqbuffer), ''.join(qualbuffer))
                qualbuffer.clear()
                seqbuffer.clear()
                header = line[1:]
                readingquality = False
                
            elif line.startswith('+'):
                readingquality = True

            elif readingquality:
                qualbuffer.append(line)

            else:
                seqbuffer.append(line)
                
        yield FastqEntry(header, ''.join(seqbuffer), ''.join(qualbuffer))


def assemblystats(inputstring, thresholds=[10000, 5000, 1000]):
    """Expects the output from following command as input:
    For SPAdes: grep '>' contigs.fasta | cut -d_ -f4 | sort -nr | uniq -c
    For MEGAHIT: grep '>' final.contigs.fa | cut -d= -f4 | sort -nr | uniq -c
    
    Returns n_contigs, largest_contig, N50, n_nucleotides_above_each_threshold,
    total nucleotides."""
    
    lines = (line.split() for line in inputstring.splitlines())
    counter = tuple((int(length), int(count)) for count, length in lines)
    
    if not 0 in thresholds:
        thresholds = sorted(list(thresholds) + [0], reverse=True)
    else:
        thresholds = sorted(thresholds, reverse=True)
    
    thresholditer = iter(thresholds)
    threshold = next(thresholditer)
    
    nucleotides = list()
    runningcount = 0
    n50 = 0
    
    assemblylength = sum(length*count for length, count in counter)
    ncontigs = sum(entry[1] for entry in counter)
    largestcontig = counter[0][0]
    
    for length, count in counter:
        while length < threshold:
            nucleotides.append(runningcount)
            threshold = next(thresholditer)
            
        runningcount += length * count

        if runningcount > assemblylength/2 and n50 == 0:
            n50 = length
            
    for remainingthreshold in [threshold] + list(thresholditer):
        if remainingthreshold < length or not nucleotides:
            nucleotides.append(runningcount)
        else:
            nucleotides.append(nucleotides[-1])
            
    return (ncontigs, largestcontig, n50, *nucleotides)


### Need to check all the new functionality of FastaEntry
# this is its hasing, __eq__, setattr stuff, __contains__, __getitem__, translate method
# check fastq entry and its iterator

if __name__ == '__main__':
    # Unit test of significant
    assert significant(543.21) == '543'
    assert significant(-5, 2) == '-5.0'
    assert significant(0.004, 2) == '0.0'
    assert significant(0.1 + 0.2) == '0.30'
    
    # Test SamLine class
    testline =  ('NS500333:93:H7NV2BGX2:4:23612:14579:20391	77	*	0	0	'
                '*	*	0	0	AAAGCACTCAGATGTTGTTTTCAAGGTCACAAATTACATCAAGAA'
                'GAATTGTGGGGACAAGCTGTCGCTCGACACACTCGCACGGGAAGTCTATCTCAGCAAATCA'
                'TATCTCAGCAGTATTTTTAAGGAAGAGACAGGCATCGGACTCA	AAAAAEEEEEEEEE'
                'EEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEA<AEEEEEEEEEAEAAEEE'
                '<A<AEEEEEAEEAEAEEEEAEEEEE/EAEEAEEEE/AAAAAEE/EE//E/AE/EE/EA/A'
                'A<AAEEAEAEE<A/	AS:i:0	XS:i:0')
    
    samline = SamLine.fromstring(testline)
    
    assert samline.__repr__() == 'NS500333:93:H7NV2BGX2:4:23612:14579:20391'
    assert samline.flag == 77
    assert samline.bits == ['paired', 'unmapped', 'partner unmapped', 'forward']
    assert samline.optional == 'AS:i:0	XS:i:0'
    assert samline.hasflag('paired')
    assert not samline.hasflag('reverse')
    assert samline.hasflag(4)
    assert samline.hasflag(12)
    
    # Test SamParser
    with open('test.sam') as samfile:
        parser = SamParser(samfile)
        
        assert next(parser).startswith('@')
        assert isinstance(next(parser), str)
        
        # Skip all the headers
        while isinstance(next(parser), str):
            pass
        
        samlines = list(parser)
        assert len(samlines) == 99 # We skipped first nonheader line
        
    # Test FastaEntry
    newentry = FastaEntry('myheader | sup dawg', 'GTAGTCGATAGAAGTAGTGCTTCTA')
    
    assert newentry.__repr__() == '<Fasta Entry myheader | sup dawg>'
    assert str(newentry) == '>myheader | sup dawg\nGTAGTCGATAGAAGTAGTGCTTCTA'
    
    complemented = newentry.reversecomplement()
    
    assert complemented.__repr__() == '<Fasta Entry myheader | sup dawg>'
    assert str(complemented) == '>myheader | sup dawg\nTAGAAGCACTACTTCTATCGACTAC'
    
    with open('test.fasta') as f:
        parser = iterfasta(f)

        fastaentries = list(parser)
    
    # Check that all the characters got included
    chars = 75 # newlines
    for entry in fastaentries:
        chars += len(entry)
        chars += len(entry.header) + 1 # include the > sign
    assert chars == 18096

