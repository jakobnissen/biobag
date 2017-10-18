"""================ Jakob's SwissProt/UniProt interface script ===============
Examples:
Make a local, searchable database with all of Swissprot:
(note: database is created in memory, then moved to disk)
    [jakni@computerome] wget [SWISSPROT FTP ADDRESS HERE] swissprot.dat.gz
    [jakni@computerome] python
    >>> import unifetch
    >>> unifetch.createdb('swissprot.gdbm', 'swissprot.dat.gz')

Search the database locally by primary citable accession number
    >>> accessions = ["D4ABB2", "Q53XC0", "Q4R3D3"]
    >>> records = [unifetch.fetch(acc, 'swissprot.gdbm') for acc in accessions]
    >>> for record in records: print(record.entry_name, record.sequence)
    
Look-up accessions remotely (any accession number)
    >>> record = unifetch.webfetch("Q75MH2")
    
Pipe through the script to annotate tabular blast output (makes extra column):
    [jakni@computerome] cat blastoutput | python unifetch.py > annotated
    [jakni@computerome] head -1 annotated | cut -f7-
    46	314	19	288	7.0e-77	288.9	Eukaryota;Metazoa;Chord [ ... ]
        
Invoke the script from command line with no flags to dump an entire entry:
    [jakni@nissen scripts]$ python unifetch.py -a A6XGL2
                [ ... SWISSPROT ENTRY ELIDED ... ]

If you use flags, it only displays specific information:
    [jakni@nissen scripts]$ python unifetch.py -a A6XGL2 -ncis
    Name: A6XGL2_HUMAN
    Data class: Unreviewed
    TaxID: 9606
    Sequence: MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRRE [ ... ]

Requirements:
Biopython
On computerome, run with anaconda3/4.4.0 (only version where dbm.gnu works)
==============================================================================
"""

__version__ = 1,0


# Maintenance of this script
# To create a new database, download the newest uniprot/swissprot text file
# (.dat.gz), then run createdb. It requires around 175 GB of RAM for uniprot.
# The reason it only works with module anaconda3/4.4.0 is a bug with anaconda
# on CentOS which causes the library dbm.gnu to not be able to be imported.
# The sysadmins fixed it for anaconda3/4.4.0 and can probably do it for other
# versions as well.


import sys as _sys
import os as _os
import io as _io
import gzip as _gzip

try:
    import dbm.gnu as _gnu
except ModuleNotFoundError as error:
    message = ('Cannot find module dbm.gnu. Make sure computerome module '
               'anaconda3/4.4.0 is loaded. If it is, contact computerome support.')
    raise ModuleNotFoundError(message) from error

from Bio import SwissProt as _SwissProt
from Bio import ExPASy as _ExPASy


SWISSPROTPATH = '/home/projects/metagenomics/data/uniprot/swissprot.gdbm'
UNIPROTPATH = '/home/projects/metagenomics/data/uniprot/uniprot.gdbm'


def createdb(outfilepath, infilepath):  
    """Creates a new database from a SwissProt/UniProt text file, gzipped or not.
    For speed, database is built in memory, then moved to disk. Takes ~11 hrs."""
    
    import shutil as _shutil
    
    if _os.path.exists(outfilepath):
        raise FileExistsError('Database already exists.')
    
    # Check whether the database is gzipped or not by searching for the two
    # signature bytes 1F8B and use gzip.open if it is.
    with open(infilepath, 'rb') as infile:
        signature = infile.read(2)
    
    if signature == b'\x1f\x8b':
        opener = _gzip.open
    else:
        opener = open
    
    # Read the content of the text file. At accession identifier, extract accession.
    # at end of record, save the current record under extracted accession ID.
    # Create a database in memory.
    accession = None
    buffer = list()
    tempfilename = '/dev/shm/temp.gdbm'
    with opener(infilepath, 'rt') as infile, _gnu.open(tempfilename, 'cf') as db: 
        for line in infile:
            buffer.append(line)
            
            if line.startswith('//'):
                assert accession is not None
                db[accession] = _gzip.compress(bytes(''.join(buffer), 'ASCII'))
                buffer.clear()
                accession = None

            elif line.startswith('AC') and accession is None:
                accession = line.split()[1][:-1]
        
        # Because I openened the database in fast mode, I need to sync before closing.
        db.sync()
        
    # Move file from memory to actual file location
    _shutil.move(tempfilename, outfilepath)


def fetch(accession, swissprot=False):
    """Given an primary accession string, returns a protein
    record if accession is present in database, else None.
    
    SwissProt as about 4x faster to access than UniProt."""
    
    dbpath = SWISSPROTPATH if swissprot else UNIPROTPATH

    try:
        with _gnu.open(dbpath) as database:
            data = database.get(accession, None)
    except _gnu.error as error:
        raise FileNotFoundError('Database not found: {}'.format(dbpath)) from error

    if data is None:
        return None

    # This part here is awkward, but Biopython only parses files, bytes/str.
    # so I use BytesIO to convert bytes to a file-like object.
    else:
        fileobject = _io.BytesIO(_gzip.decompress(data))
        return _SwissProt.read(fileobject)


def webfetch(accession):
    "Fetches a SwissProt entry from the internet given an accession string."
    
    with _ExPASy.get_sprot_raw(accession) as fileobject:
        return _SwissProt.read(fileobject)


def get_species(entry, discard={'sp', 'spp', 'bacterium', 'cf', 'archaeon',
                                'genomosp', 'parasite', 'endosymbiont'}):
    """Attempts to return the cleaned species name of an entry. Error-prone.
    
    Returns an empty string if the organism field seem to not contain
    a species name"""
    
    fields = [field.strip('.') for field in entry.organism.split()]
    
    # E.g. "Campylobacter."
    if len(fields) < 2:
        return ''
    
    # This is the case for well-characterized but uncultured bacteria.
    if fields[0] in ('Candidatus', 'Ca.'):
        genus, species = fields[1:3]
    else:
        genus, species = fields[:2]

    # if e.g. "uncultured virus sample", "environmental sample", etc.
    if not genus.istitle():
        return ''
    
    # Uncertain geni are sometimes square bracketed
    if genus.startswith('['):
        return ''
    
    # Vira have special names, like "Tobacco mosaic virus". They typically
    # end in "virus" or "phage" though. This is error-prone, and does not
    # cover e.g. the species "Salmonella virus SP6".
    for fieldno, field in enumerate(map(str.lower, fields)):
        if field.endswith('phage') or field.endswith('virus'):
            if len(fields) == fieldno + 1:
                return ' '.join(fields)
            
            else:
                if fields[fieldno + 1].isalnum():
                    return ' '.join(fields[:fieldno + 2])
                
                else:
                    return ' '.join(fields[:fieldno + 1])
    
    if species in discard:
        return ''
    
    if not species.islower(): # does not apply to vira
        return ''
    
    return species


def _shell_lookup(args):
    """This function is called when the script is used from command line:
    
    [jakni@nissen scripts]$ python unifetch.py -a A6XGL2 -ncis
    Name: A6XGL2_HUMAN
    Data class: Unreviewed
    TaxID: 9606
    Sequence: MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRRE [ ... ]
    """
    
    with _gnu.open(args.database) as database:
        data = database.get(args.accession, None)

    # If no accession is found, return "Not found."
    if data is None:
        return 'Not found.'

    fields = {'Name': [args.name],
                 'Date': [args.date],
                 'Data class': [args.dataclass],
                 'Organism': [args.organism],
                 'Taxonomy': [args.taxonomy],
                 'TaxID': [args.taxid],
                 'Sequence': [args.sequence]
              }

    # If nothing particular is specified, return the entire accession
    if not any(arr[0] for arr in fields.values()):
        text = _gzip.decompress(data).decode()
        return text
    
    else:
        # If output specified, return the relevant parts.
        fileobject = _io.BytesIO(_gzip.decompress(data))
        record = _SwissProt.read(fileobject)

        fields['Name'].append(record.entry_name)
        fields['Date'].append(record.created[0])
        fields['Data class'].append(record.data_class)
        fields['Organism'].append(record.organism)
        species = get_species(record)
        fields['Taxonomy'].append(
            ';'.join(record.organism_classification + ([species] if species else [])))
        fields['TaxID'].append(';'.join(record.taxonomy_id))
        fields['Sequence'].append(record.sequence)

        output = list()
        for title, (state, information) in fields.items():
            if state:
                output.append('{}: {}'.format(title, information))
        return '\n'.join(output)


def _shell_blast_annotate(inputhandle):
    """This generator is called when piping a blast output through the script.
    
    It is assumed that output's second column is of format "X|ACCESSION",
    where X may be any string.
    
    [jakni@computerome] cat blastoutput | python unifetch.py > annotated
    [jakni@computerome] head -1 annotated | cut -f7-
    46\t314\t19\t288\t7.0e-77\t288.9\tEukaryota;Metazoa;Chord [ ... ]"""
    
    for line in inputhandle:
        accession = line.split()[1].split('|')[-1]
        record = fetch(accession, swissprot=args.swissprot)

        if record is None:
            end = 'Unknown'
        else:
            species = get_species(record)
            end = ';'.join(record.organism_classification + ([species] if species else []))

        yield line.strip() + '\t' + end


# Use the script from command line
if __name__ == '__main__':
    # First add the relevant command line arguments
    import argparse
    
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('--swissprot', action='store_true', dest='swissprot',
                        help='use swissprot instead of uniprot')
    
    parser.add_argument('-a', dest='accession', help='accession ID to look up')
    
    parser.add_argument('-n', action='store_true', dest='name', help='return entry name')
    parser.add_argument('-s', action='store_true', dest='sequence', help='return sequence')
    parser.add_argument('-d', action='store_true', dest='date', help='return creation date')
    parser.add_argument('-c', action='store_true', dest='dataclass', help='return data class')
    parser.add_argument('-t', action='store_true', dest='taxonomy', help='return taxonomy')
    parser.add_argument('-o', action='store_true', dest='organism', help='return organism')
    parser.add_argument('-i', action='store_true', dest='taxid', help='return taxid')
    
    args = parser.parse_args()
    
    if args.swissprot:
        args.database = SWISSPROTPATH
    else:
        args.database = UNIPROTPATH
    
    # Check if database exists
    if not _os.path.exists(args.database):
        raise FileNotFoundError('Database not found: {}'.format(args.database))
    
    # If an accession is given, the script is used to lookup an accession
    if args.accession:
        print(_shell_lookup(args))
    
    # If no accession is given assume a BLAST output is piped through the script.
    else:
        for line in _shell_blast_annotate(_sys.stdin):
            print(line)


unittest = False

# Test get_species
if unittest:
    class A:
        def __init__(self, x):
            self.organism = x

    names = [ # Names to be discarded
             ('unclassified phage', ''),
             ('marine sediment metagenome', ''),
             ('Pseudomonas bacterium', ''),
             ('Ricketta sp.', ''),
             ('Toxoplasmosis cf. (from sea water)', ''),
             ('[Escerichia] coli', ''),
             ('Thermococcales archaeon 44_46', ''),
             ('Lactobacillus genomosp. strain 19', ''),
             ('Wilfordii endosymbiont of common tick', ''),
             ('Feline parasite (Nematoda sp.)', ''),
             ('Candidatus Uhrbacteria bacterium GW2011_GWE2_45_35', ''),
             ('Fusicatenibacter.', ''),
              # Viral names
             ("Frog virus SG2", "Frog virus SG2"),
             ('Human immunodeficiency virus 1.', 'Human immunodeficiency virus 1'),
             ("Tobacco Mosaic Virus", "Tobacco Mosaic Virus"),
             ('Lambda phage (SG21)', "Lambda phage"),
             ('T2 phage [2017] (E coli bacteriophage).', "T2 phage"),
              # Regular species names
             ("Homo sapiens (Human).", 'sapiens'),
             ('Lupus canis vulgaris.', 'canis'),
             ('Bos taurus.', 'taurus'),
             ('Candidatus Kryptonium thompsoni.', 'thompsoni'),
            ]
    for name, result in names:
        if not get_species(A(name)) == result:
            print('Input:', name, 'Expected:', result, 'Result:', get_species(A(name)))
            raise AssertionError

