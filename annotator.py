"""Takes a protein query file and a "bin table" file and diamondblasts the query
against swissprot, then summarizes which species each cluster hit. You can get
the two input files from the scripts parsebins.py or parsecanopy.py

Is based on a crude approach to be refined later.
"""


__author__ = 'Jakob Nybo Nissen, DTU Bioinformatics'


import os
import sys
import subprocess
import argparse
from math import pi, exp, sqrt
from collections import Counter, defaultdict
from itertools import zip_longest
from time import time


nanvalues = ('environmental samples', 'unclassified sequences', 'unclassified phages')
taxranks = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']


def timed(function):
    """Decorator adding a timer to a function,
    and prints the time elapsed in the terminal. Just eye candy."""
    
    def inner(*args, **kwargs):
        begin = time()
        result = function(*args, **kwargs)
        print('\tDone in {:,.2f} seconds'.format(time() - begin))
        return result
    
    return inner


@timed
def diamond_blast(args):
    "Runs diamondblast."
    
    print('Diamondblasting with {} cores.'.format(args.cores))
    outputfields = ('qseqid sseqid pident length mismatch gapopen qstart qend '
                    'sstart send evalue bitscore stitle')
    log = args.diamondout + '.log'
    
    command = '{} blastp -q {} -d {} -f 6 {} -o {} -k 1 -p {} > {}'
    command = command.format(args.diamondpath, args.query, args.database,
                            outputfields, args.diamondout, args.cores, log)
    
    subprocess.run(command, shell=True, check=True)
    
    return None


@timed
def create_bindict(bintablepath):
    """Reads a canopy cluster file and creates the bindict"""
    
    print('Parsing clusters.')

    bindict = defaultdict(list)
    
    with open(bintablepath) as file:
        currentbin, gene = next(file).split()
        currentbinbuffer = [(currentbin, gene)]
        
        # Iterate over all nonempty lines
        for bin, gene in filter(None, map(str.split, file)):
            if bin != currentbin:
                for bufbin, bufgene in currentbinbuffer:
                    bindict[bufbin].append(bufgene)
                
                currentbinbuffer.clear()
                currentbin = bin
            
            currentbinbuffer.append((bin, gene))
        
        for bufbin, bufgene in currentbinbuffer:
            bindict[bufbin].append(bufgene)

    return bindict


def scorer(x, rank, tau=2*pi, exp=exp, sqrt=sqrt):
    """Converts amino acid identity < 1 to a score.
    
    It's the sum of three scores:
    1) A constant base score for just getting the hit
    2) A low-homology score, which contributes in the middle range (0.4-0.8)
        Normal distribution with mean of 1 and std. deviation of sigma
    3) A high-homology score which dominates for identities > 0.9
        Sigmoid function with max x: f'(x) = 1 and f'(1) = homologyscaling
        
    It's calibrated so that scorer(x, 7) = sigmoidfunction(x), and
    scorer(1, rank) = 1 for any rank
    """
    
    # A hit with 0% identity gets this score
    base = max(0, 0.03 - 0.005 * rank)
    
    
    # The slope of high-homology part at 100% identity
    homologyscaling = 6 + rank
    
    # The relative contribution of high vs low homology at 100% identity
    highlowfraction = min(1, 0.4 + 0.1 * rank)
    
    # This sets the relative spread of low homology
    sigma = 0.3 - (0.03 * rank)
    
    # Calculate score
    if base == 0 and highlowfraction == 1:
        return 2 / (1 + exp(-4*homologyscaling * (x - 1)))
    
    sigmasquared = sigma * sigma
    
    highfraction = (1 - base) * highlowfraction
    lowfraction = 1 - base - highfraction
    
    highhomology = 2 * highfraction / (1 + exp(-4*homologyscaling * (x - 1)))
    
    lowscale = lowfraction * sqrt(tau*sigmasquared)
    lowhomology = (exp(-((x-2)*x+1) / (2*sigmasquared)) / sqrt(tau*sigmasquared)) * lowscale
    
    return lowhomology + highhomology + base


def naively_summarize_bin(taxidpairs, scorer=scorer, introspect=False, nanvalues=nanvalues):
    """Given a tuple of taxonomystring, score, returns a consensus taxonomy string
    for those genes. Nanvalues, e.g. 'unknown sequence' are filtered out.
    
    Current rules:
    1) The score is calculated by scorer, and depends on the taxonomic rank
        (higher is more specific)
    2) The score must >= 5
    3) The score must be double that of the next-best in that rank
    4) The score must be >= 0.4 times the sum of any assignment at that rank.
    
    As this function is called thousands of times, do not decorate it
    with @timed"""
    
    consensus = list()
    taxidpairs = [(t.split(';'), i) for t, i in taxidpairs]
    taxidpairs = [(t, i) for t, i in taxidpairs if t[0] not in nanvalues]
    
    # Do not take consensus of too few genes.
    if len(taxidpairs) < 10:
        return 'None'
    
    # Get a score for each hit for each level
    taxscorepairs = list()
    for taxlist, id in taxidpairs:
        taxscorepairs.append([(tax, scorer(id, rank)) for rank, tax in enumerate(taxlist)])
    
    # Get a global score for each rank (i.e. sum of all scores at that rank)
    ranks = zip_longest(*taxscorepairs, fillvalue=(None, 0))
    globalscores = [sum(score for tax, score in rank) for rank in ranks]

    # Each rank corresponds to kingdom, phylum, subphylum, superclass, etc.
    for rank, globalscore in enumerate(globalscores):
        scores = Counter()
        
        for tsp in taxscorepairs:
            tax, score = tsp[rank]
            scores[tax] += score
        
        counter = Counter([tsp[rank][0] for tsp in taxscorepairs])
        
        for nanvalue in nanvalues:
            scores.pop(nanvalue, None)
            
        mostcommons = scores.most_common()
        
        # If mostcommons is empty, i.e. only nanvalues in this rank
        # stop and do not attempt for any lower rank.
        if not mostcommons:
            break
        
        maxtaxon, maxscore = mostcommons[0]
        maxcount = counter[maxtaxon]
        secondscore = 0 if len(mostcommons) == 1 else mostcommons[1][1]
            
        # Introspect if given the command
        if introspect:
            
            
            print('{} (rankwide score: {})'.format(taxranks[rank], round(globalscore, 1)))
            for row in range(min(3, len(mostcommons))):
                taxon, score = mostcommons[row]
                shortened = taxon[:16] + '...' if len(taxon) > 19 else taxon
                print('\t{:<20}{:<8}{}'.format(shortened, round(score, 2), counter[taxon]))
            print('')
            
        # If a consensus can be found, add it to the output.
        criteria = (maxscore >= 5,
                    maxscore >= 1.8 * secondscore,
                    maxscore >= 0.4 * globalscore,
                    #maxscore / maxcount >= 0.25,
                    )
        
        if all(criteria):
            consensus.append(maxtaxon)
            
            # For next iteration, discard any taxonomies that was not part of
            # the consensus, or where we have reached the end of the taxonomy.
            taxscorepairs = [tsp for tsp in taxscorepairs
                             if tsp[rank][0] == maxtaxon and len(tsp) >= rank+2]

        # If no consensus found, do not attempt for any smaller clade
        else:
            if introspect:
                reasons = ('Maxscore too low',
                           'Second place too good',
                           'Maxscore too small fraction of whole',
                           'Average score too low')
                reasons = ', '.join(r for r, c in zip(reasons, criteria) if not c)
                print('Ended:', reasons + '.', end='\n\n')
            break
                          
    if consensus:
        return ';'.join(consensus)
    
    else:
        return 'None'


@timed
def naive_summary(bindict, diamondout, returntax=False):
    """Given the annotated blast output and a bindict as
    returned by parse_clusters, creates summary of the bin."""
    
    # Create dictionary of gene: taxonomy.
    print('Summarizing bin taxonomy.')
    print('\tExtracting taxonomy...')
    taxdict = dict()
    with open(diamondout) as blastfile:
        for line in blastfile:
            if not line: # do not crash if file ends with empty line
                continue
            
            gene, _, identity, *__, taxonomy = line.split('\t')
            
            # Remove gene name and linebreak from taxonomy
            taxonomy = taxonomy.partition(' ')[2].rstrip()
            identity = float(identity) / 100
            taxdict[gene] = (taxonomy, identity)
    
    if returntax:
        return taxdict
    
    print('\tWriting output...')
    buffer = list()
    for binname, genes in bindict.items():
        taxscorepairs = [taxdict[gene] for gene in genes if gene in taxdict]
        summary = naively_summarize_bin(taxscorepairs)
        buffer.append('{}\t{}'.format(binname, summary))
        
    with open(args.output, 'w') as outfile:
        print('\n'.join(buffer), file=outfile)
        
    return None


def meditate(bintablepath, diamondout):
    "Computes the objects necessary for introspection and keeps them in memory."
    
    bindict = create_bindict(bintablepath)
    taxdict = naive_summary(bindict, diamondout, returntax=True)
    return bindict, taxdict


def introspect(binname):
    "When having meditated, you can introspect the decision process for a bin."
    
    genes = bindict[binname]
    taxidispairs = [taxdict[gene] for gene in genes if gene in taxdict]
    return naively_summarize_bin(taxidispairs, introspect=True)


#path = '/mnt/computerome/people/jakni/archaea/archaeaout/metabat/'
#bindict, taxdict = meditate(path+'bintable.txt', path+'fields')


#genes = bindict['bin.11']
#taxes = [taxdict[gene] for gene in genes if gene in taxdict]


def main(args):
    begin_time = time()
    
    if not os.path.exists(args.diamondout):
        diamond_blast(args)
    else:
        print('Diamond output found, skipping creation.')
        
    bindict = create_bindict(args.bintable)
    naive_summary(bindict, args.diamondout)
    
    print('Complete. Total elapsed time:', round(time() - begin_time, 2), 'seconds.')


if __name__ == '__main__':
    usage = 'python annotator.py bintable output diamondout [-q query] [options]'
    
    # Parse input
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=usage)
    
    cpus = os.cpu_count()
    
    parser.add_argument('bintable', help='path to table of bins')
    parser.add_argument('output', help='path to output file')
    parser.add_argument('diamondout', help='path to put/find diamond output')
    
    parser.add_argument('-q', dest='query', metavar='query', help='path to blast query')
    parser.add_argument('-b', metavar='diamond', default='/services/tools/diamond/0.8.31/diamond',
                        dest='diamondpath', help='path to diamond blast executable. [preset]')
    parser.add_argument('-d', dest='database', metavar='database',
            default='/home/projects/metagenomics/data/uniprot/swissprot.dmnd',
                        help='path to diamond database. [preset]')
    parser.add_argument('-c', dest='cores', metavar='cores', type=int, default=cpus,
        help='no. of cores to use for diamondblasting [{}]'.format(cpus))
    
    # Print help if no arguments are given
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    args = parser.parse_args()
    
    # Check number of cores
    if args.cores < 1:
        raise argsparse.ArgumentTypeError('Zero or negative cores provided. Exiting')
        
    # Input must exist - bintable
    if not os.path.isfile(args.bintable):
        raise FileNotFoundError(args.bintable)
        
    # Input must exist - diamondout, database and diamondpath
    if not os.path.exists(args.diamondout):
        for path in (args.query, args.diamondpath, args.database):
            if path is None: # only applies to values without defaults
                raise ValueError('Must specify query, unless diamondout already exists.')
            if not os.path.isfile(path):
                raise FileNotFoundError(path)
                
        if os.path.exists(args.diamondout + '.log'):
            raise FileExistsError(args.diamondout + '.log')
        
    # Output cannot exist
    if os.path.exists(args.output):
        print('Output already exists.')
        raise FileExistsError(args.output)

    main(args)

