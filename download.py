
# coding: utf-8

# ## Problems
# The field which Study.write_config fills in 'references' is arbitrary. Fix it.
# 
# When eventually making the script that calls snakemake, make sure to run it with Anaconda3, not 2
# 
# ## Test accession numbers
# PRJEB2772 Few small files
# 
# PRJEB1786 Lots of large files
# 
# PRJEB2773 Not Illumina

# # This is the newest version of Download, made to be compatible with my own Snakefile.

# In[9]:

import os
import urllib.request
import shutil
import argparse
import json
import subprocess
import gzip

from multiprocessing import Pool
from time import sleep
from collections import OrderedDict


# In[10]:

# Globally count errors and warnings
errors = 0
warnings = 0


# In[11]:

# ENA classes
class Entry:
    """A parent class for ENA entries of any kind."""
    
    def superentry(self, parent, parentdict):
        """Define the parent entry of the entry and add it to the parentdict"""
        if self.accession in parentdict:
            raise ValueError('{} already in {}'.format(self, parent))
        
        parentdict[self.accession] = self
        return parent
    
    def __repr__(self):
        return '<{} {} at {}>'.format(
            self.__class__.__name__, self.accession, hex(id(self)))
    
    def __str__(self):
        return self.accession


# In[12]:

class Study(Entry):
    """An ENA entry corresponding to one study, AKA project."""
    
    def __init__(self, accession):
        self.accession = accession
        
        # Sub-entries
        self.samples = {}
        self.experiments = {}
        self.runs = {}
        self.files = {}
        
    def write_config(self, configpath, log):
        """Creates a JSON config file for Simon's script."""
        global errors

        samples_to_experiments = {sample.accession: list(sample.experiments)
                                 for sample in self.samples.values()}

        experiments_to_runs = {experiment.accession: list(experiment.runs)
                                 for experiment in self.experiments.values()}

        runs_to_files = {run.accession: [file.path for file in 
                                         run.files.values() if file.downloaded]
                                    for run in self.runs.values()}

        encodings = set(file.guess_phred_encoding()
                        for file in self.files.values())
        
        if len(encodings) == 1:
            encoding = encodings.pop()
        else:
            encoding = 0
        
        if encoding == 0:
            errors += 1
            print('Error: Quality encoding for study {} cannot be determined. '
                  'Set to 0'.format(self.accession), file=log)
        
        # Gather dictionaries and print in JSON file.
        jsoncontent = OrderedDict((
                ('human_reference', 'UNDEFINED'),
                ('cores', 4),
                ('phred_encoding', encoding),
                ('bwa_path', '/services/tools/bwa/0.7.15/bwa'),
                ('samtools_path', '/services/tools/bwa/0.7.15/samtools'),
                ('megahit_path', '/services/tools/megahit/1.0.4-beta/megahit'),
                ('prodigal_path', '/services/tools/prodigal-2.6.2/prodigal'),
                ('adapterremoval_path', '/services/tools/adapterremoval/2.2.0/'
                                        'bin/AdapterRemoval'),
                ('cd-hit_path', '/services/tools/cd-hit-4.6.1/bin/cd-hit-est'),
                ('canopy_path', '/home/projects/pr_99009/people/sira/share/'
                                'bin/cc.bin'),
                ('samples', samples_to_experiments),
                ('experiments', experiments_to_runs),
                ('runs', runs_to_files)))                

        with open(configpath, 'w') as configfile:
            print(json.dumps(jsoncontent, indent=4), file=configfile)

        return None


# In[13]:

class Sample(Entry):
    """An ENA entry corresponding to one sample."""
    
    def __init__(self, accession, study):         
        self.accession = accession
        
        # Super-entries
        self.study = self.superentry(study, study.samples)
        
        # Sub-entries
        self.experiments = {}
        self.runs = {}
        self.files = {}


# In[14]:

class Experiment(Entry):
    """An ENA entry corresponding to one experiment, AKA library."""
    
    def __init__(self, accession, sample, platform, library_layout):         
        self.accession = accession
        self.platform = platform
        self.se = library_layout.lower() == 'single'
        
        # Super-entries
        self.sample = self.superentry(sample, sample.experiments)
        self.study = self.superentry(sample.study, sample.study.experiments)
        
        # Sub-entries
        self.runs = {}
        self.files = {}


# In[15]:

class Run(Entry):
    """An ENA entry corresponding to one run, AKA unit."""
    
    def __init__(self, accession, experiment, ftps, filesizes):
        self.accession = accession
        
        # Super-entries
        self.experiment = self.superentry(experiment, experiment.runs) 
        self.sample = self.superentry(experiment.sample, experiment.sample.runs)
        self.study = self.superentry(experiment.study, experiment.study.runs)
        
        # Sub-entries
        self.files = {}
        
        # Instantiate Files
        ftps = ['ftp://' + ftp for ftp in ftps.split(';')]
        filesizes = [int(filesize) for filesize in filesizes.split(';')]
    
        if not (len(ftps) == 1 or (not experiment.se and len(ftps) == 2)):
            raise ValueError('Unexpected number of files in run {}.'.format(
                                                                accession))
        
        interleaved = len(ftps) == 1 and not experiment.se
        for ftpadress, filesize in zip(ftps, filesizes):
            path = os.path.abspath(os.path.join(args.destination, 
                                                self.study.accession,
                                                os.path.basename(ftpadress)))
            
            File(path, self, ftpadress, filesize, interleaved)


# In[16]:

class File(Entry):
    """A gz file. Strictly speaking this is not an ENA entry, but that makes
    my code cleaner."""
    
    def __init__(self, path, run, ftpadress, filesize, interleaved):
        self.path = path
        self.accession = path # to allow superentry method
        self.ftpadress = ftpadress
        self.filesize = filesize
        self.downloaded = False
        self.interleaved = interleaved
        
        # Super-entries
        self.run = self.superentry(run, run.files) 
        self.experiment = self.superentry(run.experiment, run.experiment.files) 
        self.sample = self.superentry(run.sample, run.sample.files)
        self.study = self.superentry(run.study, run.study.files)
        
    def validate_download(self):
        if os.path.getsize(self.path) != self.filesize:
            os.remove(self.path)
            raise OSError('Size of file {} does not match expected size.'
                          'Deleting file.'.format(self.path))

        else:
            self.downloaded = True
    
    def split(self, *newpaths):
        """If the file is split into several new files, update instances."""
        for path in newpaths:
            newfile = File(path, self.run, None, os.path.getsize(path), False)
            newfile.downloaded = True
        self.run.files.pop(self.accession)
        
    def guess_phred_encoding(self):
        """Guess the PHRED encoding (33 or 64) from the first 100 lines."""
        
        with gzip.open(self.path, 'rt') as file:
            charset = set()
            try:
                for read in range(25):
                    next(file) # Header
                    next(file) # Sequence
                    next(file) # +
                    charset.update(set(next(file))) # PHRED encoded line

            except StopIteration:
                return 0

        low = any(ch in charset for ch in """!"#$%&'()*+,-./0123456789""")
        high = any(ch in charset for ch in """KLMNOPQRSTUVWXYZ[\]^_`abcdefg""")

        if low and not high:
            return 33
        elif high and not low:
            return 64
        else:
            return 0
        
    def download(self, disable_deinterleave):
        """Downloads a file given the data from a File."""
        # This function is to be executed in parallel child processes.
        # Therefore, it only has access to copies of objects.
        # So no direct modifying of objects in this function.

        if os.path.exists(self.path):
            if self.filesize == os.path.getsize(self.path):
                return None
            else:
                raise FileExistsError('{} already exist, but its size differs '
                          'from metadata file size.'.format(self.path))

        with urllib.request.urlopen(
            self.ftpadress, timeout=30) as infile, open(self.path, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)

        # Deinterleave if appropriate
        if self.interleaved and not disable_deinterleave:
            forward, reverse = self.path[:-3]+'_1.gz', self.path[:-3]+'_2.gz'
            if os.path.exists(forward) or os.path.exists(reverse):
                raise FileExistsError('Cannot disinterleave {}, '
                        'file with same name and _1 or _2 suffix'
                        'already exist'.format(self.path))

            # Deinterleave bash command derived from Nathan Watson-Haigh's work.
            command = (r'gunzip -c {} | paste - - - - - - - - | tee >('
                       r'cut -f 1-4 | tr "\t" "\n" | gzip > {}) | '
                       r'cut -f 5-8 | tr "\t" "\n" | gzip > {}').format(
                self.path, forward, reverse)


            process = subprocess.Popen(command, shell=True,
                                       executable='/bin/bash')
            process.wait() # Unbelievably, this is required.

            if process.poll() == 0: # Error code returned, 0 means no error
                os.remove(self.path)
                return forward, reverse
            else:
                return False
        else:
            return None


# In[17]:

def mkdir(path):
    """Makes a directory at the given path."""
    
    try:
        os.mkdir(path)

    except FileExistsError:
        if not os.path.isdir(path):
            print('Non-directory file already exist at: {}'.format(path))
            raise


# In[18]:

def download_entries(accession):
    """Given a study accession returns a generator of tuples:
(study_accession, sample_accession, experiment_accession, run_accession,
platform, library_layout, fastq_ftp, fastqbytes, submitted_ftp, submbytes).
"""
    
    url = ('http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}'
          '&result=read_run&fields='
          'study_accession,sample_accession,experiment_accession,run_accession,'
          'instrument_platform,library_layout,fastq_ftp,fastq_bytes,'
          'submitted_ftp,submitted_bytes'.format(accession))
    
    with urllib.request.urlopen(url) as lines:
        rowgenerator = (line.decode().rstrip('\n').split('\t')
                        for line in lines)
        next(rowgenerator) # iterate past header line
    
        return list(rowgenerator)


# In[19]:

def generate_entries(accession_list, log):
    """Chains together multiple entries from accessions. If one fails to load,
prints error to log instead of raising an exception."""
    
    global errors
    
    for accession in set(accession_list): # Remove duplicates
        try:
            for entry in download_entries(accession):
                yield entry
        
        except urllib.error.HTTPError:
            print('Unable to fetch data from server for accession ID {}. '
                  'Is it right?'.format(accession), file=log)
            errors += 1
            
        except urllib.error.URLError:
            print('Unable to connect to server to fetch data for study {}. '
                  'Check internet connection.'.format(accession), file=log)
            errors += 1


# In[20]:

def write_metadata(entrygenerator, log):
    """Writes the content of the entries to the log"""
    
    columnnames = ('study', 'sample', 'experiment', 'run', 'platform', 'layout',
              'ftps', 'filesizes')
    print('\t'.join(columnnames), file=log)
    
    for line in entrygenerator:
        print('\t'.join(line), file=log)


# In[21]:

def instantiate_from_entries(entries, log):
    """Given an iterable of entries, instantiates all the included ENA entries."""
    # Does not instantiate files, as they do not constitute an entry.
    # Instead, initialization of a Run instantiates its associated files.
    
    global errors
    
    studies = {}
    
    for (studyacc, sampleacc, expacc, runacc, instr, layout,
             fastqftps, fastqsizes, submittedftps, submittedsizes) in entries:
        
        try:
            if fastqftps:
                ftps, filesizes = fastqftps, fastqsizes
            elif 'fastq' in submittedftps:
                ftps, filesizes = submittedftps, submittedsizes
            else:
                raise ValueError('No fastq files found for {}'.format(runacc))
            
            # Don't use 'study = studies.get(studyacc, Study(studyacc))'
            # as Study(studyacc) will be evaluated even if studyacc in studies
            study = studies.get(studyacc)
            if study is None:
                study = Study(studyacc)
                studies[studyacc] = study
                
            sample = study.samples.get(sampleacc)
            if sample is None:
                sample = Sample(sampleacc, study)
                
            experiment = sample.experiments.get(expacc)
            if experiment is None:
                experiment = Experiment(expacc, sample, instr, layout)
            
            Run(runacc, experiment, ftps, filesizes)
                
        except Exception as error:
            print('Error when instantiating: {}'.format(error), file=log)
            errors += 1
    
    return studies


# In[22]:

def download(studies, log):
    """Downloads all the files in the given dictionary of studies in parallel"""
    
    global errors, warnings

    files = set().union(*(study.files.values() for study in studies.values()))

    if not args.quiet:
        print('Found {} files, total size: {} GB. '
              'Downloading asynchronously...'.format(len(files), 
                        round(sum(file.filesize for file in files)/1e9, 3)))
    
    fileno, filestotal = 1, len(files)
    
    def callback(message):
        nonlocal fileno, filestotal
        if not args.quiet:
            end = '\n' if fileno == filestotal else '\r'
            print('    {}/{} processed...'.format(fileno, filestotal), end=end)
        fileno += 1

    pool = Pool(processes=args.cores)
    returns = []
    for file in files:
        if not os.path.exists(file.path):
            sleep(10) # Avoid overloading servers
            
        returns.append((file, pool.apply_async(file.download, (args.letleave,),
                                                callback=callback,
                                                error_callback=callback)))

    pool.close()
    pool.join()
    
    # Process returns from child processes
    for file, returnvalue in returns:
        try:
            deinterleaved = returnvalue.get()
            
        except Exception as error:
                print('Error when downloading {}: {}'.format(
                    file.path, error), file=log)
                errors += 1
        else:
            if deinterleaved: # Deinterleaved is fw, rv pair
                file.split(*deinterleaved)
                print('File {} automatically deinterleaved.'.format(file),
                      file=log)
                
            else: # If deinterleaves is None or False
                try:
                    file.validate_download()
                except OSError as error:
                    print('Error {}'.format(error), file=log)
                    errors += 1
                finally:
                    if deinterleaved is False:
                        print('Attempted automatic deinterleaving of {}, '
                              'but it failed.'.format(file), file=log)
                        errors += 1
    return None


# In[23]:

def main(accessions):
    """Executes the main script."""
    
    global errors, warnings
    
    mkdir(args.destination)
            
    logpath = os.path.join(args.destination, 'download.log')
    logexists = os.path.exists(logpath)
    with open(logpath, 'a') as log:
        if logexists:
            print('\n\nNew run:', file=log)
    
        entries = generate_entries(accessions, log)
        
        if args.metadata:
            write_metadata(entries, log)
            if not args.quiet:
                print('Added metadata to log file')
            return
        
        studies = instantiate_from_entries(entries, log)
        
        if not args.metadata:
            for accession in studies:
                mkdir(os.path.join(args.destination, accession))
                
        download(studies, log)
    
        if args.config:
            for study in studies.values():
                path = os.path.join(args.destination, study.accession,
                                    'config.json')
                study.write_config(path, log)
            if not args.quiet:
                print('Made JSON config files. Check their values.')
    
    if not args.quiet and not (warnings or errors):
        print('Done.')
    
    elif not args.quiet:
        message = 'Completed with {} warnings and {} errors. Check log file.'
        print(message.format(warnings, errors))
        
    return None


# In[17]:

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Automatic ENA file downloader.
Given a list of ENA study accessions, does the following:

- Downloads metadata from the ENA accession IDs.
- Downloads all the associated fastq files and saves them in a directory.
- If metadata suggests paired end, but only one file is seen, deinterleaves.
- Optionally generates a configfile for use with Simon's metagenomic pipeline

Created by Jakob Nybo Nissen, jakni@dtu.dk, 2016-12-13
PLEASE CONTACT ME WITH BUGS AND REQUESTS!
""")
    
    parser.add_argument('destination', type=str,
        help='Destination directory')
    
    parser.add_argument('accessions', type=str, nargs='+',
        help='ENA accession ID(s) to download from')
    
    parser.add_argument('-c', type=int, default=1, dest='cores',
        help='No. of cores to use [1].')
    
    parser.add_argument('--quiet', action='store_true',
        help = 'Suppress non-error outputs.')
    
    parser.add_argument('--metadata', action='store_true',
        help = 'Only write metadata to log, do not download files.')
    
    parser.add_argument('--config', action='store_true',
        help = 'Make config files.')
    
    parser.add_argument('--letleave', action='store_true',
        help = 'Disable automatic deinterleaving.')

    args = parser.parse_args()
    
    if type(args.cores) != int or args.cores < 1:
        print('Zero or negative cores provided. Exiting')
    
    elif args.metadata and args.config:
        print('Cannot produce configuration files from only metadata. '
             'Choose either --metadata or --config. Exiting.')
        
    else:
        if not args.letleave:
            if os.name != 'posix':
                args.letleave = True
                print('Non-UNIX operating system found. Deinterleaving disabled.')
            
            if not os.path.isfile('/bin/bash'):
                args.letleave = True
                print('/bin/bash not detected. Deinterleaving disabled.')
    
        main(args.accessions)   


# # DEBUGGING

# In[24]:

debug = False


# In[25]:

if debug:
    class A:
        def __init__(self):
            self.destination = 'deletemedir'
            self.config = True
            self.cores = 4
            self.metadata = False
            self.quiet = False
            self.letleave = False

    args = A()


# In[27]:

if debug:
    with open('deletemelog' ,'a') as log:
        #main(['PRJEB2772'])
        print(list(generate_entries(['PRJEB2772'], log)))

