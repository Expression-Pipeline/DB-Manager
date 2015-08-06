import os
import sys
import errno
import argparse
from argparse import RawTextHelpFormatter


import time
import urllib
import urllib
import gzip, tarfile, zipfile, re
from ftplib import FTP
from subprocess import call

# Global variables

download_progress = 0
pathsep = "/"
date = time.strftime("%d-%m-%Y")


# The crap databases

crapURL = "ftp://ftp.thegpm.org/fasta/cRAP/"
crapFile = "crap.fasta"

maxQuantcrapURL = "http://maxquant.org/"
maxQuantcrapFile = "contaminants.zip"

# General Uniprot-SwissProt databases

uniprotSprotURL = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/"
uniprotSprotFile = "uniprot_sprot.fasta.gz"

uniprotSprotDatURL = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/"
uniprotSprotDatFile = "uniprot_sprot.dat.gz"

uniprotSprotTaxonomiesDumpURL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
uniprotSprotTaxonomiesDumpFile = "taxdump.tar.gz"

uniprotSprotSpecListURL = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/docs/"
uniprotSprotSpecListFile = "speclist.txt"

humanUniprotProteomeURL = 'http://www.uniprot.org/uniprot/?query=taxonomy:9606+AND+keyword:"Complete+proteome"&force=yes&format=fasta&include=yes&compress=yes'
humanUniprotProteomeFile = 'humanUniprotFile.fasta.gz'

humanSwissProtURL  = 'http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=reviewed:yes+AND+taxonomy:9606&format=fasta&force=yes&include=yes'
humanSwissProtFile = 'humanSwissFile.fasta.gz'

# Human ReqSeq: The human ReqSeq contains species-specific RefSeq directories provide a cumulative set of records for transcripts and proteins for those species.

humanRefSeqURL = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/'


# Human ENSEMBL Peptides database: the super-set of all translations resulting from Ensembl known or novel gene predictions.

humanENSEMBLURL  = 'ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/'
humanENSEMBLFILE = 'human_pep.all.fa.gz'

# Spectral Libraries

gpmdbURL = "ftp://ftp.thegpm.org/projects/xhunter/libs/"
nistURL = "ftp://chemdata.nist.gov/download/peptide_library/libraries/"



# Databases path
latest = pathsep + "latest"
crapDir = pathsep + "crap"
maxQuantCrap = pathsep + "mq-crap"
path = "."
swissprot = pathsep + "swissprot"
humanCompleteUniprot = pathsep + "human-complete"
humanSwissProt  = pathsep + "human-swissprot"
humanENSEMBL    = pathsep + "human-ensembl"
humanRefSeq     = pathsep + "human-refseq"
gpmdbSpectraLib = pathsep + "gpmdb-splib"
nistSpectraLib  = pathsep + "nist-splib"


def runDecoy(dir, decoyCommand):
    files = os.listdir(dir)
    print("Running the Decoy commandline" + decoyCommand)
    for f in files:
        if f.endswith(".fasta") and not f.endswith("_decoy.fasta"):
            finalCommand = decoyCommand%(dir + pathsep + f,dir + pathsep + f[:-6] + "_decoy.fasta")
            print("Running the Decoy commandline: " + finalCommand)
            try:
                call(finalCommand, shell=True)
            except:
                print ("Unexpected error from systems:", sys.exc_info())



def gunzipAddMaster(fileName, masterFile):
    inF = gzip.open(fileName, 'rb')

    # uncompress the gzip_path INTO THE 's' variable
    s = inF.read()
    inF.close()

    # store uncompressed file data from 's' variable
    open(masterFile, 'ab').write(s)

def retrieveFromPathAndMerge(url, dir, sufix, date):
    urlpath = urllib.urlopen(url)
    string = urlpath.read().decode('utf-8')
    filelist = string.split()

    for a in filelist:
        if a.endswith(sufix):
            urllib.urlretrieve(url + a, filename=dir + a)

    fasta_master_file = 'refseq-'+ date + '-protein.fasta'

    if os.path.isfile(dir + fasta_master_file):
        os.remove(dir + fasta_master_file)

    files = os.listdir(dir)


    for f in files:
        if f.endswith(".gz"):
            gunzipAddMaster(dir+f, dir+fasta_master_file)


def retrieveFromPath(url, file, sufix):
    urlpath = urllib.urlopen(url)
    string = urlpath.read().decode('utf-8')
    filelist = string.split()

    for a in filelist:
        if a.endswith(sufix):
            urllib.urlretrieve(url + a, filename=file)


def report(block_no, block_size, file_size):
    global download_progress
    download_progress += block_size
    print("Downloaded block %i, %i/%i bytes recieved." % (block_no, download_progress, file_size))
    sys.stdout.flush()

def downloadFiles(url, fileName):
    urllib.urlretrieve(url, fileName, reporthook=report)
    print("File Donwload Complete: " + url)

def rsyncFiles(url,path):
    subprocess.call("wget " + url + "--recursive -np", shell=True)

def gunzipBig(fileName):
    try:
        call(["gunzip -f ", fileName])
    except:
        print ("Unexpected error from systems:", sys.exc_info())

def symbolink(src, desc):
    files = os.listdir(src)
    for f in files:
        fileSource = src + "/" + f
        fileDesc = desc + "/" + f
        if os.path.isfile(fileSource):
            try:
                os.symlink(fileSource, fileDesc)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    os.remove(fileDesc)
                    os.symlink(fileSource, fileDesc)


def gunzip(fileName):
    inF = gzip.open(fileName, 'rb')

    # uncompress the gzip_path INTO THE 's' variable
    s = inF.read()
    inF.close()

    # get original filename (remove 3 characters from the end: ".gz")
    fname = fileName[:-3]

    # store uncompressed file data from 's' variable
    open(fname, 'wb').write(s)


def untar(fileName):
    tfile = tarfile.open(fileName, 'r:gz')

    rootPath = os.path.dirname(fileName)

    tfile.extractall(rootPath)


def unzip(fileName):
    zfile = zipfile.ZipFile(fileName)
    for name in zfile.namelist():
        (dirname, filename) = os.path.split(name)
        print ("Decompressing " + filename + " on " + dirname)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        zfile.extract(name, dirname)

def downloadWgetFiles(url, dir):
    try:
        call("wget " + url + " --recursive -np -P " + dir, shell=True)
    except:
        print ("Unexpected error from systems:", sys.exc_info())

def cleanDirectories(url , dir):
    url = url.replace("ftp://", "")
    url = url.replace("http://", "")
    os.rename(dir + url, dir + pathsep)


def main():
    updateAll = False

    parser = argparse.ArgumentParser(
        description='This python script update databases from different providers including contaminants databases, uniprot, etc:\n'
                    '\n\t The structure of the data is in this way:\t'
                    '\n\t\t- crap'
                    '\n\t\t\t- version-date'
                    '\n\t\t\t\t- crap'
                    '\n\t\t\t\t- crap+decoy'
                    '\n\t\t\t- current (a link to the latest version)'
                    '\n\t\t- swissprot'
                    '\n\t\t\t- version-date'
                    '\n\t\t\t\t- uniprotKB'
                    '\n\t\t\t\t- uniprotKB+decoy'
                    '\n\t\t\t- current (a link to the latest version)'
                    '\n\t\t...', formatter_class=RawTextHelpFormatter)

    parser.add_argument('-s', '--start', action='store_true',
                        help='This parameter will go thorough all the databases and download an create the decoy versions')
    parser.add_argument('-p', '--path',
                        help='This variable control the place were the user want to generate the database structure, if not information is provided the structure will be generated in .')
    parser.add_argument('-u', '--update', nargs='+',
                        help='This variable enables the update of a particular database including the the generation of decoys, (supported databases: crap, swissprot, human-uniprot)')
    parser.add_argument('-d', '--decoy', help="The command line option to generate the decoy, using different tools (eg: 'Decoy -shuffle -in s -out s')")

    args = parser.parse_args()

    if args.decoy is None and args.path is None and args.start is False and args.update is None:
        parser.print_help()

    if args.path:
        path = args.path

    if not os.path.isdir(path):
        print ("The path for the directory do not exists")
        print parser.print_help()
        exit(1)

    if args.start:
        print('The database structure will be generated in the following path: ' + os.path.dirname(path))
        updateAll = True

    if not os.path.isdir(path + crapDir):
        os.mkdir(path + crapDir)
        if not os.path.isdir(path + crapDir + latest):
            os.mkdir(path + crapDir + latest)

    if not os.path.isdir(path + maxQuantCrap):
        os.mkdir(path + maxQuantCrap)
        if not os.path.isdir(path + maxQuantCrap + latest):
            os.mkdir(path + maxQuantCrap + latest)

    if not os.path.isdir(path + swissprot):
        os.mkdir(path + swissprot)
        if not os.path.isdir(path + swissprot + latest):
            os.mkdir(path + swissprot + latest)

    if not os.path.isdir(path + humanSwissProt):
        os.mkdir(path + humanSwissProt)
        if not os.path.isdir(path + humanSwissProt + latest):
            os.mkdir(path + humanSwissProt + latest)

    if not os.path.isdir(path + humanENSEMBL):
        os.mkdir(path + humanENSEMBL)
        if not os.path.isdir(path + humanENSEMBL + latest):
            os.mkdir(path + humanENSEMBL + latest)

    if not os.path.isdir(path + humanRefSeq):
        os.mkdir(path + humanRefSeq)
        if not os.path.isdir(path + humanRefSeq + latest):
            os.mkdir(path + humanRefSeq + latest)

    if not os.path.isdir(path + humanCompleteUniprot):
        os.mkdir(path + humanCompleteUniprot)
        if not os.path.isdir(path + humanCompleteUniprot + latest):
            os.mkdir(path + humanCompleteUniprot + latest)

    if not os.path.isdir(path + gpmdbSpectraLib):
        os.mkdir(path + gpmdbSpectraLib)
        if not os.path.isdir(path + gpmdbSpectraLib + latest):
            os.mkdir(path + gpmdbSpectraLib + latest)

    if not os.path.isdir(path + nistSpectraLib):
        os.mkdir(path + nistSpectraLib)
        if not os.path.isdir(path + nistSpectraLib + latest):
            os.mkdir(path + nistSpectraLib + latest)


    databaseUpdated = args.update
    if databaseUpdated is None:
        databaseUpdated = []

    if updateAll or ('crap' in databaseUpdated):

        print("Updating CRAP database and all the dependencies...")

        print("Downloading the CRAP database..")
        # urllib.request.urlretrieve(crapURL + crapFile, filename=path + crapDir + pathsep + crapFile)
        try:
            downloadFiles(crapURL + crapFile, path + crapDir + pathsep + crapFile)
            print("Moving and unzip files in the current file... ")
            if not os.path.isdir(path + crapDir + pathsep + date):
                os.mkdir(path + crapDir + pathsep + date)
            os.rename(path + crapDir + pathsep + crapFile, path + crapDir + pathsep + date + pathsep + crapFile)

            if args.decoy:
                runDecoy(path + crapDir + pathsep + date, args.decoy)

            symbolink(path + crapDir + pathsep + date, path + crapDir + latest)

        except:
            print("The download of CRAP database wasn't successful...")

    # Updating information for MaxQuant Contaimant database


    if updateAll or ('max-crap' in databaseUpdated):

        print("Updating MaxQuant CRAP database and all the dependencies...")

        try:
            print("Downloading the MaxQuant CRAP database..")
            downloadFiles(maxQuantcrapURL + maxQuantcrapFile,path + maxQuantCrap + pathsep + maxQuantcrapFile)

            print("Moving and unzip files in the current file... ")

            if not os.path.isdir(path + swissprot + pathsep + date):
                os.mkdir(path + maxQuantCrap + pathsep + date)

            os.rename(path + maxQuantCrap + pathsep + maxQuantcrapFile,
                  path + maxQuantCrap + pathsep + date + pathsep + maxQuantcrapFile)
            unzip(path + maxQuantCrap + pathsep + date + pathsep + maxQuantcrapFile)

            if args.decoy:
                runDecoy(path + maxQuantCrap + pathsep + date, args.decoy)

            symbolink(path + maxQuantCrap + pathsep + date, path + maxQuantCrap + latest)

        except:
            print("The download of CRAP database wasn't successful...")

    # Updating Uniprot

    if updateAll or ('swissprot' in databaseUpdated):

        print("Updating Swissprot/Uniprot databases and all the dependencies... ")

        print("Downloading the SwissProt database..")
        urllib.urlretrieve(uniprotSprotURL + uniprotSprotFile,
                                   filename=path + swissprot + pathsep + uniprotSprotFile)

        print("Downloading the SwissProt Dat...")
        downloadFiles(uniprotSprotDatURL + uniprotSprotDatFile, path + swissprot + pathsep + uniprotSprotDatFile)
        #
        print("Downloading the SwissProt Taxonomy Dump... ")
        downloadFiles(uniprotSprotTaxonomiesDumpURL + uniprotSprotTaxonomiesDumpFile,path + swissprot + pathsep + uniprotSprotTaxonomiesDumpFile)

        print("Downloading the SwissProt species list ...")
        downloadFiles(uniprotSprotSpecListURL + uniprotSprotSpecListFile,path + swissprot + pathsep + uniprotSprotSpecListFile)

        print("Creating the release folder and moving the files to it....")

        if not os.path.isdir(path + swissprot + pathsep + date):
            os.mkdir(path + swissprot + pathsep + date)

        print("Moving and unzip files in the current file... ")

        os.rename(path + swissprot + pathsep + uniprotSprotFile,path + swissprot + pathsep + date + pathsep + uniprotSprotFile)
        gunzip(path + swissprot + pathsep + date + pathsep + uniprotSprotFile)


        os.rename(path + swissprot + pathsep + uniprotSprotDatFile, path + swissprot + pathsep + date + pathsep + uniprotSprotDatFile)
        gunzipBig(path + swissprot + pathsep + date + pathsep + uniprotSprotDatFile)

        os.rename(path + swissprot + pathsep + uniprotSprotTaxonomiesDumpFile, path + swissprot + pathsep + date + pathsep + uniprotSprotTaxonomiesDumpFile)
        untar(path + swissprot + pathsep + date + pathsep + uniprotSprotTaxonomiesDumpFile)

        os.rename(path + swissprot + pathsep + uniprotSprotSpecListFile,
                  path + swissprot + pathsep + date + pathsep + uniprotSprotSpecListFile)

        print("Creating the links for latest version to the current folder... ")

        if args.decoy:
            runDecoy(path + swissprot + pathsep + date, args.decoy)

        symbolink(path + swissprot + pathsep + date, path + swissprot + latest)

    if updateAll or ('complete-human' in databaseUpdated):

        print("Updating Uniprot Complete Human database and all the dependencies... ")

        print("Downloading the Uniprot Complete Human database..")

        if not os.path.isdir(path + humanCompleteUniprot + pathsep + date):
            os.mkdir(path + humanCompleteUniprot + pathsep + date)

        urllib.urlretrieve(humanUniprotProteomeURL, filename=path + humanCompleteUniprot + pathsep + humanUniprotProteomeFile)
        os.rename(path + humanCompleteUniprot + pathsep + humanUniprotProteomeFile, path + humanCompleteUniprot + pathsep + date + pathsep + humanUniprotProteomeFile)

        gunzip(path + humanCompleteUniprot + pathsep + date + pathsep + humanUniprotProteomeFile)

        if args.decoy:
            runDecoy(path + humanCompleteUniprot + pathsep + date, args.decoy)

        symbolink(path + humanCompleteUniprot + pathsep + date, path + humanCompleteUniprot + latest)

        print("Human Complete Uniprot updated!!")

    if updateAll or ('human-swissprot' in databaseUpdated):

        print("Updating Swissprot Human database and all the dependencies... ")

        print("Downloading the Swissprot Human database..")

        urllib.urlretrieve(humanSwissProtURL, filename=path + humanSwissProt + pathsep + humanSwissProtFile)

        if not os.path.isdir(path + humanSwissProt + pathsep + date):
            os.mkdir(path + humanSwissProt + pathsep + date)

        os.rename(path + humanSwissProt + pathsep + humanSwissProtFile, path + humanSwissProt + pathsep + date + pathsep + humanSwissProtFile)
        gunzip(path + humanSwissProt + pathsep + date + pathsep + humanSwissProtFile)

        if args.decoy:
            runDecoy(path + humanSwissProt + pathsep + date, args.decoy)

        symbolink(path + humanSwissProt + pathsep + date, path + humanSwissProt + latest)

        print("Human Complete SwissProt updated!!")

    if updateAll or ('human-ensembl' in databaseUpdated):

        print("Updating human-ensembl and all the dependencies... ")

        print("Downloading the Human ENSEMBL database..")

        sufix = '.pep.all.fa.gz'

        retrieveFromPath(humanENSEMBLURL, path + humanENSEMBL + pathsep + humanENSEMBLFILE, sufix)

        if not os.path.isdir(path + humanENSEMBL + pathsep + date):
            os.mkdir(path + humanENSEMBL + pathsep + date)

        os.rename(path + humanENSEMBL + pathsep + humanENSEMBLFILE, path + humanENSEMBL + pathsep + date + pathsep + humanENSEMBLFILE)
        gunzip(path + humanENSEMBL + pathsep + date + pathsep + humanENSEMBLFILE)

        if args.decoy:
            runDecoy(path + humanENSEMBL + pathsep + date, args.decoy)

        symbolink(path + humanENSEMBL + pathsep + date, path + humanENSEMBL + latest)

        print("Human ENSEMBL updated!!")

    if updateAll or ('human-refseq' in databaseUpdated):

        print("Updating Human RefSeq and all the dependencies... ")

        print("Downloading the Human RefSeq database..")

        if not os.path.isdir(path + humanRefSeq + pathsep + date):
            os.mkdir(path + humanRefSeq + pathsep + date)

        sufix = "protein.faa.gz"

        retrieveFromPathAndMerge(humanRefSeqURL, path+humanRefSeq+ pathsep+date+pathsep, sufix, date)

        if args.decoy:
            runDecoy(path + humanRefSeq + pathsep + date, args.decoy)

        symbolink(path + humanRefSeq + pathsep + date, path + humanRefSeq + latest)

        print("Human ReqSeq updated!!")

    if updateAll or ('gpmdb-splib' in databaseUpdated):

        print("Updating the gpmdb spectrum library... ")

        if not os.path.isdir(path + gpmdbSpectraLib + pathsep + date):
            os.mkdir(path + gpmdbSpectraLib + pathsep + date)

        downloadWgetFiles(gpmdbURL, path+gpmdbSpectraLib +pathsep+date+pathsep)

        cleanDirectories(gpmdbURL, path + gpmdbSpectraLib + pathsep + date + pathsep)

        symbolink(path + gpmdbSpectraLib + pathsep + date, path + gpmdbSpectraLib + latest)

    if updateAll or ('nist-splib' in databaseUpdated):

        print("Updating the NIST spectrum library..")

        if not os.path.isdir(path + nistSpectraLib + pathsep + date):
            os.mkdir(path + nistSpectraLib + pathsep + date)

        downloadWgetFiles(nistURL, path+nistSpectraLib + pathsep + date + pathsep)

        cleanDirectories(nistURL, path+nistSpectraLib + pathsep + date + pathsep)

        symbolink(path + nistSpectraLib + pathsep + date, path + nistSpectraLib + latest)





if __name__ == '__main__':
    main()


