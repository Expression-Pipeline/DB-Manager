import os
import sys
import errno
import argparse
from argparse import RawTextHelpFormatter

import time
import urllib.request
import urllib.parse
import gzip, tarfile, zipfile, re
from ftplib import FTP
from subprocess import call


download_progress = 0


def gunzipAddMaster(fileName, masterFile):
    inF = gzip.open(fileName, 'rb')

    # uncompress the gzip_path INTO THE 's' variable
    s = inF.read()
    inF.close()

    # get original filename (remove 3 characters from the end: ".gz")
    fname = fileName[:-3]

    # store uncompressed file data from 's' variable
    open(masterFile, 'ab').write(s)

def retrieveFromPathAndMerge(url, dir, sufix, date):
    urlpath = urllib.request.urlopen(url)
    string = urlpath.read().decode('utf-8')
    filelist = string.split()

    for a in filelist:
        if a.endswith(sufix):
            urllib.request.urlretrieve(url + a, filename=dir + a)

    files = os.listdir(dir)
    fasta_master_file = 'refseq-'+ date + '-protein.fasta'
    for f in files:
        gunzipAddMaster(dir+f, dir+fasta_master_file)


def retrieveFromPath(url, file, sufix):
    urlpath = urllib.request.urlopen(url)
    string = urlpath.read().decode('utf-8')
    filelist = string.split()

    for a in filelist:
        if a.endswith(sufix):
            urllib.request.urlretrieve(url + a, filename=file)



def report(block_no, block_size, file_size):
    global download_progress
    download_progress += block_size
    print("Downloaded block %i, %i/%i bytes recieved." % (block_no, download_progress, file_size))
    sys.stdout.flush()

def downloadFiles(url, fileName):
    urllib.request.urlretrieve(url, fileName, reporthook=report)
    print("File Donwload Complete: " + url)

def rsyncFiles(url,path):
    subprocess.call(["rsync", "-va", "--progress", url, " ", path])

def gunzipBig(fileName):
    try:
        call(["gunzip", fileName])
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


def main():
    # Get the date for any update, etc
    date = time.strftime("%d-%m-%Y")

    # <!-- This section contains the URL of all the databases, we can use in the future a more useful system based on text files, but for now is fine -->

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


    ## Human Databases: This databases are used for Human ##

    # Human Complete Proteome Uniprot: A UniProt complete proteome consists of the set of proteins thought to be expressed by an organism whose genome has been completely sequenced (SwissProt+TreMBL)

    humanUniprotProteomeURL = 'http://www.uniprot.org/uniprot/?query=taxonomy:9606+AND+keyword:"Complete+proteome"&force=yes&format=fasta&include=yes&compress=yes'
    humanUniprotProteomeFile = 'humanUniprotFile.fasta.gz'

    humanSwissProtURL  = 'http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=reviewed:yes+AND+taxonomy:9606&format=fasta&force=yes&include=yes'
    humanSwissProtFile = 'humanSwissFile.fasta.gz' \
                         ''
    # Human ReqSeq: The human ReqSeq contains species-specific RefSeq directories provide a cumulative set of records for transcripts and proteins for those species.

    humanRefSeqURL = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/'


    # Human ENSEMBL Peptides database: the super-set of all translations resulting from Ensembl known or novel gene predictions.

    humanENSEMBLURL  = 'ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/'
    humanENSEMBLFILE = 'human_pep.all.fa.gz'

    # <!-- End of Databases --!>

    # <!--- Conf Paths -->
    pathsep = "/"

    latest = pathsep + "latest"
    crapDir = pathsep + "crap"
    maxQuantCrap = pathsep + "mq-crap"
    path = "."
    swissprot = pathsep + "swissprot"
    humanCompleteUniprot = pathsep + "complete-human"
    humanSwissProt  = pathsep + "human-swissprot"
    humanENSEMBL    = pathsep + "human-ensembl"
    humanRefSeq     = pathsep + "human-refseq"

    # <!-- End of the Conf Paths -->

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

    args = parser.parse_args()

    if args.path:
        path = args.path

    if not os.path.isdir(path):
        print ("The path for the directory do not exists")
        parser.print_help()
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


    databaseUpdated = args.update

    # Updating the crap database


    if updateAll or ('crap' in databaseUpdated):

        print("Updating CRAP database and all the dependencies...")

        print("Downloading the CRAP database..")
        # urllib.request.urlretrieve(crapURL + crapFile, filename=path + crapDir + pathsep + crapFile)
        downloadFiles(crapURL + crapFile, path + crapDir + pathsep + crapFile)
        # rsyncFiles(crapURL + crapFile, path + crapDir + pathsep)

        print("Moving and unzip files in the current file... ")

        if not os.path.isdir(path + crapDir + pathsep + date):
            os.mkdir(path + crapDir + pathsep + date)

        os.rename(path + crapDir + pathsep + crapFile, path + crapDir + pathsep + date + pathsep + crapFile)
        symbolink(path + crapDir + pathsep + date, path + crapDir + latest)

    # Updating information for MaxQuant Contaimant database


    if updateAll or ('max-crap' in databaseUpdated):

        print("Updating MaxQuant CRAP database and all the dependencies...")

        print("Downloading the MaxQuant CRAP database..")
        downloadFiles(maxQuantcrapURL + maxQuantcrapFile,path + maxQuantCrap + pathsep + maxQuantcrapFile)

        print("Moving and unzip files in the current file... ")

        if not os.path.isdir(path + swissprot + pathsep + date):
            os.mkdir(path + maxQuantCrap + pathsep + date)

        os.rename(path + maxQuantCrap + pathsep + maxQuantcrapFile,
                  path + maxQuantCrap + pathsep + date + pathsep + maxQuantcrapFile)
        unzip(path + maxQuantCrap + pathsep + date + pathsep + maxQuantcrapFile)

        symbolink(path + maxQuantCrap + pathsep + date, path + maxQuantCrap + latest)


    # Updating uniprot

    if updateAll or ('swissprot' in databaseUpdated):

        print("Updating Swissprot/Uniprot databases and all the dependencies... ")

        print("Downloading the SwissProt database..")
        urllib.request.urlretrieve(uniprotSprotURL + uniprotSprotFile,
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

        symbolink(path + swissprot + pathsep + date, path + swissprot + latest)

    if updateAll or ('complete-human' in databaseUpdated):

        print("Updating Uniprot Complete Human database and all the dependencies... ")

        print("Downloading the Uniprot Complete Human database..")

        urllib.request.urlretrieve(humanUniprotProteomeURL, filename=path + swissprot + pathsep + humanUniprotProteomeFile)
        os.rename(path + swissprot + pathsep + humanUniprotProteomeFile, path + swissprot + pathsep + date + pathsep + humanUniprotProteomeFile)
        gunzip(path + swissprot + pathsep + date + pathsep + humanUniprotProteomeFile)

        print("Human Complete Uniprot updated!!")

    if updateAll or ('human-swissprot' in databaseUpdated):

        print("Updating Swissprot Human database and all the dependencies... ")

        print("Downloading the Swissprot Human database..")

        urllib.request.urlretrieve(humanSwissProtURL, filename=path + humanSwissProt + pathsep + humanSwissProtFile)

        if not os.path.isdir(path + humanSwissProt + pathsep + date):
            os.mkdir(path + humanSwissProt + pathsep + date)

        os.rename(path + humanSwissProt + pathsep + humanSwissProtFile, path + humanSwissProt + pathsep + date + pathsep + humanSwissProtFile)
        gunzip(path + humanSwissProt + pathsep + date + pathsep + humanSwissProtFile)

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

        symbolink(path + humanENSEMBL + pathsep + date, path + humanENSEMBL + latest)

        print("Human ENSEMBL updated!!")

    if updateAll or ('human-refseq' in databaseUpdated):

        print("Updating Human RefSeq and all the dependencies... ")

        print("Downloading the Human RefSeq database..")

        if not os.path.isdir(path + humanRefSeq + pathsep + date):
            os.mkdir(path + humanRefSeq + pathsep + date)

        sufix = "protein.faa.gz"

        retrieveFromPathAndMerge(humanRefSeqURL, path+humanRefSeq+ pathsep+date+pathsep, sufix, date)

        symbolink(path + humanRefSeq + pathsep + date, path + humanRefSeq + latest)

        print("Human ReqSeq updated!!")













if __name__ == '__main__':
    main()


