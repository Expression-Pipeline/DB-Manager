import os
import argparse
from argparse import RawTextHelpFormatter

import time
import urllib.request
import gzip
import tarfile
import zipfile

def symbolink(src, desc):
    files = [f for f in os.listdir(src) if os.path.isfile(f)]
    for f in files:
        os.symlink(f, desc)

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

    crapURL  = "ftp://ftp.thegpm.org/fasta/cRAP/"
    crapFile = "crap.fasta"

    maxQuantcrapURL = "http://maxquant.org/"
    maxQuantcrapFile = "contaminants.zip"

    # General Uniprot-SwissProt databases

    uniprotSprotURL  = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/"
    uniprotSprotFile = "uniprot_sprot.fasta.gz"

    uniprotSprotDatURL  = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/"
    uniprotSprotDatFile = "uniprot_sprot.dat.gz"

    uniprotSprotTaxonomiesDumpURL  = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
    uniprotSprotTaxonomiesDumpFile = "taxdump.tar.gz"

    uniprotSprotSpecListURL  = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/docs/"
    uniprotSprotSpecListFile = "speclist.txt"


    ## Human Databases: This databases are used for Human ##

    # Human Complete Proteome Uniprot: A UniProt complete proteome consists of the set of proteins thought to be expressed by an organism whose genome has been completely sequenced (SwissProt+TreMBL)

    humanUniprotProteomeURL = 'http://www.uniprot.org/uniprot/?query=taxonomy:9606+AND+keyword:"Complete+proteome"&force=yes&format=fasta&include=yes'
    humanSwissProtURL       = 'http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=&fil=reviewed:yes+AND+taxonomy:9606&format=fasta&force=yes&include=yes'

    # Human ReqSeq: The human ReqSeq contains species-specific RefSeq directories provide a cumulative set of records for transcripts and proteins for those species.

    humanReqSeqRootURL = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/'

    # Human ENSEMBL Peptides database: the super-set of all translations resulting from Ensembl known or novel gene predictions.

    humanENSEMBLURL    = 'ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/'


    #  <!-- End of Databases --!>

    # <!--- Conf Paths -->
    pathsep = "/"

    latest = pathsep + "latest"
    crapDir = pathsep + "crap"
    maxQuantCrap = pathsep + "mq-crap"
    path="."
    swissprot= pathsep + "swissprot"

    # <!-- End of the Conf Paths -->

    updateAll = False

    parser = argparse.ArgumentParser(description='This python script update databases from different providers including contaminants databases, uniprot, etc:\n'
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

    parser.add_argument('-s', '--start', action='store_true', help='This parameter will go thorough all the databases and download an create the decoy versions')
    parser.add_argument('-p', '--path' , help='This variable control the place were the user want to generate the database structure, if not information is provided the structure will be generated in .')
    parser.add_argument('-u', '--update', nargs='+', help='This variable enables the update of a particular database including the the generation of decoys, (supported databases: crap, swissprot, human-uniprot)')



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

    databaseUpdated = args.update

    # Updating the crap database


    if updateAll or ('crap' in databaseUpdated):

        print("Updating CRAP database and all the dependencies...")

        print("Downloading the CRAP database..")
        urllib.request.urlretrieve(crapURL + crapFile, filename=path + crapDir + pathsep + crapFile)

        print("Moving and unzip files in the current file... ")

        os.rename(path + crapDir + pathsep + crapFile, path + crapDir + pathsep + date + pathsep + crapFile)
        symbolink(path + crapDir + pathsep + date, f, path + crapDir + latest)

    # Updating information for MaxQuant Contaimant database


    if updateAll or ('max-crap' in databaseUpdated):

        print("Updating MaxQuant CRAP database and all the dependencies...")

        print("Downloading the MaxQuant CRAP database..")
        urllib.request.urlretrieve(maxQuantcrapURL + maxQuantcrapFile, filename=path + maxQuantCrap + pathsep + maxQuantcrapFile)

        print("Moving and unzip files in the current file... ")

        os.rename(path + maxQuantCrap + pathsep + maxQuantcrapFile, path + maxQuantCrap + pathsep + date + pathsep + maxQuantcrapFile)
        unzip(path + maxQuantCrap + pathsep + date + pathsep + maxQuantcrapFile)

        symbolink(path + maxQuantCrap + pathsep + date, f, path + maxQuantCrap + latest)


    # Updating uniprot

    if updateAll or ('swissprot' in databaseUpdated):

        print("Updating Uniprot databases and all the dependencies... ")

        print("Downloading the SwissProt database..")
        urllib.request.urlretrieve(uniprotSprotURL + uniprotSprotFile, filename=path + swissprot + pathsep + uniprotSprotFile)

        print("Downloading the SwissProt Dat...")
        urllib.request.urlretrieve(uniprotSprotDatURL + uniprotSprotDatFile, filename=path + swissprot + pathsep + uniprotSprotDatFile)

        print("Downloading the SwissProt Taxonomy Dump... ")
        urllib.request.urlretrieve(uniprotSprotTaxonomiesDumpURL + uniprotSprotTaxonomiesDumpFile, filename=path + swissprot + pathsep + uniprotSprotTaxonomiesDumpFile)

        print("Downloading the SwissProt species list ...")
        urllib.request.urlretrieve(uniprotSprotSpecListURL + uniprotSprotSpecListFile, filename=path + swissprot + pathsep + uniprotSprotSpecListFile)

        print("Creating the release folder and moving the files to it....")

        if not os.path.isdir(path + swissprot + pathsep + date):
            os.mkdir(path + swissprot + pathsep + date)

        print("Moving and unzip files in the current file... ")

        os.rename(path + swissprot + pathsep + uniprotSprotFile, path + swissprot + pathsep + date + pathsep+ uniprotSprotFile)
        gunzip(path + swissprot + pathsep + date + pathsep+ uniprotSprotFile)

        os.rename(path + swissprot + pathsep + uniprotSprotDatFile, path + swissprot + pathsep + date + pathsep + uniprotSprotDatFile)
        gunzip(path + swissprot + pathsep + date + pathsep+ uniprotSprotDatFile)

        os.rename(path + swissprot + pathsep + uniprotSprotTaxonomiesDumpFile, path + swissprot + pathsep + date + pathsep + uniprotSprotTaxonomiesDumpFile)
        untar(path + swissprot + pathsep + date + pathsep + uniprotSprotTaxonomiesDumpFile)

        os.rename(path + swissprot + pathsep + uniprotSprotSpecListFile, path + swissprot + pathsep + date + pathsep + uniprotSprotSpecListFile)

        print("Creating the links for latest version to the current folder... ")

        symbolink(path + swissprot + pathsep + date, f, path + swissprot + latest)

if __name__ == '__main__':
    main()


