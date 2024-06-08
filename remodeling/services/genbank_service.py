import gzip
import os
import re
import shutil
import urllib.request
from ftplib import FTP
from urllib.parse import urlparse
from config.celery_config import celery_app
from Bio import Entrez, SeqIO
from dotenv import load_dotenv

load_dotenv(override=True)

Entrez.email = os.getenv('EMAIL_ENTREZ')


def search_fungi_id_by_name(especie: str):
    """
    The function `search_fungi_id_by_name` searches for the assembly IDs of fungi
    species based on their names.

    :param especie: The parameter "especie" is a string that represents the name of
    a species of fungi
    :type especie: str
    :return: The function `search_fungi_id_by_name` returns two values:
    `genbank_assembly_accession` and `refseq_assembly_accession`.
    """
    # Configurando o email para evitar problemas com o NCBI
    Entrez.email = "seu@email.com"

    # Construindo o termo de pesquisa
    search_term = f'"{especie}"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter])'

    # Fazendo a pesquisa no NCBI
    handle = Entrez.esearch(db='assembly', term=search_term, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    count = int(record['Count'])

    if count > 0:
        assembly_id = record['IdList'][0]

        # Obtendo os detalhes do ID de montagem
        handle = Entrez.esummary(db='assembly', id=assembly_id)
        record = Entrez.read(handle)
        handle.close()

        # Inicializando os acessos ao GenBank e RefSeq
        genbank_assembly_accession = None
        refseq_assembly_accession = None

        # Iterando sobre os documentos
        for doc in record['DocumentSummarySet']['DocumentSummary']:
            for synonym in doc['Synonym']:
                if 'Genbank' in synonym:
                    genbank_assembly_accession = doc['AssemblyAccession']
                if 'RefSeq' in synonym:
                    refseq_assembly_accession = doc['AssemblyAccession']

        return genbank_assembly_accession, refseq_assembly_accession

    else:
        return None, None



def check_ftp_file_existence(ftp_url: str, file_name: str):
    """
    The function `check_ftp_file_existence` checks if a file exists on an FTP
    server given the FTP URL and file name.

    :param ftp_url: The `ftp_url` parameter is a string that represents the URL of
    the FTP server. It should follow the format
    `ftp://<server_address>/<file_path>`, where `<server_address>` is the address
    of the FTP server and `<file_path>` is the path to the directory where the file
    :param file_name: The name of the file you want to check for existence on the
    FTP server
    :return: a boolean value indicating whether the specified file exists on the
    FTP server.
    """
    url_parts = urlparse(ftp_url)
    server_address = url_parts.netloc
    file_path = url_parts.path
    try:
        ftp = FTP(server_address)
        ftp.login()

        ftp.cwd(file_path)
        file_list = ftp.nlst()
        file_exists = file_name in file_list

        ftp.quit()

        return file_exists

    except Exception as e:
        print(f'Error checking FTP file existence: {str(e)}')
        return False

@celery_app.task
def download_fasta_sequence_by_id(acession_number: str, file_path):
    """
    The function `download_fasta_sequence_by_id` downloads a FASTA sequence file
    from an FTP server using an accession number and returns the file path,
    organism name, and any errors encountered.

    :param acessiore
    """
    file_path = f'{file_path}'
    print(file_path)
    try:
        temp_genomes_path = os.path.join(
            os.getcwd(), f'input/{file_path}'
        )

        if not os.path.exists(temp_genomes_path):
            os.makedirs(temp_genomes_path)

        handle = Entrez.esearch(
            db='assembly',
            term=f'{acession_number}[Assembly Accession]',
            retmax=1,
        )

        record = Entrez.read(handle)
        handle.close()
        if not record['IdList']:
            raise ValueError(
                'Assembly not found for the given accession number.'
            )

        assembly_id = record['IdList'][0]

        handle = Entrez.esummary(db='assembly', id=assembly_id)
        record = Entrez.read(handle)
        handle.close()

        full_organism_name = record['DocumentSummarySet']['DocumentSummary'][
            0
        ]['SpeciesName']

        ftp_url = record['DocumentSummarySet']['DocumentSummary'][0].get(
            'FtpPath_RefSeq'
        )
        if not ftp_url:
            ftp_url = record['DocumentSummarySet']['DocumentSummary'][0].get(
                'FtpPath_GenBank'
            )
            if not ftp_url:
                raise ValueError(
                    'No FTP path found for the given accession number.'
                )

        fasta_file = ftp_url.split('/')[-1] + '_protein.faa.gz'
        fasta_url = ftp_url + '/' + fasta_file

        if not check_ftp_file_existence(ftp_url=ftp_url, file_name=fasta_file):
            raise ValueError(f'File {fasta_file} not found on the FTP server.')

        fasta_folder_path = os.path.join(
            temp_genomes_path, f'{acession_number}.fasta.gz'
        )

        if not os.path.exists(fasta_folder_path):
            urllib.request.urlretrieve(fasta_url, fasta_folder_path)

        fasta_output_path = os.path.join(
            temp_genomes_path, f'{acession_number}.fasta'
        )

        with gzip.open(fasta_folder_path, 'rb') as f_in:
            with open(fasta_output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        
        return f'{acession_number}.fasta', False

    except Exception as error:
        return False, False, str(error)


def get_protein_name(genbank_id: str):
    """
    The function `get_protein_name` retrieves the protein name from a GenBank ID.

    :param genbank_id: The `genbank_id` parameter is a string that represents the
    unique identifier of a protein in the GenBank database
    :type genbank_id: str
    :return: the protein name associated with the given GenBank ID.
    """

    handle = Entrez.efetch(
        db='protein', id=genbank_id, rettype='gb', retmode='text'
    )
    record = SeqIO.read(handle, 'genbank')
    handle.close()
    protein_name = record.description
    return protein_name



# pasta= 'reyassis@gmail.com'
# download_fasta_sequence_by_id('GCA_029542845.1', pasta)