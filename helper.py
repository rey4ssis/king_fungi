import os
import re
import uuid
from Bio import SeqIO
import traceback
from services.genbank_service import download_fasta_sequence_by_id
from services.genbank_create_id import generate_fungi_by_name

def validate_email(email: str) -> bool:
    pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
    return bool(re.match(pattern, email))

def create_user_dirs(email: str) -> bool:
    base_path = os.path.abspath('.')
    input_path = os.path.join(base_path, 'input', email)
    output_path = os.path.join(base_path, 'output', email)

    os.makedirs(input_path, exist_ok=True)
    os.makedirs(output_path, exist_ok=True)

    return (input_path, output_path)

def delete_user_dirs(email: str) -> bool:
    base_path = os.path.abspath('.')
    input_path = os.path.join(base_path, 'input', email)
    output_path = os.path.join(base_path, 'output', email)

    try:
        os.rmdir(input_path)
        os.rmdir(output_path)
        return True
    except FileNotFoundError:
        # Caso o diretório não exista, retorna False
        return False
    except OSError:
        # Caso haja outro erro ao remover o diretório, retorna False
        return False
    
def generate_unique_filepath(directory: str, extension: str) -> str:
    filename = f"{uuid.uuid4()}.{extension}"
    full_path = os.path.join(directory, filename)
    
    while os.path.exists(full_path):
        filename = f"{uuid.uuid4()}.{extension}"
        full_path = os.path.join(directory, filename)
    
    return full_path


def verification_fasta(file_path):
    try:
        with open(file_path, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.description and record.seq:
                    print(record.seq)
                    return (file_path), True
    except Exception as e:
        # Trate a exceção de forma estruturada
        print(f"Erro ao verificar o arquivo {file_path}: {e}")
        return (None), 404

            
def verification_input(input_text: str):
    # Verificar se é um acession number
    if re.match(r'^[A-Za-z0-9_]+\.\d+$', input_text):
        return 0
    
    # Verificar se é um nome de fungo
    if re.match(r'^[A-Za-z\s]+$', input_text):
        return 1
    
    return 'Input inválido'

def download_sequence_file(is_sequence_ref, path_out):
    if verification_input(is_sequence_ref) == 0:  # caso for accession number
        if is_sequence_ref:
            file_down, _ = download_fasta_sequence_by_id(is_sequence_ref, path_out)
    elif verification_input(is_sequence_ref) == 1:  # caso for nome
        rt_generate, status_code = generate_fungi_by_name(is_sequence_ref)
        genbank_acession = rt_generate.get('genbank_assembly_accession')
        refseq_acession = rt_generate.get('refseq_assembly_accession')

        if refseq_acession == ' ':
            file_down, _ = download_fasta_sequence_by_id(genbank_acession, path_out)
        else:
            file_down, _ = download_fasta_sequence_by_id(refseq_acession, path_out)
    
    return file_down


def fasta_to_proteins(input_file):
    proteins = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
        
        header = None
        sequence = ''
        
        for line in lines:
            line = line.strip()
            
            if line.startswith('>'):
                if sequence:
                    proteins.append((header, sequence))
                    sequence = ''
                header = line
            else:
                sequence += line
                
        # Adicionar a última sequência
        if sequence:
            proteins.append((header, sequence))
    
    # Nome do arquivo de saída
    output_file = os.path.join(os.path.dirname(input_file), 'proteins.fasta')
    
    with open(output_file, 'w') as f:
        for header, sequence in proteins:
            f.write(f"{header}\n{sequence}\n")