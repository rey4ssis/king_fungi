# Execs MyMetal: Identifies the sequence nutrients
# input: read a fasta sequence file from input dir
# output: save a csv file in output dir

from Bio import SeqIO
import os
import uuid
import argparse
import os
import shutil
import zipfile
import biolib
from config.celery_config import celery_app
from billiard import current_process



def check_input_file(fasta_file: str, input_dir: str) -> str:
    fasta_path = os.path.join(input_dir, fasta_file)
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f'O arquivo {fasta_file} não foi encontrado em {input_dir}')
    if not fasta_file.endswith('.fasta'):
        raise ValueError(f'O arquivo {fasta_file} não possui a extensão .fasta')

    return fasta_path

@celery_app.task
def processDeep(fasta_file, input_dir, output_dir) -> str:
    # Clear the output directory
    output_dir = output_dir
    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(output_dir):
        file_path = os.path.join(output_dir, file)
        # if os.path.isfile(file_path):
        #     # os.remove(file_path)
    fasta_path = check_input_file(fasta_file, input_dir)
    deeptmhmm = biolib.load('DTU/DeepTMHMM')
    deeptmhmm_job = deeptmhmm.cli(args=f'--fasta {fasta_path}', machine='local', timeout=10000000000000000000000)
    deeptmhmm_job.save_files(output_dir)
    # print(deeptmhmm_job.list_output_files())
    # file_3line = deeptmhmm_job.get_output_file('/predicted_topologies.3line')

def main():
    parser = argparse.ArgumentParser(description='Processa arquivo FASTA.')
    parser.add_argument('--input-file', required=True, help='Caminho para o arquivo FASTA de entrada')
    args = parser.parse_args()
   
    input_file = args.input_file
    output_dir= f'{output_dir}/Deep'

    # Process the input file using processMetal method
 
    result_file = processDeep(input_file, input_dir='input', output_dir=output_dir)
    print(f'Result file saved at: {result_file}')