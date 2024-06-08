from Bio import SeqIO
from mymetal import mbp, iof
import shutil
import argparse
from io import StringIO
import os
from config.celery_config import celery_app
from celery import Celery
from typing import Tuple
from billiard import current_process



def process_mymetal_filter(arquivo_tsv, arquivo_fasta, input_dir):
    tmp_dir = input_dir
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    
    nome_arquivo_fasta = os.path.basename(arquivo_fasta)
    arquivo_saida_fasta = nome_arquivo_fasta  # Use the same name for the new fasta file
    termos_guardados = set()

    with open(arquivo_tsv, 'r') as f:
        # Skip the first 8 lines of the TSV file
        for _ in range(8):
            next(f)
        
        for linha in f:
            colunas = linha.strip().split('\t')
            coluna_0 = colunas[0]
            if any(float(valor) >= 0.8 for valor in colunas[1:]):
                termos_guardados.add(coluna_0)
    
    sequencias_selecionadas = []

    with open(arquivo_fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in termos_guardados:
                sequencias_selecionadas.append(record)
    
    arquivo_novo_fasta = os.path.join(tmp_dir, arquivo_saida_fasta)
    SeqIO.write(sequencias_selecionadas, arquivo_novo_fasta, "fasta")
    
    return arquivo_novo_fasta

@celery_app.task
def run_mymetal_filter_task(arquivo_fasta, input_dir, output_dir):
    # Ensure output directory is `output/Deep`
    output_dir = os.path.join(output_dir, 'Mymetal')
    os.makedirs(output_dir, exist_ok=True)
            
    # Ensure the current process is not daemonized
    current_process()._config['daemon'] = False
        
    caminho_arquivo_fasta = os.path.join(input_dir, arquivo_fasta)
        
    if not os.path.exists(caminho_arquivo_fasta):
        raise FileNotFoundError(f'The fasta file {arquivo_fasta} was not found in {input_dir}')
        
    # Run the mymetal model on the fasta file
    mbp_preb = mbp.predict(caminho_arquivo_fasta)
    # Save the model output as a CSV file in the `output/Deep` directory
    arquivo_tsv = os.path.join(output_dir, 'input_metal.tsv')

    print(iof.save_out_csv(mbp_preb, arquivo_tsv))
        
    # Process the generated CSV file
    arquivo_saida_fasta = process_mymetal_filter(arquivo_tsv, caminho_arquivo_fasta, input_dir)
        
    return arquivo_saida_fasta 


def main():
    parser = argparse.ArgumentParser(description='Process a FASTA file.')
    parser.add_argument('--input-file', required=True, help='Path to the input FASTA file')
    parser.add_argument('--input-dir', required=True, help='Directory containing the input FASTA file')
    parser.add_argument('--output-dir', required=True, help='Directory to save the output files')
    
    args = parser.parse_args()
    
    input_file = args.input_file
    input_dir = args.input_dir
    output_dir = args.output_dir

    result_file, _ = run_mymetal_filter_task(input_file, input_dir=input_dir, output_dir=output_dir)
    # print(f'Result file saved at: {result_file}')


# file_fasta='Novo.fasta'
# input_dir='/home/rey/Documentos/KF/kingfungi/input/reyassis7@gmail.com'
# output_dir='/home/rey/Documentos/KF/kingfungi/output/reyassis7@gmail.com/Deep'
# result_fasta, tsv= run_mymetal_filter_task(file_fasta, input_dir, output_dir)