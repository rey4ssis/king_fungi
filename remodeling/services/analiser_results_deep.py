from Bio import SeqIO, Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
import os
import time
from config.celery_config import celery_app



@celery_app.task 
def extract_tm_seqs_and_filter_tsv(input_dir, output, *args, **kwargs):
    """
    Função para extrair sequências com TM do arquivo .3line, verificar esses códigos em um arquivo .tsv,
    e salvar as linhas correspondentes em um novo arquivo .tsv, excluindo a coluna 1.

    Args:
        input_dir: Diretório de entrada onde o arquivo predicted_topologies.3line está localizado.

    Returns:
        Caminho para o novo arquivo .tsv filtrado.
    """
    output_dir = os.path.join(input_dir, 'Deep')
    output_dir_my = os.path.join(input_dir, 'Mymetal')
    tsv_file = os.path.join(output_dir_my, 'input_metal.tsv')
    input_file = os.path.join(output_dir, 'predicted_topologies.3line')
    
    output_dir_my = os.path.join(input_dir, 'Mymetal')
    output_tsv = 'filtered_output.tsv'
    output_tsv_path = os.path.join(output_dir_my, output_tsv)
    
    # Coleta dos códigos de acesso com termos específicos
    codes_with_terms = set()
    with open(input_file, 'r') as f_input:
        for line in f_input:
            line = line.strip()
            if line.startswith('>'):
                line = line[1:]
                terms = line.split('|')
                if terms[-1].strip() in {'TM','SP+TM','BETA'}:
                    seq_id = terms[-2].strip()
                    codes_with_terms.add(seq_id)
    
    # Leitura do arquivo .tsv ignorando as 8 primeiras linhas
    tsv_data = pd.read_csv(tsv_file, sep='\t', skiprows=8)
    
    # Filtragem das linhas com os códigos encontrados
    filtered_tsv_data = tsv_data[tsv_data.iloc[:, 0].isin(codes_with_terms)]
    
    # Excluindo a coluna 1 (segunda coluna)
    filtered_tsv_data = filtered_tsv_data.drop(filtered_tsv_data.columns[1], axis=1)
    
    # Salvando o novo arquivo .tsv
    filtered_tsv_data.to_csv(output_tsv_path, sep='\t', index=False)
    
    # Removendo o arquivo input_metal.tsv
    os.remove(tsv_file)
    
    return output_tsv_path
