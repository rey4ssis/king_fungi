from Bio import SeqIO, Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
import config.celery_config as celery_config 
import os
import time
from config.celery_config import celery_app
from dotenv import load_dotenv
import hashlib


# @celery_app.task
# def extract_headers(input_dir, arquivo_proteinas, name_file):
#     output_folder = os.path.join('input', input_dir)
#     ids_arquivo = os.path.join(output_folder, f'{name_file}.fasta')
    
#     cabecalhos = {}  # Dicionário para armazenar os cabeçalhos das proteínas

#     # Dicionário para armazenar os IDs de proteínas
#     ids_proteinas = {}

#     # Lê os IDs de proteínas do arquivo TSV e armazena em um dicionário
#     with open(ids_arquivo, 'r') as f:
#         # Pula o cabeçalho do arquivo TSV se houver
#         next(f)
#         for linha in f:
#             id_proteina = linha.strip().split('\t')[0]  # Obtém o ID da primeira coluna
#             ids_proteinas[id_proteina] = True

#     # Lê os cabeçalhos do arquivo de proteínas e armazena no dicionário de cabeçalhos
#     with open(arquivo_proteinas, 'r') as f_proteinas:
#         cabecalho = None
#         for linha_proteina in f_proteinas:
#             if linha_proteina.startswith('>'):  # Se a linha é um cabeçalho de proteína
#                 cabecalho = linha_proteina.strip()  # Armazena o cabeçalho
#                 id_proteina = cabecalho.split()[0][1:]  # Extrai o ID da proteína do cabeçalho
#                 if id_proteina in ids_proteinas:  # Verifica se o ID da proteína está nos IDs fornecidos
#                     cabecalhos[id_proteina] = cabecalho  # Armazena o cabeçalho correspondente

#     return cabecalhos, ids_arquivo

# @celery_app.task
# def update_ids(ids_arquivo, arquivo_proteinas, name_file):
#     ids_arquivo = os.path.join(f'{name_file}')
#     # Dicionário para armazenar os cabeçalhos e IDs correspondentes
    
#     cabecalhos_ids = {}

#     # Extrai os cabeçalhos do arquivo de saída e armazena no dicionário
#     cabecalhos, caminho_arquivo = extract_headers(os.path.dirname(ids_arquivo), arquivo_proteinas, os.path.basename(ids_arquivo).split('.')[0])

#     # Atualiza o arquivo de IDs com os cabeçalhos correspondentes
#     linhas_atualizadas = []
#     with open(ids_arquivo, 'r') as f_ids:
#         linhas_ids = f_ids.readlines()

#     for linha in linhas_ids:
#         id_proteina = linha.strip().split('\t')[0]
#         # Verifica se o ID existe no dicionário de cabeçalhos e IDs
#         if id_proteina in cabecalhos:
#             cabecalho_correspondente = cabecalhos[id_proteina]  # Obtém o cabeçalho correspondente ao ID
#             linha_atualizada = linha.replace(id_proteina, cabecalho_correspondente)  # Substitui o ID pelo cabeçalho
#             linhas_atualizadas.append(linha_atualizada)
#         else:
#             linhas_atualizadas.append(linha)

#     return caminho_arquivo


@celery_app.task
def create_result(file_path, organism):
    """
    Recebe um arquivo .tsv e adiciona comentários específicos no início dele.

    Parâmetros:
    file_path (str): O caminho para o arquivo .tsv existente.
    organism (str): O nome do arquivo para ser usado como base para o hash.

    Retorna:
    str: O caminho para o arquivo modificado.
    """
    # Verifica se o arquivo .tsv existe
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f'O arquivo {file_path} não foi encontrado.')

    # Comentários a serem adicionados
    comments = [
        f"Resultados da análise da genoma do organismo: {organism}.",
        "Todas as proteínas contidas nesse arquivo são proteínas validadas pelo Software DeepTMHMM.",
        "Os resultados sobre ligação com íons metálicos são criados de forma binária, onde 1 representa que a proteína possui ligação com íon metálico e 0 se não tem ligação."
    ]

    # Calcula o hash do nome do arquivo
    hash_filename = hashlib.sha256(organism.encode()).hexdigest()

    # Extrai o diretório do arquivo de entrada
    directory = os.path.dirname(file_path)

    # Cria o caminho completo para o arquivo de saída
    output_filename = os.path.join(directory, hash_filename + '.tsv')

    # Lê o conteúdo do arquivo .tsv
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Insere os comentários no início do arquivo
    with open(output_filename, 'w') as f:
        for comment in comments:
            f.write("# " + comment + "\n")
        f.writelines(lines)

    return output_filename


# Exemplo de uso:
# Suponha que você tenha um DataFrame chamado 'df' e o nome do arquivo seja 'meu_arquivo.tsv':
# dataframe_to_tsv_with_comments(df, 'meu_arquivo')
