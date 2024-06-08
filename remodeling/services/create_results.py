from Bio import SeqIO, Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
import config.celery_config as celery_config 
import os
import time
from config.celery_config import celery_app
from dotenv import load_dotenv
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import zipfile


# def plot_graph_all(file):
#     # Extrair o diretório do arquivo de entrada
#     output_dir = os.path.dirname(file)
    
#     df = pd.read_csv(file, sep='\t', header=0)
#     num_prot= len(df)
#     df = df.drop(columns=['ID protein'])
#     df = df.apply(pd.to_numeric, errors='coerce')
#     sums = df.sum()
    
#     # Criar cores para todas as barras
#     colors = ['#06C' if i in [0, 1, 2, 3, 4, 5, 7, 9] else '#8BC1F7' for i in range(len(sums))]
    
#     plt.figure(figsize=(10, 6))
    
#     # Ajusta o espaçamento para acomodar a legenda
#     plt.subplots_adjust(right=0.75)
    
#     bars = plt.bar(sums.index, sums.values, color=colors)
#     for bar in bars:
#         yval = bar.get_height()
#         plt.text(bar.get_x() + bar.get_width()/2, yval, round(yval, 2), va='bottom')
    
#     plt.title('Quantidade de proteínas transmenbranas que se ligam a metais')
#     plt.xlabel(f'Número de proteínas transmenbranas encontradas: {num_prot}')
#     plt.ylabel('Quantidade de proteínas')
#     plt.xticks(rotation=45)
    
#     # Define a legenda para as cores azuis
#     legendas_azuis = [plt.Rectangle((0,0),1,1, color='#06C')]
#     leg_labels_azuis = ['Micronutrientes']
    
#     # Define a legenda para as cores vermelhas
#     legendas_vermelhas = [plt.Rectangle((0,0),1,1, color='#8BC1F7')]
#     leg_labels_vermelhas = ['Elementos essenciais']
    
#     # Define o tamanho da legenda e sua posição
#     plt.legend(legendas_azuis + legendas_vermelhas, leg_labels_azuis + leg_labels_vermelhas, loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
    
#     plt.tight_layout()
#     plt_file = os.path.join(output_dir, 'grafico.png')  # Salvar na mesma pasta do arquivo de entrada
#     plt.savefig(plt_file)
#     plt.close()
#     return plt_file

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
    # #gerar gráficos
    # chart = plot_graph_all(file_path)
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
    os.remove(file_path)
    
    return output_filename


# Exemplo de uso:
# Suponha que você tenha um DataFrame chamado 'df' e o nome do arquivo seja 'meu_arquivo.tsv':
# dataframe_to_tsv_with_comments(df, 'meu_arquivo')
