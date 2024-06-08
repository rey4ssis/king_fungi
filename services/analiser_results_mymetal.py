from Bio import SeqIO, Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
import config.celery_config as celery_config 
import os
import time
from config.celery_config import celery_app

@celery_app.task
def tranformation_results(input_dir):
    output_folder = os.path.join(input_dir,'Mymetal')
    input_file = os.path.join(output_folder, 'result_mymetal.tsv')
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")

    dados_processados = []

    with open(input_file, 'r') as f_in:
        for _ in range(8):
            next(f_in)

        for linha in f_in:
            colunas = linha.strip().split('\t')
            id_protein = colunas[0]
            valores_processados = []
            for valor in colunas[2:]:
                valor_processado = '1' if float(valor) >= 0.8 else '0'
                valores_processados.append(valor_processado)
            nova_linha = [id_protein] + valores_processados
            
            if any(nova_linha):
                dados_processados.append(nova_linha)

    arquivo_saida = os.path.join(output_folder, 'result_final.tsv')

    df = pd.DataFrame(dados_processados, columns=["ID protein", "Cálcio(Ca)", "Cobalto(Co)", "Cobre(Cu)", "Ferro(Fe)", "Pótassio(K)", "Magnésio(Mg)", "Manganês(Mn)", "Sódio(Na)", "Níquel(Ni)", "Zinco(Zn)"])

    dataframe_result =df.to_csv(arquivo_saida, sep='\t', index=False)

    return arquivo_saida

# Exemplo de uso da função
# arquivo_entrada = 'kingfungi/output/reyassis7@gmail.com/Mymetal'
# dados_processados, arquivo_saida = tranformation_results.delay(arquivo_entrada)
# print(dados_processados)
# print(arquivo_saida)
