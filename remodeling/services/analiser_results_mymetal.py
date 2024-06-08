from Bio import SeqIO, Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
import os
import time
from config.celery_config import celery_app

@celery_app.task
def tranformation_results(input_dir):
    output_folder = os.path.join(input_dir, 'Mymetal')
    input_file = os.path.join(output_folder, 'filtered_output.tsv')
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")

    dados_processados = []
    header = ["ID protein", "Cálcio(Ca)", "Cobalto(Co)", "Cobre(Cu)", "Ferro(Fe)", "Pótassio(K)", "Magnésio(Mg)", "Manganês(Mn)", "Sódio(Na)", "Níquel(Ni)", "Zinco(Zn)"]

    with open(input_file, 'r') as f_in:
        # Lê o cabeçalho do arquivo e cria um DataFrame vazio
        df = pd.read_csv(f_in, sep='\t', nrows=0)
        df_result = pd.DataFrame(columns=header)  # DataFrame para os resultados
        
        # Reinicia a leitura do arquivo
        f_in.seek(0)
        
        # Processa cada linha do arquivo
        for linha in f_in:
            colunas = linha.strip().split('\t')
            id_protein = colunas[0]
            valores_processados = []
            for valor in colunas[1:]:
                try:
                    valor_float = float(valor)
                    if valor_float >= 0.8:
                        valores_processados.append('1')
                    else:
                        valores_processados.append('0')
                except ValueError:
                    print(f"Valor inválido encontrado: {valor}. Será tratado como '0'.")
                    valores_processados.append('0')
            nova_linha = [id_protein] + valores_processados
            dados_processados.append(nova_linha)

    # Cria o DataFrame com os dados processados
    df_result = pd.DataFrame(dados_processados, columns=header)
    
    # Salva o DataFrame em um novo arquivo .tsv
    arquivo_saida = os.path.join(output_folder, 'result_final.tsv')
    df_result.to_csv(arquivo_saida, sep='\t', index=False)
    # Removendo o arquivo input_metal.tsv
    # os.remove(input_file)
    return arquivo_saida




# input='/home/rey/Documentos/SALA DE TESTES/TESTE_reverse'
# tranformation_results(input)

