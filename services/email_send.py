from flask import Flask, request, render_template
import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.mime.base import MIMEBase
from email import encoders
from dotenv import load_dotenv
from config.celery_config import celery_app
from jinja2 import Template
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import zipfile

# Configurações de e-mail
EMAIL_ADDRESS = os.getenv('EMAIL_ADDRESS')
EMAIL_PASSWORD = os.getenv('EMAIL_PASSWORD')



# def count_ones(column):
#     return column.eq(1).sum()

# def create_graphs(df, file_path):
#     counts = df.iloc[:, 1:].apply(count_ones)
#     total_rows = len(df)
#     paths = []

#     folder_path = os.path.dirname(file_path)  # Obtenha o caminho da pasta do arquivo
#     for col, count in counts.items():
#         plt.figure(figsize=(8, 6))
#         plt.bar(['Quantidade de proteínas', f'Total de {col}'], [total_rows, count])
#         plt.text(0, total_rows, total_rows, ha='center', va='bottom')
#         plt.text(1, count, f'{count} ({count/total_rows*100:.2f}%)', ha='center', va='bottom')
#         plt.title(f'Relação entre o total de linhas e o número de 1s em {col}')
#         plt.ylabel('Quantidade')
#         graph_file = os.path.join(folder_path, f'{col}_graph.png')  # Caminho do arquivo de gráfico
#         plt.savefig(graph_file)
#         paths.append(graph_file)
#         plt.close()

#     zip_file_path = os.path.join(folder_path, 'graphs.zip')
#     with zipfile.ZipFile(zip_file_path, 'w') as zipf:
#         for path in paths:
#             zipf.write(path, os.path.basename(path))

#     return zip_file_path, paths

def plot_graph_all(file):
    # Extrair o diretório do arquivo de entrada
    output_dir = os.path.dirname(file)
    
    df = pd.read_csv(file, sep='\t', skiprows=3, header=0)
    num_prot= len(df)
    df = df.drop(columns=['ID protein'])
    df = df.apply(pd.to_numeric, errors='coerce')
    sums = df.sum()
    
    # Criar cores para todas as barras
    colors = ['#06C' if i in [0, 1, 2, 3, 4, 5, 7, 9] else '#8BC1F7' for i in range(len(sums))]
    
    plt.figure(figsize=(10, 6))
    
    # Ajusta o espaçamento para acomodar a legenda
    plt.subplots_adjust(right=0.75)
    
    bars = plt.bar(sums.index, sums.values, color=colors)
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval, round(yval, 2), va='bottom')
    
    plt.title('Quantidade de proteínas transmenbranas que se ligam a metais')
    plt.xlabel(f'Número de proteínas transmenbranas encontradas: {num_prot}')
    plt.ylabel('Quantidade de proteínas')
    plt.xticks(rotation=45)
    
    # Define a legenda para as cores azuis
    legendas_azuis = [plt.Rectangle((0,0),1,1, color='#06C')]
    leg_labels_azuis = ['Micronutrientes']
    
    # Define a legenda para as cores vermelhas
    legendas_vermelhas = [plt.Rectangle((0,0),1,1, color='#8BC1F7')]
    leg_labels_vermelhas = ['Elementos essenciais']
    
    # Define o tamanho da legenda e sua posição
    plt.legend(legendas_azuis + legendas_vermelhas, leg_labels_azuis + leg_labels_vermelhas, loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
    
    plt.tight_layout()
    plt_file = os.path.join(output_dir, 'grafico.png')  # Salvar na mesma pasta do arquivo de entrada
    plt.savefig(plt_file)
    plt.close()
    return plt_file



@celery_app.task
def send_email(file_path, nome, email):
    nome= nome
    email= email
    msg = MIMEMultipart()
    msg['From'] = EMAIL_ADDRESS
    msg['To'] = email
    msg['Subject'] = 'Resultados da Análise genomica - KingFungi(G2BC)'
   

    html_template = '''
        <!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <title>E-mail G2BC</title>
        </head>
        <body>

            <div style="margin-left: 15%; margin-right: 15%; text-align: center; margin-top: 1%;">
                <p style="font-size: 20px; font-family: Verdana;"><b><i>INFORMATIVO DA ANÁLISE</i></b></p>
            </div>
            <div style="background-color: #4682B4; text-align: center; margin-left: 5%; margin-right: 5%; margin-top: 5px; padding: 10px;">
                <p style="font-size: 18px; color: white; font-family: Verdana;"><b>Olá, Sr. {{ nome }}</b></p>
                <p style="font-size: 14px; color: white; font-family: Verdana;">Estamos entrando em contato para avisar que os resultados da sua análise para identificar proteínas transmembrana que têm capacidade de se ligar a íons metálicos foram concluídos. 
                <br>
                O anexo com os dados está incluído no corpo deste e-mail.</p>
                <p style="font-size: 14px; color: white; font-family: Verdana;">
                Muito obrigado por utilizar o nosso software.
                </p>

                <p style="text-align: center; font-family: Verdana; font-size: 14px; color: #ffffff">
                    G2BC - Grupo de Pesquisa em Bioinformática e Biologia Computacional
                    <br>
                    Este <b>e-mail foi enviado automaticamente</b> pelo sistema. <b>Não é necessário respondê-lo</b>. Em caso de dúvidas, responda este e-mail ou entre em contato pelo telefone.
                    <br>
                    E-mail: <b>kingfungi.uneb@gmail.com</b>
                    <br>
                    Telefone: <b>(71) 3117-2274</b>
                </p>
                
            </div>
        </body>
        </html>
        '''
    html_rendered = Template(html_template).render(nome=nome)
    msg.attach(MIMEText(html_rendered, 'html'))


    # Criar gráficos e zipar
    # zip_file_path, paths = create_graphs(pd.read_csv(file_path, sep='\t', skiprows=3, header=0), file_path)
    file_graph = plot_graph_all(file_path)

    # Adicionar anexos ao e-mail
    # anexos = [zip_file_path, file_path, file_graph]
    anexos = [file_path, file_graph]
    for anexo in anexos:
        part = MIMEBase('application', 'octet-stream')
        with open(anexo, 'rb') as arquivo:
            part.set_payload(arquivo.read())
        encoders.encode_base64(part)
        nome_arquivo = os.path.basename(anexo)
        part.add_header('Content-Disposition', f"attachment; filename= {nome_arquivo}")
        msg.attach(part)

    # Conexão com o servidor SMTP

    with smtplib.SMTP('smtp.gmail.com', 587) as server:
        server.starttls()
        server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
        server.sendmail(EMAIL_ADDRESS, email, msg.as_string())
    os.remove(file_graph)

# # Exemplo de uso:
# arquivo_csv = "ids.tsv"
# nome = "Dr. Fulano"
# email = "reyassis7@gmail.com"
# send_email(arquivo_csv, nome, email)

