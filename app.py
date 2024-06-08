import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.mime.base import MIMEBase
from email import encoders
from flask import Flask, request, send_file, render_template
from services.email_send import send_email
from Bio import SeqIO
import traceback
import json
from config.celery_config import celery_app
from celery import chain
from helper import *
from services.analiser_results_deep import extract_tm_seqs
from services.deeptmhmm import run_deep_container
from services.run_container_mymetal import run_mymetal_container
from services.genbank_service import download_fasta_sequence_by_id
from services.genbank_create_id import generate_fungi_by_name
from services.analiser_results_mymetal import tranformation_results
from services.create_results import create_result
from kingfungi.services.deeptmhmm import processDeep
from config.run_commands import *
import time  


app = Flask(__name__)

#file_t = "/home/rey/Documentos/Git-tcc/input.fasta"
@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        if 'file_sequence' not in request.files:
            return 'Nenhum arquivo enviado!'

        # fast_sequence = request.form.get("ref_busca")        
        nome = request.form.get("name_sequence")
        email = request.form.get("email_sequence")
        file = request.files['file_sequence']
        is_sequence_ref = request.form.get('ref_busca')
        
        # Medição do tempo de execução
        start_time = time.time()
        #run_rabbitmq_celery()
        #validates user email
        validate_email(email)
        #check if user sent a sequence file or sequence ref
        #   if is a sequence ref downlaod sequence file from genbank
        input_dir, output_dir = create_user_dirs(email)
        path_out = email
        
        file_down= []
        #Verificar se o input é o nome ou acession number
        if verification_input(is_sequence_ref)== 0: #caso for acession number
            if is_sequence_ref:
                (
                    file_down,_
                ) = download_fasta_sequence_by_id(is_sequence_ref, path_out)
        if verification_input(is_sequence_ref)== 1: #caso for nome ---- PROBLEMA , Não baixar a porra do genoma certo
            rt_generate, status_code = generate_fungi_by_name(is_sequence_ref)
            genbank_acession= rt_generate.get('genbank_assembly_accession')
            refseq_acession = rt_generate.get('refseq_assembly_accession')
            
            if refseq_acession == ' ':
                file_down , _ = download_fasta_sequence_by_id(genbank_acession, path_out)
            
            else:
                file_down , _ = download_fasta_sequence_by_id(refseq_acession, path_out)
                
        #validates file extension
            verification_fasta(file)
                
        #validates file extension
            verification_fasta(file)
        
        #create dirs to stores user files
        # input_dir, output_dir = create_user_dirs(email)
       
        #save the user request fasta file
        #fasta_filepath = generate_unique_filepath(input_dir, 'fasta')
        file_path = os.path.join(input_dir, file.filename)
        if file and file.filename.endswith('.fasta' or '.fna' or 'faa'):
            
            # print(file_path)
            file.save(file_path)
            file_down = file.filename
                        
        #   run_deep_container.si(file_down, input_dir, f'{output_dir}/Deep'),
        chain(
            processDeep.s(file_down, input_dir,f'{output_dir}/Deep'),
            extract_tm_seqs.s(f'{output_dir}/Deep/', output_dir),
            run_mymetal_container.s(output=output_dir),
            tranformation_results.si(output_dir),
            # extract_headers.s(input_dir,file_down),
            # update_ids.s(input_dir,file_down),
            create_result.s(file_down),
            send_email.s(nome, email)
        )()
        
        # Envia o arquivo por e-mail
        delete_user_dirs(email)
        # Remove o arquivo da pasta "uploads"
        #os.remove(file_path)
        # Medição do tempo de execução
        end_time = time.time()
        total_time = end_time - start_time
        print("Tempo total de execução:", total_time)

        return 'Arquivo enviado com sucesso e removido da pasta uploads!'
    return render_template('index.html')


if __name__ == "__main__":
    app.run(debug=True)