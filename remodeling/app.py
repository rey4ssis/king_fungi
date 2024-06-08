import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.mime.base import MIMEBase
from email import encoders
from flask import Flask, request, send_file, render_template, jsonify
from services.email_send import send_email
from Bio import SeqIO
import traceback
import json
from config.celery_config import celery_app
from celery import chain
from helper import *
from services.analiser_results_deep import extract_tm_seqs_and_filter_tsv
from services.run_container_mymetal import run_mymetal_container
from services.genbank_service import download_fasta_sequence_by_id
from services.genbank_create_id import generate_fungi_by_name
from services.analiser_results_mymetal import tranformation_results
from services.create_results import create_result
from services.deeptmhmm import processDeep
from services.pipeline import run_mymetal_filter_task
from config.run_commands import *
import time  

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        if 'file_sequence' not in request.files:
            return jsonify({'success': False, 'message': 'Nenhum arquivo enviado!'})

        nome = request.form.get("name_sequence")
        email = request.form.get("email_sequence")
        file = request.files['file_sequence']
        is_sequence_ref = request.form.get('ref_busca')
        
        # Medição do tempo de execução
        start_time = time.time()
        
        # Validações e processamento
        try:
            # Validates user email
            validate_email(email)

            # Create user directories
            input_dir, output_dir = create_user_dirs(email)
            path_out = email

            # Handle sequence reference
            file_down = []
            if verification_input(is_sequence_ref) == 0:
                if is_sequence_ref:
                    file_down, _ = download_fasta_sequence_by_id(is_sequence_ref, path_out)
            elif verification_input(is_sequence_ref) == 1:
                rt_generate, status_code = generate_fungi_by_name(is_sequence_ref)
                genbank_acession = rt_generate.get('genbank_assembly_accession')
                refseq_acession = rt_generate.get('refseq_assembly_accession')
                if refseq_acession == ' ':
                    file_down, _ = download_fasta_sequence_by_id(genbank_acession, path_out)
                else:
                    file_down, _ = download_fasta_sequence_by_id(refseq_acession, path_out)

            # Validate and save the file
            if file and file.filename.endswith(('.fasta', '.fna', 'faa')):
                file_path = os.path.join(input_dir, file.filename)
                file.save(file_path)
                file_down = file.filename

            # Run the analysis tasks
            chain(
                run_mymetal_filter_task.s(file_down, input_dir, output_dir),
                processDeep.s(input_dir, f'{output_dir}/Deep'),
                extract_tm_seqs_and_filter_tsv.si(output_dir, output_dir),
                tranformation_results.si(output_dir),
                create_result.s(file_down),
                send_email.s(nome, email)
            )()

            # Cleanup
            delete_user_dirs(email)

            # Measure execution time
            end_time = time.time()
            total_time = end_time - start_time
            print("Tempo total de execução:", total_time)

            return jsonify({'success': True, 'message': 'Arquivo enviado com sucesso e removido da pasta uploads!'})
        except Exception as e:
            print(traceback.format_exc())
            return jsonify({'success': False, 'message': str(e)})

    return render_template('index.html')

if __name__ == "__main__":
    app.run(debug=True)
