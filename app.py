from flask import Flask, request, send_file
from tasks import process_metal_task, process_tmhmm_task

app = Flask(__name__)

@app.route('/processMetal', methods=['POST'])
def processMetal():
    if 'file' not in request.files:
        return "Nenhum arquivo enviado", 400

    file = request.files['file']
    if file.filename.endswith('.fasta'):
        file.save('input_metal.fasta')
        process_metal_task.delay('input_metal.fasta')
        return "Processamento iniciado em segundo plano", 202

    return "Arquivo inválido", 400

@app.route('/processTMHMM', methods=['POST'])
def processTMHMM():
    if 'file' not in request.files:
        return "Nenhum arquivo enviado", 400

    file = request.files['file']
    if file.filename.endswith('.fasta'):
        file.save('input.fasta')
        process_tmhmm_task.delay('input.fasta')
        return "Processamento iniciado em segundo plano", 202

    return "Arquivo inválido", 400

if __name__ == "__main__":
    app.run(debug=True)
