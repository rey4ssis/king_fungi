import os
import docker
from config.celery_config import celery_app



@celery_app.task
def run_mymetal_container( deep_analiser_result, output):
    file, input = deep_analiser_result
    output_dir = f'{output}/Mymetal'
    os.makedirs(output_dir, exist_ok=True)
    client = docker.from_env()
    # print(file)
    # print(input)
    # print(file)
    container_params = {
        'image': 'mymetal:1.5',
        'volumes': {
            os.path.abspath(input): {'bind': f'/app/input', 'mode': 'rw'},
             os.path.abspath(output_dir): {'bind': f'/app/output/Mymetal', 'mode': 'rw'}
        },
        'command': [
            '--input-file',file],
        'remove': True,
    }
    try:
        print(f'[DEBUG] - Tentando executar container com parâmetros: {container_params}')
        client.containers.run(**container_params)
        # return True
    except Exception as e:
        print(f'[DEBUG] - Erro durante a execução do container: {str(e)}')
        return False, f'[MYMETAL STEP] - Unexpected error: {str(e)}'

# output = '/home/rey/Documentos/Git-tcc/kingfungi/output'
# input = "/home/rey/Documentos/Git-tcc/kingfungi/input"
# result = run_mymetal_container('input.fasta', input, output)
# print(result)


