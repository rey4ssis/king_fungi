import subprocess
import psutil
from time import sleep
from config.celery_config import celery_app

def check_rabbitmq_running():
    try:
        # Verifica se o container RabbitMQ está em execução
        subprocess.run(["docker", "inspect", "rabbitmq_container"], check=True, stdout=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False

def run_rabbitmq_celery():
    if not check_rabbitmq_running():
        print("Iniciando RabbitMQ...")
        subprocess.run(["docker", "run", "-d", "-p", "5672:5672", "--name", "rabbitmq_container", "rabbitmq"])
        sleep(10)
    
    print("Iniciando Celery Worker e Celery Flower...")
    subprocess.Popen(["celery", "-A", "config.celery_config", "worker", "--loglevel=info"])
    subprocess.Popen(["celery", "-A", "config.celery_config", "flower", "--loglevel=info"])

def stop_services():
    print("Parando Celery Worker e Celery Flower...")
    for proc in psutil.process_iter():
        if "celery" in proc.name():
            proc.terminate()
            
    print("Parando e removendo o container RabbitMQ...")
    subprocess.run(["docker", "stop", "rabbitmq_container"])
    subprocess.run(["docker", "rm", "rabbitmq_container"])

if __name__ == '__main__':
    run_rabbitmq_celery()
    # Aqui você pode adicionar uma lógica para executar algumas tarefas ou testes
    sleep(60)  # Por exemplo, aguarde 1 minuto antes de parar os serviços
    stop_services()
