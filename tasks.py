from celery import Celery
import biolib
from mymetal import mbp, iof
import shutil

celery = Celery(__name__, broker='amqp://guest@localhost//')  # Configure conforme necessário

@celery.task
def process_metal_task(file_name):
    mbp_preb = mbp.predict(file_name)
    iof.save_out_csv(mbp_preb, 'input_metal.csv')

@celery.task
def process_tmhmm_task(file_name):
    deeptmhmm = biolib.load('DTU/DeepTMHMM')
    deeptmhmm_job = deeptmhmm.cli(args=f'--fasta {file_name}')
    deeptmhmm_job.save_files('output')
    shutil.make_archive('results', 'zip', '.', 'output')
