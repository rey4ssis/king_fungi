import os
from celery import Celery
import biolib
from mymetal import mbp, iof
import shutil
from celery import Celery, chain
from dotenv import load_dotenv


load_dotenv(override=True)

celery_app = Celery(
    'celery_config',
    broker=os.getenv('BROKER_URL'),
    backend=os.getenv('CELERY_RESULT_BACKEND'),
    include=[
        'services.analiser_results_deep',
        'services.analiser_results_mymetal',
        'services.genbank_service',
        'services.email_send',
        'services.create_results',
        'services.deeptmhmm',
        'services.pipeline'
    ],
)
celery_app.conf.update(
    CELERYD_FORCE_EXECV=True
)
