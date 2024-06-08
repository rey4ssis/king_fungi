from Bio import SeqIO, Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
import config.celery_config as celery_config 
import os
import time
from config.celery_config import celery_app


@celery_app.task
def extract_tm_seqs(tmp, input_dir, output_dir):
  """
  Função para extrair sequências com TM do arquivo .3line e gerar novo arquivo .fasta.

  Args:
    input_file: Caminho para o arquivo .3line.
    output_file: Caminho para o novo arquivo .fasta.

  Returns:
    None.
  """

  input_file = os.path.join(input_dir,'predicted_topologies.3line')
  output_file = ('deep_out.fasta')
  output_file_path = os.path.join(f'{output_dir}/Deep',output_file)
  with open(input_file, 'r') as f_input, open(output_file_path, 'w') as f_output:
        seq_id = None
        sequence = None
        for line in f_input:
            line = line.strip()
            if line.startswith('>'):
                line=line[1:]
                terms = line.split('|')
                if terms[-1].strip() == 'TM' or terms[-1].strip() == 'SP+TM':
                    seq_id = terms[-2].strip()
            else:
                if seq_id is not None:
                    sequence = line.strip()
                    f_output.write(f'>{seq_id}\n{sequence}\n')
                    seq_id = None
                    sequence = None
  return output_file, f'{output_dir}/Deep'

# if __name__ == "__main__":
#   # Caminho para o arquivo .3line
#   input_file = "output/Deep/predicted_topologies.3line"
#   # Caminho para o novo arquivo .fasta
#   output_file = "output.fasta"

#   extract_tm_seqs(input_file, output_file)







