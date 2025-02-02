from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import subprocess

# Criando um parser de argumentos
parser = argparse.ArgumentParser(description='Processa um alinhamento múltiplo e gera sequências individuais não alinhadas em formato FASTA.')
parser.add_argument('arquivo_alinhamento', help='O caminho para o arquivo de alinhamento no formato FASTA.')

# Obtendo o arquivo de alinhamento a partir dos argumentos
args = parser.parse_args()
arquivo_alinhamento = args.arquivo_alinhamento

# Carrega o alinhamento múltiplo em formato FASTA
alignment = AlignIO.read(arquivo_alinhamento, "fasta")

# Cria uma lista de objetos SeqRecord para armazenar as sequências
sequences = []
for record in alignment:
    seq_record = SeqRecord(Seq(str(record.seq).replace("-", "")), id=record.id, description="")
    # seq_record = SeqRecord(record.seq.replace("-", ""), id=record.id, description="")
    sequences.append(seq_record)

# Salva as sequências individuais não alinhadas em formato FASTA
SeqIO.write(sequences, "sequencias.fasta", "fasta")


# Função para gerar o arquivo de configuração "input.p3" para cada sequência
def generate_primer3_input(output_file, sequence_id, sequence_template):
    with open(output_file, "a") as f:
        f.write(f"SEQUENCE_ID={sequence_id}\n")
        f.write(f"SEQUENCE_TEMPLATE={sequence_template}\n")
        f.write("PRIMER_PRODUCT_SIZE_RANGE=300-400\n")
        f.write("PRIMER_NUM_RETURN=100\n")
        f.write("PRIMER_OPT_SIZE=20\n")         # Tamanho ótimo do primer
        f.write("PRIMER_MIN_SIZE=18\n")         # Tamanho mínimo do primer
        f.write("PRIMER_MAX_SIZE=25\n")         # Tamanho máximo do primer
        f.write("PRIMER_OPT_TM=60.0\n")         # Temperatura ótima de anelamento
        f.write("PRIMER_MIN_TM=57.0\n")         # Temperatura mínima de anelamento
        f.write("PRIMER_MAX_HAIRPIN_TH=0.0\n")
        f.write("PRIMER_RIGHT_HAIRPIN_TH=0.0\n")
        f.write("PRIMER_LEFT_HAIRPIN_TH=0.0\n")
        f.write("PRIMER_INTERNAL_MAX_HAIRPIN_TH=0.0\n")
        f.write("PRIMER_MAX_TM=63.0\n")         # Temperatura máxima de anelamento
        f.write("PRIMER_OPT_GC_PERCENT=50.0\n") # Conteúdo ótimo de GC (%)
        f.write("PRIMER_MIN_GC=40.0\n")         # Conteúdo mínimo de GC (%)
        f.write("PRIMER_MAX_GC=60.0\n")         # Conteúdo máximo de GC (%)
        f.write("PRIMER_MAX_SELF_ANY_TH=1.0\n")# Máxima energia de emparelhamento do primer
        f.write("PRIMER_MAX_SELF_END_TH=1.0\n")# Máxima energia de emparelhamento no extremo 3'
        f.write("PRIMER_PAIR_MAX_COMPL_ANY_TH=1.0\n") # Máxima energia de emparelhamento do par de primers
        f.write("PRIMER_PAIR_MAX_COMPL_END_TH=1.0\n") # Máxima energia de emparelhamento do par de primers no extremo 3'
        f.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/blue/munoz/deandradesilvae/colaborations/breno/transposable_elements/primer3/src/primer3_config/\n")
        f.write("=\n")

# Caminho para o arquivo FASTA contendo as sequências
fasta_file = "sequencias.fasta"

# Caminho para o arquivo de configuração "input.p3" que será gerado
output_file = "input.p3"

# Loop para processar cada sequência no arquivo FASTA e gerar o arquivo de configuração
with open(fasta_file) as f:
    for record in SeqIO.parse(f, "fasta"):
        sequence_id = record.id
        sequence_template = str(record.seq)
        generate_primer3_input(output_file, sequence_id, sequence_template)

# Comando para executar o primer3_core com redirecionamento de entrada e saída
# command = "primer3_core < input.p3 > output.p3"

# try:
    # Executa o comando usando o subprocess
    # subprocess.run(command, shell=True, check=True)
    # print("primer3_core executado com sucesso.")
# except subprocess.CalledProcessError as e:
    # print(f"Erro ao executar o primer3_core: {e}")


def format_to_tabular(input_file, output_file, output_file2):
    # Abrir o arquivo de entrada e ler todas as linhas
    with open(input_file, 'r') as f:
        lines = f.readlines()

    data = {}  # Dicionário para armazenar os dados formatados
    data_formated = {}  # Dicionário para armazenar os dados formatados por chave
    current_key = ""  # Variável para armazenar a chave atual sendo processada
    for line in lines:
        line = line.strip()  # Remover espaços em branco do início e do fim da linha
        if not line.startswith("="):
            if line.startswith("SEQUENCE_ID"):
                current_key = line.split("=")[1]  # Obter a chave atual a partir da linha
                data[current_key] = {"SEQUENCE_ID": current_key}  # Criar uma entrada para a chave no dicionário 'data'
            else:
                key, value = line.split("=")  # Obter a chave e o valor da linha
                data[current_key][key] = value  # Adicionar a chave e o valor ao dicionário da chave atual

    # 'data' é o dicionário que contém os dados fornecidos

    # Loop para criar um novo dicionário (data_formated) com as informaç~eos de interesse
    for key in data:
        for k, v in data[key].items():
            splited = k.split("_")  # Dividir a chave em partes usando o caractere '_'
            if(len(splited) > 2):  # Verificar se a chave tem mais de duas partes
                if(splited[2].isnumeric() == True ):  # Verificar se a terceira parte é numérica
                    formated_key = key + "_" + str(splited[2])  # Criar uma chave formatada
                    data_formated.setdefault(formated_key, []).append(v)  # Adicionar valor à chave formatada no dicionário 'data_formated'

    # Cabeçalho das colunas no arquivo de saída
    header = "SEQUENCE_ID\tPRIMER_PAIR_PENALTY\tPRIMER_LEFT_PENALTY\tPRIMER_RIGHT_PENALTY\tPRIMER_LEFT_SEQUENCE\tPRIMER_RIGHT_SEQUENCE\tPRIMER_LEFT\tPRIMER_RIGHT\tPRIMER_LEFT_TM\tPRIMER_RIGHT_TM\tPRIMER_LEFT_GC_PERCENT\tPRIMER_RIGHT_GC_PERCENT\tPRIMER_LEFT_SELF_ANY_TH\tPRIMER_RIGHT_SELF_ANY_TH\tPRIMER_LEFT_SELF_END_TH\tPRIMER_RIGHT_SELF_END_TH\tPRIMER_LEFT_HAIRPIN_TH\tPRIMER_RIGHT_HAIRPIN_TH\tPRIMER_LEFT_END_STABILITY\tPRIMER_RIGHT_END_STABILITY\tPRIMER_PAIR_COMPL_ANY_TH\tPRIMER_PAIR_COMPL_END_TH\tPRIMER_PAIR_PRODUCT_SIZE\n"

    # Escrever o cabeçalho no arquivo de saída
    with open(output_file, 'w') as out_f, open(output_file2, 'w') as out_f2:
        out_f.write(header)
        
        # Escrever os dados formatados no arquivo de saída
        for k, v in data_formated.items():
            line = [k]
            line2 = []
            line.extend(v)
            line2.append(line[0])
            line2.append(line[4])
            line2.append(line[5])
            line2.append(line[22])
            out_f.write('\t'.join(line) + '\n')
            out_f2.write('\t'.join(line2) + '\n')

if __name__ == "__main__":
    input_file = "output.p3"  # Coloque o nome do arquivo de entrada aqui
    output_file = "output_p3.tsv"  # O resultado será salvo nesse arquivo
    output_file2 = "input_ePCR.tsv"  # O resultado será salvo nesse arquivo no formato de entrada para a e-PCR
    format_to_tabular(input_file, output_file, output_file2)
    comando_epcr = "./e-PCR input_ePCR.tsv sequencias.fasta M=400 N=3 T=3 G=3 U=+ O=primerBlast.txt"
    try:
        # Executa o comando usando o subprocess
        subprocess.run(comando_epcr, shell=True, check=True)
        print("e-PCR executado com sucesso.")
    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar o e-PCR: {e}")
