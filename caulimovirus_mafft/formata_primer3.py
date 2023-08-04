from Bio import SeqIO

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
        f.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/emasilva/primer3/src/primer3_config/\n")
        f.write("=\n")

# Caminho para o arquivo FASTA contendo as sequências
fasta_file = "sequencias_cl1.fasta"

# Caminho para o arquivo de configuração "input.p3" que será gerado
output_file = "input_cl1.p3"

# Loop para processar cada sequência no arquivo FASTA e gerar o arquivo de configuração
with open(fasta_file) as f:
    for record in SeqIO.parse(f, "fasta"):
        sequence_id = record.id
        sequence_template = str(record.seq)
        generate_primer3_input(output_file, sequence_id, sequence_template)

