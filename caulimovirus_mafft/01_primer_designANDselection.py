from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import subprocess
import sys
import os

file_name = 'sequencias.fasta'

num_sequences = 0

path = os.getcwd()

file_path = os.path.join(path, file_name)

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
               
    f.close()
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
        
        print(output_file, output_file2)

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

def run_primer3_epcr(input_file, output_file, output_file2):
    arquivo_alinhamento = sys.argv[1]

    # Carrega o alinhamento múltiplo em formato FASTA
    alignment = AlignIO.read(arquivo_alinhamento, "fasta")

    # Cria uma lista de objetos SeqRecord para armazenar as sequências
    sequences = []
    
    for record in alignment:
    
        seq_record = SeqRecord(record.seq.replace("-", ""), id=record.id, description="")
        
        sequences.append(seq_record)

    # Salva as sequências individuais não alinhadas em formato FASTA
    SeqIO.write(sequences, "sequencias.fasta", "fasta")


    # Função para gerar o arquivo de configuração "input.p3" para cada sequência
    def generate_primer3_input(output_file, sequence_id, sequence_template):
        with open(output_file, "a") as f:
            f.write(f"SEQUENCE_ID={sequence_id}\n")
            f.write(f"SEQUENCE_TEMPLATE={sequence_template}\n")
            f.write("PRIMER_PRODUCT_SIZE_RANGE=500-600\n")
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
            f.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/emasilva/primer3/src/primer3_config/\n") # caminho para o arquivo de parâmetros termodinâmicos, esse caminho precisa ser modificado para o local do computador do usuário
            f.write("=\n")
            
    # Caminho para o arquivo FASTA contendo as sequências
    fasta_file = "sequencias.fasta"

    # Caminho para o arquivo de configuração "input.p3" que será gerado
    output_file = "input.p3"

    num_sequences = 0

    # Loop para processar cada sequência no arquivo FASTA e gerar o arquivo de configuração
    with open(fasta_file) as f:
    
        for record in SeqIO.parse(f, "fasta"):
        
            sequence_id = record.id
            
            sequence_template = str(record.seq)
            
            generate_primer3_input(output_file, sequence_id, sequence_template)

    # Comando para executar o primer3_core com redirecionamento de entrada e saída
    command = "~/primer3/src/primer3_core < input.p3 > output.p3"

    try:
    
        # Executa o comando usando o subprocess
        subprocess.run(command, shell=True, check=True)
        
        print("primer3_core executado com sucesso.")
        
    except subprocess.CalledProcessError as e:
    
        print(f"Erro ao executar o primer3_core: {e}")

if os.path.exists(file_path):

    print("##########################################################################")
    answer = input("WORNING! Se você responder 'sim' para a pergunta a seguir, o script considera que você já executou esse script antes e possui dois princiapis arquivos nesse diretório: 'primerBlast.txt' e 'sequencias.fasta'. Verifique se esses arquivos estão neste diretório.\n\nVocê já executou esse script anteriormente? Responta Yes ou No: ")
    
    if(answer == "Yes" or answer == "yes" or answer == "YES"): 
    
    # Loop para processar cada sequência no arquivo FASTA e gerar o arquivo de configuração
    
        with open(file_path) as f:
        
            for record in SeqIO.parse(f, "fasta"):
            
                sequence_id = record.id
                
                sequence_template = str(record.seq)
                
                num_sequences += 1
                
        file_epcr = open(sys.argv[2], 'r').readlines()

        print("##########################################################################")
        cutt_off = input(f"Quantos primers você deseja ter de retorno? Digite um valor numerico inteiro: ")
        print("##########################################################################")        
        

        def conta_primer(fragment, dic_epcr, cutt_off):
            for i in fragment:
                try:
                    i = i.strip()
                    i = i.split("\t")
                    if int(i[6]) <= 2 and int(i[7]) <= 2:
                        dic_epcr.setdefault(i[1], []).append(i[0])
                except Exception as e:
                    print("Ocorreu um erro:", str(e))

        dic_epcr = {}

        fragment_size = len(file_epcr) // 3

        fragment_1 = file_epcr[:fragment_size]
        
        conta_primer(fragment_1, dic_epcr, cutt_off)
        
        del fragment_1

        file_epcr = file_epcr[fragment_size:]

        fragment_size = len(file_epcr) // 2

        fragment_2 = file_epcr[:fragment_size]
        
        conta_primer(fragment_2, dic_epcr, cutt_off)

        del fragment_2

        fragment_3 = file_epcr[fragment_size:]

        conta_primer(fragment_3, dic_epcr, cutt_off)

        del file_epcr

        del fragment_3

        dic_contagem = {}

        for k, v in dic_epcr.items():
        
            dic_contagem[k] = len(v)

        dic_contagem = dict(sorted(dic_contagem.items(), key=lambda item: item[1], reverse = True))
        
        primers = open("output_p3.tsv", 'r').readlines()

        best_primers = []
        
        for indice, (k, v) in enumerate(dic_contagem.items()):
            
            if int(indice) > int(cutt_off):
            
                break

            best_primers.append((k, v))
            
#            print(k, v)
        
        with open("best_primers.tsv", 'w') as best_primer_out:
            for i in best_primers:

                for x in primers:
                
                    linha_primers = []

                    x = x.strip()
                
                    line = x.split("\t")
                    
                    linha_primers.append(i[0])
                    
                    linha_primers.append(i[1])

                    if i[0] == line[0]:
                    
                        for k in range(len(line)):
                            
                            linha_primers.append(line[k])

                        best_primer_out.write("\t".join(map(str, linha_primers))+"\n")
                
                        linha_primers = []

    elif(answer == "no" or answer == "No" or answer == "NO"):
    
        subprocess.run("rm input.p3 input_ePCR.tsv output.p3 primerBlast.txt sequencias.fasta", shell=True)
        
        if __name__ == "__main__":
        
            input_file = "output.p3"  # Coloque o nome do arquivo de entrada aqui
            
            output_file = "output_p3.tsv"  # O resultado será salvo nesse arquivo
            
            output_file2 = "input_ePCR.tsv"  # O resultado será salvo nesse arquivo no formato de entrada para a e-PCR
            
            run_primer3_epcr(input_file, output_file, output_file2)
            
            format_to_tabular(input_file, output_file, output_file2)
            
            comando_epcr = "e-PCR input_ePCR.tsv sequencias.fasta M=400 N=3 T=3 G=3 U=+ O=primerBlast.txt"
            
            try:
            
                # Executa o comando usando o subprocess
                
                subprocess.run(comando_epcr, shell=True, check=True)
                
                print("e-PCR executado com sucesso.")
                
            except subprocess.CalledProcessError as e:
            
                print(f"Erro ao executar o e-PCR: {e}")

else:

    if __name__ == "__main__":
    
        input_file = "output.p3"  # Coloque o nome do arquivo de entrada aqui
        
        output_file = "output_p3.tsv"  # O resultado será salvo nesse arquivo
        
        output_file2 = "input_ePCR.tsv"  # O resultado será salvo nesse arquivo no formato de entrada para a e-PCR
        
        run_primer3_epcr(input_file, output_file, output_file2)
        
        format_to_tabular(input_file, output_file, output_file2)
        
        comando_epcr = "e-PCR input_ePCR.tsv sequencias.fasta M=400 N=3 T=3 G=3 U=+ O=primerBlast.txt"
        
        try:
        
            # Executa o comando usando o subprocess
            
            subprocess.run(comando_epcr, shell=True, check=True)
            
            print("e-PCR executado com sucesso.")
            
        except subprocess.CalledProcessError as e:
        
            print(f"Erro ao executar o e-PCR: {e}")
