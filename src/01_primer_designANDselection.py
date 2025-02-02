from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import multiprocessing
import argparse
import subprocess
import sys
import os

file_name = '../output/sequencias.fasta'

num_sequences = 0

path = os.getcwd()

file_path = os.path.join(path, file_name)

################### PARALELIZM BLOCK

def dividir_arquivo(input_file, num_partes):
    """Divide o arquivo input_ePCR.tsv em várias partes"""
    with open(input_file, "r") as f:

        linhas = f.readlines()

    total_linhas = len(linhas)

    tamanho_parte = total_linhas // num_partes

    arquivos_divididos = []

    for i in range(num_partes):

        parte = linhas[i * tamanho_parte: (i + 1) * tamanho_parte] if i < num_partes - 1 else linhas[i * tamanho_parte:]

        parte_nome = f"../output/input_ePCR_part_{i}.tsv"

        with open(parte_nome, "w") as f_part:

            f_part.writelines(parte)

        arquivos_divididos.append(parte_nome)

    return arquivos_divididos

def executar_epcr(arquivo_parte, resultado_parte):

    """Executa o e-PCR em uma parte do arquivo"""
    comando_epcr = f"../bin/e-PCR {arquivo_parte} ../output/sequencias.fasta M=400 N=3 T=3 G=3 U=+ O={resultado_parte}"

    try:

        subprocess.run(comando_epcr, shell=True, check=True)

        print(f"Processo finalizado: {arquivo_parte}")

    except subprocess.CalledProcessError as e:

        print(f"Erro ao executar o e-PCR para {arquivo_parte}: {e}")

def rodar_epcr_em_paralelo(input_file, num_processos):

    """Divide o input e executa o e-PCR em paralelo"""
    if not os.path.exists(input_file):

        print(f"Erro: {input_file} não encontrado.")

        return

    # 1. Dividir arquivo
    arquivos_divididos = dividir_arquivo(input_file, num_processos)

    # 2. Criar lista de processos
    processos = []

    resultados_parciais = []

    for i, arquivo_parte in enumerate(arquivos_divididos):

        resultado_parte = f"../output/primerBlast_part_{i}.txt"

        resultados_parciais.append(resultado_parte)

        p = multiprocessing.Process(target=executar_epcr, args=(arquivo_parte, resultado_parte))

        processos.append(p)

    # 3. Iniciar processos
    for p in processos:

        p.start()

    # 4. Esperar processos terminarem
    for p in processos:

        p.join()

    # 5. Unir resultados em um único arquivo
    with open("../output/primerBlast.txt", "w") as resultado_final:

        for resultado_parte in resultados_parciais:

            if os.path.exists(resultado_parte):

                with open(resultado_parte, "r") as f:

                    resultado_final.writelines(f.readlines())

                os.remove(resultado_parte)  # Remover arquivos parciais para limpeza

    print("Execução paralela do e-PCR finalizada!")


####################

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

        seq_record = SeqRecord(Seq(str(record.seq).replace("-", "")), id=record.id, description="")
    
        # seq_record = SeqRecord(record.seq.replace("-", ""), id=record.id, description="")
        
        sequences.append(seq_record)

    # Salva as sequências individuais não alinhadas em formato FASTA
    SeqIO.write(sequences, "../output/sequencias.fasta", "fasta")


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
            f.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=../bin/primer3/src/primer3_config/\n") # caminho para o arquivo de parâmetros termodinâmicos, esse caminho precisa ser modificado para o local do computador do usuário
            f.write("=\n")
            
    # Caminho para o arquivo FASTA contendo as sequências
    fasta_file = "../output/sequencias.fasta"

    # Caminho para o arquivo de configuração "input.p3" que será gerado
    output_file = "../output/input.p3"

    num_sequences = 0

    # Loop para processar cada sequência no arquivo FASTA e gerar o arquivo de configuração
    with open(fasta_file) as f:
    
        for record in SeqIO.parse(f, "fasta"):
        
            sequence_id = record.id
            
            sequence_template = str(record.seq)
            
            generate_primer3_input(output_file, sequence_id, sequence_template)

    # Comando para executar o primer3_core com redirecionamento de entrada e saída
    command = "../bin/primer3/src/primer3_core < ../output/input.p3 > ../output/output.p3"

    try:
    
        # Executa o comando usando o subprocess
        subprocess.run(command, shell=True, check=True)
        
        print("primer3_core executado com sucesso.")
        
    except subprocess.CalledProcessError as e:
    
        print(f"Erro ao executar o primer3_core: {e}")

##############################################################################################################################
################################################ BEST PRIMER SELECTION #######################################################
##############################################################################################################################
if os.path.exists(file_path):

    print("##########################################################################")
    answer = input("WORNING! Se você responder 'sim' para a pergunta a seguir, o script considera que você já executou esse script antes e possui dois princiapis arquivos nesse diretório: 'primerBlast.txt' e 'sequencias.fasta'. Verifique se esses arquivos estão neste diretório.\n\nVocê já executou esse script anteriormente? Responta Yes ou No: ")
    
    if(answer.lower() == "yes"): 
    
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

            # print(fragment, dic_epcr)
            for i in fragment:
                # print(i)
                try:
                    i = i.strip()
                    i = i.split("\t")

                    #### PARTE DO CODIGO RESPONSAVEL POR FILTRAR PRIMERS COM BASE NO NUMERO DE GAR E MISMATCH. AQUI O DICIONARIO  dic_epcr{}
                    #### SERA CRIADO COM OS ALVOS DE CADA PRIMER E POSTERIORMENTE, MAIS A FRENTE DO CODIGO, O RECUEPRADO O NUMERO DE ALVOS
                    #### COM A FUNCAO len()
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
        
        # ITERAGE NO DICIONARIO dic_epcr{} PARA SBER O COMPRIMEITO DO VALOR DE CADA CHAVE. COM ISSO SERA OBTIDO O NUMERO DE ALVOS
        for k, v in dic_epcr.items():

            # print(k, v, len(v))
        
            dic_contagem[k] = len(v)
        # DICIONARIO QUERECEBE COMO VALOR, A CONTAGEM DE ALVOS PARA CADA PRIMER DE MODO CRESCENTE
        dic_contagem = dict(sorted(dic_contagem.items(), key=lambda item: item[1], reverse = True))
        
        primers = open("../output/output_p3.tsv", 'r').readlines()

        best_primers = []
        
        for indice, (k, v) in enumerate(dic_contagem.items()):
            
            if int(indice) > int(cutt_off):
            
                break

            best_primers.append((k, v))
            
            # print(k, v)
        #### AQUI ESTAMOS RELACIONANDO O NUMERO DE ALVOS DE CADA PRIMER, COM AS INFORMACOES DO ARQUIVO DE PRIEMRS DO PRIMERR 3
        #### A IDEIA E TER UMA VARIAVEL QUE GUARDE O NOME DO PRIMER E O NUMERO DE ALVOS E SABER TAMBEM OS PARAMETROS DE QUALIDADE
        #### OBTIDO NO PRIMER3 
        with open("../output/best_primers.tsv", 'w') as best_primer_out:

            for i in best_primers:

                for x in primers:

                    print(i, x)

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

    elif(answer.lower() == "no"):

        files = ["../output/input.p3", "../output/input_ePCR.tsv", "../output/output.p3", "../output/primerBlast.txt", "../output/sequencias.fasta"]

        for i in files:

            if os.path.exists(i):

                os.remove(i)

                print(f"Removido: {i}")

            else:

                print(f"Arquivo não encontrado: {i}")

        
        # subprocess.run("rm input.p3 input_ePCR.tsv output.p3 primerBlast.txt sequencias.fasta", shell=True)
        
        if __name__ == "__main__":
        
            input_file = "../output/output.p3"  # Coloque o nome do arquivo de entrada aqui
            
            output_file = "../output/output_p3.tsv"  # O resultado será salvo nesse arquivo
            
            output_file2 = "../output/input_ePCR.tsv"  # O resultado será salvo nesse arquivo no formato de entrada para a e-PCR
            
            run_primer3_epcr(input_file, output_file, output_file2)
            
            format_to_tabular(input_file, output_file, output_file2)

            num_processos = min(multiprocessing.cpu_count(), 4)  # Define o número de processos (4 ou número de CPUs)

            rodar_epcr_em_paralelo("../output/input_ePCR.tsv", num_processos)

            # comando_epcr = "bin/e-PCR output/input_ePCR.tsv output/sequencias.fasta M=400 N=3 T=3 G=3 U=+ O=output/primerBlast.txt"
            
            # try:
            
                # Executa o comando usando o subprocess
                
                # subprocess.run(comando_epcr, shell=True, check=True)
                
                # print("e-PCR executado com sucesso.")
                
            # except subprocess.CalledProcessError as e:
            
                # print(f"Erro ao executar o e-PCR: {e}")

else:

    if __name__ == "__main__":
    
        input_file = "../output/output.p3"  # Coloque o nome do arquivo de entrada aqui
        
        output_file = "../output/output_p3.tsv"  # O resultado será salvo nesse arquivo
        
        output_file2 = "../output/input_ePCR.tsv"  # O resultado será salvo nesse arquivo no formato de entrada para a e-PCR
        
        run_primer3_epcr(input_file, output_file, output_file2)
        
        format_to_tabular(input_file, output_file, output_file2)

        num_processos = min(multiprocessing.cpu_count(), 4)  # Define o número de processos (4 ou número de CPUs)

        rodar_epcr_em_paralelo("../output/input_ePCR.tsv", num_processos)
        
        # comando_epcr = "bin/e-PCR output/input_ePCR.tsv output/sequencias.fasta M=400 N=3 T=3 G=3 U=+ O=output/primerBlast.txt"
        
        # try:
        
            # Executa o comando usando o subprocess
            
            # subprocess.run(comando_epcr, shell=True, check=True)
            
            # print("e-PCR executado com sucesso.")
            
        
        # except subprocess.CalledProcessError as e:
        
            # print(f"Erro ao executar o e-PCR: {e}")




