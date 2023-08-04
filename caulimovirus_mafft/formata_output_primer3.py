def format_to_tabular(input_file, output_file):
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
    with open(output_file, 'w') as out_f:
        out_f.write(header)

        # Escrever os dados formatados no arquivo de saída
        for k, v in data_formated.items():
            line = [k]
            line.extend(v)
            out_f.write('\t'.join(line) + '\n')


if __name__ == "__main__":
    input_file = "output_cl1.p3"  # Coloque o nome do arquivo de entrada aqui
    output_file = "output_output_cl1.p3"  # O resultado será salvo nesse arquivo

    format_to_tabular(input_file, output_file)
