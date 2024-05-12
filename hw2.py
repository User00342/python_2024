import pandas as pd
import matplotlib.pyplot as plt
import logging


# подготовка материала

gene_map = pd.read_csv('ЯКарта.csv')

# преобразование файла, нужно будет сохранить и использовать один и тот же, названия столбцов поменять
gene_map.rename(columns = {'Unnamed: 2':'AK'}, inplace = True )
for i in range(len(gene_map.Кодон)):
    for j in range(len(gene_map.Кодон)):
        if gene_map.Аминокислота[i] == gene_map['Аминокислота.1'][j]:
            gene_map.AK[i] = gene_map.Сокращение[j]


def is_ak(seqs):
# Проверка являются ли последовательности аминокислотами
    alphabet_AK = set('ARNDCQEGHILKMFPSTWYV')
    answer = []
    for seq in seqs:
        seq = set(seq.upper())
        if seq <= alphabet_AK:
            answer.append('jup')
        else:
            answer.append('no, don`t do that')
    return answer



# вспомогательные функции

def DNA_karta(seq): 
# функция восстановления последовательности, подбор триплетов с помощью генетической карты (gene_map)
    all_variants = 1
    DNA_seq = ''
    for i in seq:
        variants = list(gene_map.Кодон[gene_map.AK == i])
        all_variants *= len(variants)
        DNA_seq = DNA_seq + variants[0]
    print('Данная аминокислота могла получиться из', all_variants, "различных последовательностей")
    print('Получившаяся последовательность:', DNA_seq)
    return DNA_seq


def letter_pie(SEQ):
    # Подсчет количества аминокислот в заданной последователсти seq и визуализация результата
    # составление словарика(AK_counts_dict) с количеством каждой аминокислоты в последовательности
    AK_counts_dict = {}
    for letter in SEQ:
        if letter in AK_counts_dict:
            AK_counts_dict[letter] += 1
        else:
            AK_counts_dict[letter] = 1

    message = ('В последовательности ' + SEQ + ' - ' + str(len(AK_counts_dict)) + ' видов аминокислот.')

    # мини-визуализация

    plt.show(plt.pie(list(AK_counts_dict.values()), labels=list(AK_counts_dict.keys()), autopct='%1.1f%%', shadow=True,
            wedgeprops={'width': 0.2, 'lw': 1, 'ls': '--', 'edgecolor': "k"}))

    return message

def palindrom(seq):
    # функция проверяет, является ли последовательность палиндромом
    ka = seq[:(len(seq)-1)//2:-1]
    ak = seq[:len(seq)//2]
    if ak == ka:
        result = 'palindrom'
    else:
        result = "isn't palindrom"
    return result

def Needlman_Wunsch(ak1, ak2):
    #     Воссоздание знаменитого алгоритма Нидлмана-Вунша, его нечитабельная версия.
    #     al_frame - матрица схожести; path - словарь пути алгоритма

    #     Функция словно не работала. Затем при проверке работала. А потом опять не работала. Я не понимаю, работает она или нет.
    #     Код слишком сумбурен, его следует разделить на несколько функций и оптимизировать. Думала разукрасить ячейки пути, но не получилось :(

    ak1 = ak1.lower()
    ak2 = ak2.lower()
    word1 = ['']
    word2 = ['']
    path = {'00': '00'}
    result_al1 = ''
    result_al2 = ''
    
    for letter in ak1:
        word1.append(letter)
    for letter in ak2:
        word2.append(letter)

    al_frame = pd.DataFrame(columns=word1, index=word2)
    al_frame.iloc[0][0] = 0
    for i in range(1, len(word1)):
        al_frame.iloc[0, i] = al_frame.iloc[0, i - 1] - 1
        path[str(0) + str(i)] = (str(0) + str(i - 1) + 'l')
    for j in range(1, len(word2)):
        al_frame.iloc[j, 0] = al_frame.iloc[j - 1, 0] - 1
        path[str(j) + str(0)] = (str(j - 1) + str(0) + 'v')

    for i in range(1, len(word1)):
        for j in range(1, len(word2)):
            a = al_frame.iloc[j, i - 1]
            b = al_frame.iloc[j - 1, i]
            c = al_frame.iloc[j - 1, i - 1]
            if word1[i] == word2[j]:
                al_frame.iloc[j, i] = max(a - 1, b - 1, c + 1)
                if c + 1 == max(a - 1, b - 1, c + 1):
                    path[str(j) + str(i)] = (str(j - 1) + str(i - 1) + 'd')
                elif b - 1 == max(a - 1, b - 1, c + 1):
                    path[str(j) + str(i)] = (str(j - 1) + str(i) + 'v')
                else:
                    path[str(j) + str(i)] = (str(j) + str(i - 1) + 'l')
            else:
                al_frame.iloc[j, i] = max(a, b, c) - 1
                if c - 1 == max(a - 1, b - 1, c - 1):
                    path[str(j) + str(i)] = (str(j - 1) + str(i - 1) + 'd')
                elif b - 1 == max(a - 1, b - 1, c - 1):
                    path[str(j) + str(i)] = (str(j - 1) + str(i) + 'v')
                else:
                    path[str(j) + str(i)] = (str(j) + str(i - 1) + 'l')

    point = str(j) + str(i)

    while path[point] != '00':
        path_point = path[point][2]
        j = point[0]
        i = point[1]
        if path_point == 'd':
            result_al1 = word1[int(i)] + result_al1
            result_al2 = word2[int(j)] + result_al2
        elif path_point == 'v':
            result_al1 = '_' + result_al1
            result_al2 = word2[int(j)] + result_al2
        else:
            result_al1 = word1[int(i)] + result_al1
            result_al2 = '_' + result_al2
        point = path[point][:2]
    return al_frame, result_al1, result_al2





# Основная функция

# ты говорил, что операции можно вызывать с помощью словаря, но я не совсем поняла, как с его помощью компактно обращаться к другим функциям
def AK_functions(*args):
    *seqs, operation = args

    full_result = []
    
    i=0
    isAKseq = is_ak(seqs)
    if (operation == 'Needlman_Wunsch' and len(seqs) >=2):
        for seq1 in seqs:
            for seq2 in seqs:
                if seq1 != seq2:
                    res = Needlman_Wunsch(seq1, seq2)
                    full_result.append(res)
        return full_result
    
    for seq in seqs:
        if isAKseq[i] == ('no, don`t do that'):
            result = (seq + "- isn't sequence")
        else:
            if operation == 'DNA_karta':
                result = DNA_karta(seq)
            elif operation == 'palindrom':
                result = palindrom(seq)
            elif operation == 'letter_pie':
                result = letter_pie(seq)
        full_result.append(result)
        i += 1

    return full_result