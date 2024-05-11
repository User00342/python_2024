import pandas as pd
import matplotlib.pyplot as plt


print('')

gene_map = pd.read_csv('ЯКарта.csv')

# преобразование файла, нужно будет сохранить и использовать один и тот же, названия столбцов поменять
gene_map.rename(columns = {'Unnamed: 2':'AK'}, inplace = True )
for i in range(len(gene_map.Кодон)):
    for j in range(len(gene_map.Кодон)):
        if gene_map.Аминокислота[i] == gene_map['Аминокислота.1'][j]:
            gene_map.AK[i] = gene_map.Сокращение[j]



# Проверка и подготовка:

def is_ak(*seqs, alphabet_AK = set('ARNDCQEGHILKMFPSTWYV')):
    answer = []
    for seq in seqs:
        print(seq)
        seq = set(seq.upper())
        if seq <= alphabet_AK:
            answer.append('jup')
        else:
            answer.append('no, don`t do that')
    return answer




# восстановление последовательности
def DNA_karta(seq):
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

    # составление словарика с количеством аминокислот в последовательности
    AK_counts_dict = {}
    for letter in SEQ:
        if letter in AK_counts_dict:
            AK_counts_dict[letter] += 1
        else:
            AK_counts_dict[letter] = 1

    print('В представленной последовательности', len(AK_counts_dict), 'видов аминокислот.')

    # мини-визуализация
    plt.pie(list(AK_counts_dict.values()), labels=list(AK_counts_dict.keys()), autopct='%1.1f%%', shadow=True,
            wedgeprops={'width': 0.2, 'lw': 1, 'ls': '--', 'edgecolor': "k"})
    plt.show()
    return AK_counts_dict

def palindrom(seq):
    # функция проверяет, является ли последовательность палиндромом
    ka = seq[:len(seq)//2:-1]
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

    ak1 = ak1.lower
    ak2 = ak2.lower
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

letter_pie('asda')
print(is_ak('as3s','aaa'))





# Основная функция

