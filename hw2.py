import pandas as pd
import matplotlib.pyplot as plt


# подготовка материала
gene_map = pd.read_csv('ЯНормальнаяКарта.csv')

def is_ak(seqs):

    ''' 
    Проверка являются ли последовательности аминокислотами. Кавычки теперь работают.
    '''

    alphabet_AK = set('ARNDCQEGHILKMFPSTWYV')
    answer = [] # список результатов работы функции
    for seq in seqs:
        seq = set(seq.upper())
        if seq <= alphabet_AK:
            answer.append('jup')
        else:
            answer.append('no, don`t do that')
    return answer



# вспомогательные функции

def DNA_karta(seq): 
    ''' 
    функция восстановления последовательности, подбор триплетов
    ''' 

    all_variants = 1 # количество возможных вариантов получения заданной ак-последовательности
    DNA_seq = ''
    map_seq = seq.upper()
    for i in map_seq:
        variants = list(gene_map.Кодон[gene_map.AK == i]) # gene_map - датафрейм, генетическая карта; variants - возможные триплеты
        all_variants *= len(variants)
        DNA_seq = DNA_seq + variants[0]
    return all_variants, DNA_seq


def letter_pie(SEQ, plot):
    '''  
    подсчет количества аминокислот в заданной последователсти seq и визуализация результата
    ''' 

    AK_counts_dict = {} # AK_counts_dict - словарь с указанием количества каждой аминокислоты в последовательности
    for letter in SEQ:
        if letter in AK_counts_dict:
            AK_counts_dict[letter] += 1
        else:
            AK_counts_dict[letter] = 1

    message = ('В последовательности ' + SEQ + ' - ' + str(len(AK_counts_dict)) + ' видов аминокислот.')

# мини-визуализация

    if plot:
        plt.show(plt.pie(list(AK_counts_dict.values()), labels=list(AK_counts_dict.keys()), autopct='%1.1f%%', shadow=True, 
                         wedgeprops={'width': 0.2, 'lw': 1, 'ls': '--', 'edgecolor': "k"}))

    return message

def is_palindrom(seq):
    ''' 
    является ли последовательность палиндромом
    '''
    ka = seq[:(len(seq)-1)//2:-1]
    ak = seq[:len(seq)//2]
    if ak == ka:
        result = seq + ' - palindrom'
    else:
        result = seq + " - isn't palindrom"
    return result

def Needlman_Wunsch(ak1, ak2):
    '''
    Воссоздание знаменитого алгоритма Нидлмана-Вунша, его нечитабельная версия.

    '''

    ak1 = ak1.lower()
    ak2 = ak2.lower()
    word1 = ['']
    word2 = ['']
    path = {(0,0): (0,0)} # словарь пути алгоритма, все еще без кортежей, этот момент я буду откладывать до последнего
    result_al1 = ''
    result_al2 = ''
     
     # разбиение каждой последовательности
    for letter in ak1:
        word1.append(letter)
    for letter in ak2:
        word2.append(letter)

    # создание и заполнение первых строк датафрейма. Запись пути заполнения в словарь path
    al_frame = pd.DataFrame(columns=word1, index=word2) #  матрица схожести
    al_frame.iloc[0][0] = 0
    for i in range(1, len(word1)):
        al_frame.iloc[0, i] = al_frame.iloc[0, i - 1] - 1
        path[(0,i)] = ((0, i - 1, 'l'))
    for j in range(1, len(word2)):
        al_frame.iloc[j, 0] = al_frame.iloc[j - 1, 0] - 1
        path[(j, 0)] = ((j - 1, 0,'v'))
        
    # заполнение датафрейма по алгоритму Нидлмана-Вунша. Запись пути заполнения в словарь path
    for i in range(1, len(word1)):
        for j in range(1, len(word2)):
            a = al_frame.iloc[j, i - 1]
            b = al_frame.iloc[j - 1, i]
            c = al_frame.iloc[j - 1, i - 1]
            if word1[i] == word2[j]:
                al_frame.iloc[j, i] = max(a - 1, b - 1, c + 1)
                if c + 1 == max(a - 1, b - 1, c + 1):
                    path[(j,i)] = ((j - 1,i - 1,'d'))
                elif b - 1 == max(a - 1, b - 1, c + 1):
                    path[(j,i)] = ((j - 1,i,'v'))
                else:
                    path[(j,i)] = ((j,i - 1,'l'))
            else:
                al_frame.iloc[j, i] = max(a, b, c) - 1
                if c - 1 == max(a - 1, b - 1, c - 1):
                    path[(j,i)] = ((j - 1,i - 1,'d'))
                elif b - 1 == max(a - 1, b - 1, c - 1):
                    path[(j,i)] = ((j - 1,i,'v'))
                else:
                    path[(j,i)] = ((j,i - 1,'l'))

    point = (j,i)
#    return al_frame, path

    # Обратный ход по пути, запись выровненных последовательностей
    while path[point] != (0,0):
        path_point = path[point][-1]
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

def break_all(seq):
    '''
    Наш план - помешать тайной корпорации в создании генномодифицированного кузнечика. Злодеи целого мира сплотились, чтобы получить белок Х, 
    который способен поработить прыгунов!
    В нашем арсенале не так много рестриктаз... 
    Способны ли мы остановить негодяев...? Настал час узнать ответ на этот вопрос...
    '''

    DNA_seq = DNA_karta(seq)[-1]
    arsenal = {'restrict1': 'GAACT',
               'restrict2': 'CAGT',
               'restrict3': 'TTT',
               'restrict4': 'GTGCTC'}
    sites_places = {}
    for site in arsenal.values():
        site_length = len(site)
        k = 1
        for i in range(len(DNA_seq) - site_length + 1):
            if DNA_seq[i:i+site_length] == site:
                if arsenal.get(site) == '':
                    sites_places[site] = i
                else:
                    sites_places[str(k) + '-й ' + site] = i
                k += 1
    if sites_places == {}:
        sites_places = 'нам стоило подготовиться получше... '
                
                
    return sites_places





# Основная функция


def AK_functions(*args):
    '''
    DNA_karta - восстановление последовательности ДНК
    is_palindrom - проверка является ли палиндромом
    letter_pie - подсчет количества различных аминокислот
    break_all - поиск заданных сайтов рестрикции
    Needlman_Wunsch - парное выравнивание
    '''

    *seqs, operation = args


    full_result = [] # результат работы всех вызванных вспомогательных функций
    result = '' # результат работы одной вспомогательной функции
    work_seqs = []
    func_names = ['DNA_karta', 'palindrom', 'letter_pie', 'restriction', 'Needlman_Wunsch']
    functions = [DNA_karta, is_palindrom, letter_pie, break_all, Needlman_Wunsch]
    functions_dict = dict(zip(func_names,functions))
    isAKseq = is_ak(seqs) # список являются ли последовательности аминокислотными
    
    for k in range(len(seqs)):
        if isAKseq[k] == 'jup':
            work_seqs.append(seqs[k])
        else: 
            full_result.append(seqs[k] + " - isn't sequence")

    if operation == 'Needlman_Wunsch':  
        if len (work_seqs) < 2: 
            full_result.append('Недостаточно последовательностей')
            return full_result    
        else: 
            for seq1 in work_seqs:
                work_seqs.remove(seq1)
                for seq2 in work_seqs:
                    if seq1 != seq2:
                        result = Needlman_Wunsch(seq1, seq2)
                        full_result.append(result)
            return full_result

    
    for seq in work_seqs:
        try: result = functions_dict[operation](seq)
        except Exception: result = 'Данная операция не найдена. Возможные операции: DNA_karta, palindrom, letter_pie, restriction, Needlman_Wunsch'

        if operation == 'DNA_karta':
            result = ('аминокислотная последовательность ' + seq + ' могла получиться из ' + str(result[0]) + " различных нуклеотидных последовательностей. " + "Пример нуклеотидной последовательности: " + str(result[1]))

        full_result.append(result)



    return full_result


AK_functions('aaa','adf')
