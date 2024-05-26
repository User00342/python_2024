import pandas as pd


def two_from_one(limit):
    '''
    создание нижней и верхней границы переменной
    '''   
    if type(limit) == int or type(limit) == float:
        limit = (0,limit)
    return limit


def del_def(seqs, del_names):
    '''
    удаление из словаря строк по ключу
    '''
    for name in del_names:
        del seqs[name]
    return seqs

def add_verification(seqs, del_names, alphabet):    
    nucleotides = set(alphabet)
    
    for name in seqs.keys():    
        if len(seqs[name][0]) != len(seqs[name][1]):
            del_names.append(name)
        elif (set(seqs[name][0]) <= nucleotides) == False:
            del_names.append(name)
    return del_names


def check_gc(seqs, del_names, gc_bounds):
    gc_bounds = two_from_one(gc_bounds)  
    
    for name in seqs.keys():
        seq = seqs[name][0].upper()
        percent = 100*(seq.count('C') + seq.count('G'))/len(seq)
        if gc_bounds[0] > percent or percent > gc_bounds[1]:
            del_names.append(name)
            
    return del_names


def check_length(seqs, del_names, lenght_bounds):
    lenght_bounds = two_from_one(lenght_bounds)  

    for name in seqs.keys():  
        if lenght_bounds[0] > len(seqs[name][0]) or len(seqs[name][0]) > lenght_bounds[1]:
            del_names.append(name)  
    
    return del_names
 
def check_quality(seqs, del_names, quality_threshold):
    ASCII_df = pd.read_csv('ASCII.csv')
    ASCII_df.Dec = ASCII_df.Dec - 33
    
    for name in seqs.keys():
        sred_value = 0
        for letter in seqs[name][1]:
            sred_value = sred_value + int(ASCII_df.Dec[ASCII_df.Character == letter])
        quality = sred_value/len(seqs[name][1])
        if quality < quality_threshold:
            print(name)
            del_names.append(name)
    
    return del_names

def main(seqs, gc_bounds = (0,100), lenght_bounds = (0,2**32), quality_threshold = 0, add_verif = False):

    '''
    Фунция проверки fastq-последовательности по трем критериям:
    - процентное содержание gc нуклеотидов
    - длина последовательности
    - качество последовательности 

    Дополнительная возможность проверки: 
    - является ли последовательность нуклеотидной
    - одинаковой ли длины строки последовательности и ее качества
    --------------------------------
    Аргументы:
    seqs - словарь, состоящий из fastq-сиквенсов. Ключ - строка, имя последовательности. Значение - кортеж из двух строк: последовательность и качество
    gc_bounds - интервал GC состава (в процентах) для фильтрации (по-умолчанию равен (0, 100)). Если в аргумент передать одно число, то считается, что это верхняя граница.
    lenght_bounds - интервал длины для фильтрации, аналогично gc_bounds, по-умолчанию равен (0, 2**32).
    quality_threshold - пороговое значение среднего качества рида для фильтрации, по-умолчанию равно 0 (шкала phred33).
    add_verif - проводить ли дополнительную проверку, по умолчанию "False"
    --------------------------------
    Результат:
    Функция возвращает словарь, состоящих только из тех сиквенсов, которые прошли все условия.
    '''
    del_names = [] # список ключей строк, подлежащих удалению
    if add_verif:
        del_names = add_verification(seqs, del_names, alphabet='ATGCatgc')
    del_names = check_gc(seqs, del_names, gc_bounds)
    del_names = check_length(seqs,del_names, lenght_bounds)
    del_names = check_quality(seqs, del_names, quality_threshold)

    result = del_def(seqs, del_names)
    return result
