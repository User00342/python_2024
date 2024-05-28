def two_from_one(limit):
    '''
    создание нижней и верхней границы переменной
    '''   
    if type(limit) == int or type(limit) == float:
        limit = (0,limit)
    return limit


def del_def(seqs, del_names):
    '''
    создание нового словаря на основе старого
    '''
    result_dict = {}
    for name in seqs.keys():
        if name not in del_names:
            result_dict[name] = seqs[name]
    return result_dict

def check_nucl(seqs, del_names):
    nucleotides = set('ATGCatgc')
    for name in seqs.keys():
        if (set(seqs[name][0]) <= nucleotides) == False:
            del_names.append(name)
    return del_names
    
def check_equal_lengths(seqs, del_names):            
    for name in seqs.keys():    
        if len(seqs[name][0]) != len(seqs[name][1]):
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
    
    for name in seqs.keys():
        sred_value = 0
        for letter in seqs[name][1]:
            sred_value = sred_value + ord(letter) - 33
        quality = sred_value/len(seqs[name][1])
        if quality < quality_threshold:
            del_names.append(name)
    
    return del_names

def main(seqs, gc_bounds = (0,100), lenght_bounds = (0,2**32), quality_threshold = 0, additional_verifiaction = False):

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
    if additional_verifiaction:
        del_names = check_nucl(seqs, del_names)
        del_names = check_equal_lengths(seqs, del_names)
    del_names = check_gc(seqs, del_names, gc_bounds)
    del_names = check_length(seqs,del_names, lenght_bounds)
    del_names = check_quality(seqs, del_names, quality_threshold)

    result = del_def(seqs, del_names)
    return result
