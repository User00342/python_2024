import datetime
from abc import ABCMeta, abstractmethod, ABC 
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
import statistics as stat
import numpy as np

class User():
    def __init__(self, name):
        self.name = name
        
class EquipmentError(Exception):
    pass 

class DateError(Exception):
    pass

class BookingError(Exception):
    pass

class Booking():
    '''
    хранение информации о настоящем бронировании?
    Мне кажется этот класс лучше смотрится с наследованием. Однако, как я понимаю из дисклеймера 3 - 
    его в этом задании использовать не нужно. Мне кажется что-то я в этом задании сделала не так, оно смотрится некрасиво, но учитывая требования из описания...
    Вероятно, завтра буду думать как его переделать
    
    date - дата брони 
    start_time - время начала брони
    end_time - время окончания брони
    equipment - наименование инструмента
    '''
    def __init__(self, device, name, date, start_time, end_time, schedule):
        self.equipment = device
        self.user = name
        self.date = date
        self.start_time = start_time
        self.end_time = end_time
        

        if device not in schedule.keys():
            raise EquipmentError('no such equipment')
        elif date not in schedule[device].keys():
            raise DateError('no such date')
            
    def is_intersect(self, other_booking):
        '''
        функция проверки свободного времени
        '''
        if not ((other_booking.start_time < self.start_time and other_booking.end_time <= self.start_time) or (other_booking.start_time >= self.end_time and other_booking.end_time > self.end_time)):
            return False
        return True
        
    
class LabEquipment():
    '''
    Класс добавления брони. Выполняет проверку возможности бронирования инвентаря в заданные дату и время. 
    
    schedule - расписание
    device - выбор инвентаря для бронирования
    date - выбор даты для бронирования
    
    
    '''
    
    def __init__(self, schedule):
        self.schedule = schedule
        
    def __new_cases(self, device, date):
        '''
        Функция добавления новых значений даты и инвентаря, если они отсутствуют
        '''
        
        if device not in self.schedule.keys():
            self.schedule[device] = {}
            print ('we have new device: ', device)
        if date not in self.schedule[device].keys():
            self.schedule[device][date] = {'start_time': [],
                                           'end_time':[],
                                           'user': []}
            print ('we have new date: ', date)
        
        return self
    
    def check_existence error_masseges(self, device, date):
        
        if device not in self.schedule.keys():
            raise EquipmentError('no such equipment')
        elif date not in self.schedule[device].keys():
            raise DateError('no such date')        
            
    def is_available(self, start_time, end_time, device, date):
        '''
        Проверка является ли время для брони свободным. В случае, если время занято, возращает имя забронировавшего человека.
        '''

        for i in range(len(self.schedule[device][date]['start_time'])):
            if not ((self.schedule[device][date]['start_time'][i] < start_time and self.schedule[device][date]['end_time'][i] <= start_time) or (self.schedule[device][date]['start_time'][i] >= end_time and self.schedule[device][date]['end_time'][i] > end_time)):
                busy_dude = self.schedule[device][date]['user'][i]
                return busy_dude
        return True
       
    def add_time(self, start_time, end_time, device, date, name):
        '''
        Добавление новой брони в расписание, если время свободно. 
        '''
        just_ones = self.is_available(start_time, end_time, device, date)
        if just_ones == True:
            self.schedule[device][date]['start_time'].append(start_time)
            self.schedule[device][date]['end_time'].append(end_time)
            self.schedule[device][date]['user'].append(name)
        elif just_ones == name:
            raise BookingError('Du hast das schon gemacht.')
        else:
            raise BookingError('you can speak about this with ' + just_ones)
        
        return self
      
    def book(self, start_time, end_time, device, date, name, add_new_device_or_date = False):
        '''
        Проверка добавлять ли новые значения инвентаря и времени, затем использование функции  добавления новой брони.
        
        add_new_device_or_date - по умолчанию False, новые дата и инвентарь не добавляются. В случае использования таковых будет возвращена ошибка.
        
        '''
        if add_new_device_or_date == True:
            self = self.__new_cases(device, date)
        else:
            self.error_masseges(device, date)
        
        self.add_time(start_time, end_time, device, date, name) 
  

    
class GenCodeInterpreter():

    def __init__(self):
        self.memory = [0]*5000
    
    def LetterError(Exception):
        pass
    def MemoryError(Exception):
        pass

    def eval(self, code):
        '''
        Принимает на вход строку с программой. При вызове `eval` происходит следующее:
        - Буфер очищается
        - Интерпретируется переданная программа
        - Результат возвращается пользователю
        
        Правила интерпретации следующие: 
        'A' - перемещение указателя вперед, 
        'T' - перемещение указателя назад, 
        'G' - инкремент значения в ячейке памяти, 
        'C' - декремент значения в ячейке памяти, 
        'N' - добавление значения в ячейке памяти в буфер на возврат пользователю
        '''
        bufer = ''
        our_place = 0
        errors = ''

        try:
            for letter in code:
                match letter:
                    case "A": our_place += 1
                    case "T": our_place -= 1
                    case "G": self.memory[our_place] += 1
                    case "C": self.memory[our_place] -= 1
                    case "N": bufer += str(chr(self.memory[our_place]))
                    case _ : raise LetterError('unknown command')
        except  IndexError: 
            raise MemoryError('Index out of memory')
        bufer = self.angl_rus(bufer)
        return bufer
    
    def angl_rus(self,bufer):
        '''
        Функция замены английского словаря на русский.
        '''
        rus_word = ''
        angl_alphabet = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ' # эта строка создана с поMощью чата гпт
        rus_alphabet = 'абцдефгхийклмнопкрстуфвсызАБЦДЕФГХИЙКЛМНОПКРСТУфВСЫЗ' # эта строка создана с помощью чата гпт

        for letter in bufer:
            try: 
                rus_word += rus_alphabet[angl_alphabet.index(letter)]
            except ValueError:
                rus_word += letter
        return rus_word



def meet_the_dunders():
    '''
    Повеселились, мне больше нечего здесь сказать.
    '''
    res = (0).__int__()
    matrix = [].__new__(list)

    for idx in range(0, 100, 10):
        matrix.__iadd__([list(range(idx, idx.__add__(10)))])

    def func_1(x):
        return x in range(1, 5, 2)
        
    def func_2(x):
        return [x[col] for col in selected_columns_indices]

        
    selected_columns_indices = list(filter(func_1.__call__, range(matrix.__len__())))
    selected_columns = map(func_2.__call__, matrix)

    arr = np.array(list(selected_columns))

    mask = arr[:, 1].__mod__(3) == 0
    new_arr = arr[mask]

    product = new_arr.__matmul__(new_arr.T)

    if (product[0].__lt__(1000)).all().__and__((product.__getitem__(2).__gt__(1000)).any()):
        res = ((product.mean().__floordiv__(10)).__mod__(100)).__int__()
    return res





def two_from_one(limit):
    '''
    создание нижней и верхней границы переменной
    '''   
    if type(limit) == int or type(limit) == float:
        limit = (0,limit)
    return limit 

def filter_fastq(file_name, gc_bounds = (0,100), lenght_bounds = (0,2**32), quality_threshold = (0,100)):
    '''
     Фунция проверки fastq-последовательности по трем критериям:
    - процентное содержание gc нуклеотидов
    - длина последовательности
    - качество последовательности 
    
    gc_bounds - интервал GC состава (в процентах) для фильтрации (по-умолчанию равен (0, 100)). Если в аргумент передать одно число, то считается, что это верхняя граница.
    lenght_bounds - интервал длины для фильтрации, аналогично gc_bounds, по-умолчанию равен (0, 2**32).
    quality_threshold - интервал среднего качества рида для фильтрации, по-умолчанию равен (0,100) (шкала phred33).
    
    '''
    good_sequences = [] 
    gc_bounds, lenght_bounds, quality_threshold = list(map(two_from_one, [gc_bounds, lenght_bounds, quality_threshold]))
    
    def check_lenght(record, lenght_bounds):
        return lenght_bounds[0] <=  len(record.seq) <= lenght_bounds[1]
    def check_gc(record, gc_bounds):
        return gc_bounds[0] <=  100 * gc_fraction(record.seq) <= gc_bounds[1]
    def check_quality(record, quality_threshold):
        return quality_threshold[0] <= stat.mean((record.letter_annotations["phred_quality"])) <= quality_threshold[1]
    
    with open('fastq.fastq') as file:
        for record in SeqIO.parse(file, "fastq"):
            if check_lenght(record, lenght_bounds) and check_gc(record, gc_bounds) and check_quality(record, quality_threshold):
                good_sequences.append(record)
    return "Found %i good sequences" % len(good_sequences)



class AlphabetError(Exception):
    pass
class ComplementError(Exception):
    pass

class BiologicalSequence(ABC):   #Вау, абстрактные методы и абстрактные классы это разные вещи. Этот класс существует только для теоретической халочки?
        
    @abstractmethod 
    def __len__(self):
        pass
    
    @abstractmethod
    def __getitem__(self):
        pass
    
    @abstractmethod
    def slice(self):
        pass
   
    @abstractmethod
    def __str__(self):
        pass
    
    @abstractmethod
    def __repr__(self):
        pass
    
    @abstractmethod
    def check_alphabet(self):
        pass

    
    
    
class NucleicAcidSequence(BiologicalSequence):
    '''
    Класс обработки нуклетидных последовательностей. Включает в себя две главные функции:
    
    - gc_perc - подсчет процентного содержания гуанина и цитозина
    - seqs_complement - создание комплементарной цепи
    - check_alphabet - функция проверки наличия введенных букв в алфавите
    - slice - функция среза

    Дааа, dunder методы не все работают так, как следовало бы, потому что список. Но зато результат возвращает удобно, иначе... Мы списались, я не дописала
    '''
    alphabet = ''
    complement_alphabet = ''
    def __init__(self, seqs):
        self.seqs = seqs
        if all(map(self.check_alphabet, self.seqs)) == True:
            return True
        raise AlphabetError('Isn,t in alphabet')
        
    def __len__(self):
        return list(map(len, self.seqs)) 
    
    def __getitem__(self, idx):
        return [i[idx] if len(i) > idx else ' ' for i in self.seqs ]
    
    def slice(self, from_letter, to_letter):
        return [seq[from_letter:to_letter] if to_letter < len(seq) else seq[from_letter:] for seq in self.seqs]
   
    def __str__(self):
        return str(self.seqs)
    
    def __repr__(self):
        return ' '.join(self.seqs)
    
    def check_alphabet(self, seq):
        return set(seq) <= set(self.alphabet)
        
    def gc_perc(self):
        return [100 * gc_fraction(seq) for seq in self.seqs]
        
    def seq_complement(self,seq):
        res_com_seq = ''
        for i in range(len(seq)):
            res_com_seq += self.complement_alphabet[self.alphabet.index(seq[i])]
        return res_com_seq
    
    def seqs_complement(self):
        return [self.seq_complement(seq) for seq in self.seqs]



class DNASequence(NucleicAcidSequence):
    '''
    Класс обработки последовательности ДНК. 
    Главная функция - transcribe (транскрибция цепи)
    '''
    alphabet = 'agctAGCT'
    complement_alphabet = 'tcgaTCGA'
    
    def __init__(self,seqs):
        super().__init__(seqs)
        
    def replace_trans(self,seq):
        seq = seq.replace('T', 'U')
        seq = seq.replace('t', 'u')
        return seq

    def transcribe(self):
        return [self.replace_trans(seq) for seq in self.seqs]
    
    
    
class RNASequence(NucleicAcidSequence):
    '''
    Класс обработки последовательности РНК
    '''
    alphabet = 'agcuAGCU'
    complement_alphabet = 'ucgaUCGA'
    
    def __init__(self,seqs):
        super().__init__(seqs)

  
class AminoAcidSequence(NucleicAcidSequence):
    '''
    Класс обработки аминокислотной последовательности. 
    '''
    alphabet = 'ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv'
    
    def seqs_complement(self):
        raise ComplementError('just for rna and dna')
        
    def __init__(self,seqs):
        super().__init__(seqs)
        
    def copy_palindrom(self,seq):
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
    
    def is_palindrom(self):
        return [self.copy_palindrom(seq) for seq in self.seqs]
  
