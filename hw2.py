import pandas as pd

print('Hello, Ivan!')
AK = input('coise')
gene_map = pd.read_csv('ЯКарта.csv')

# преобразование файла, нужно будет сохранить и использовать один и тот же, названия столбцов поменять
gene_map.rename(columns = {'Unnamed: 2':'AK'}, inplace = True )
for i in range(len(gene_map.Кодон)):
    for j in range(len(gene_map.Кодон)):
        if gene_map.Аминокислота[i] == gene_map['Аминокислота.1'][j]:
            gene_map.AK[i] = gene_map.Сокращение[j]

seq = 'fli'
SEQ = seq.upper()


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


DNA_karta(SEQ)

