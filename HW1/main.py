##### Вспомогательные функции #####


def kort_rev(seqs):
    stop_k = len(seqs)-1
    seqs = seqs[0:stop_k]
    return seqs

def rna_dna(seqs,alphabet_dna = set('AGCTagct'), alphabet_rna = set('AGCUagcu')):
    uniq_seq = set(seqs)
    if uniq_seq <= alphabet_dna:
        who = 'dna'
    elif uniq_seq <= alphabet_rna:
        who = 'rna'
    else:
        who = 'isnt seq'
    return who

def transcribe(seqs):
    seqs = seqs.replace('T', 'U')
    seqs = seqs.replace('t', 'u')
    return seqs

def reverse(seqs):
    return seqs[::-1]

def complement(seqs, who, alphabet_first = 'AaGgCcTt', alphabet_second = 'TtCcGgAa',alphabet_rna_first = 'AaGgCcUu', alphabet_rna_second = 'UuCcGgTt'):
    res_com_seq = ''
    if who == 'dna':
        for i in range(len(seqs)):
            res_com_seq += alphabet_second[alphabet_first.index(seqs[i])]
    else:
        for i in range(len(seqs)):
            res_com_seq += alphabet_rna_second[alphabet_rna_first.index(seqs[i])]
    return res_com_seq

def reverse_complement(seqs, who):
    res_seq = reverse(complement(seqs, who))
    return res_seq



##### Основная функция #####

def run_dna_rna_tools(*seqs):
    command = seqs[-1]
    seqs = kort_rev(seqs)
    result_seq = []

    for seq in seqs:
        uniq_seq = set(seq)
        who = rna_dna(uniq_seq)

        if (who == 'rna' or who == 'dna'):
            if command == 'transcribe':
                result_seq.append(transcribe(seq))
            elif command == 'reverse':
                result_seq.append(reverse(seq))
            elif command == 'complement':
                result_seq.append(complement(seq, who))
            elif command == 'reverse_complement':
                result_seq.append(reverse_complement(seq, who))
            else:
                print(command, ' - unknown command')


        else:
            print(seq, ' - isnt seq')
    if result_seq != []:
        print(result_seq)
        return result_seq



