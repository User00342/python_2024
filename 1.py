class GenCodeInterpreter():
    def __init__(self):
        self.memory = [0]*15

    def eval(self, code):
        bufer = ''
        our_place = 0
        errors = ''
        print('init')
        for letter in code:
            if letter == "A": our_place += 1
            elif letter == "T": our_place -= 1
            elif letter == "G": self.memory[our_place] += 1
            elif letter == "C": self.memory[our_place] -= 1
            elif letter == "N": bufer += str(chr(self.memory[our_place]))
            else: errors = errors + ' ' + letter + ' -  unknown command. '
        bufer = self.angl_rus(bufer)
        return bufer, errors
    
    def angl_rus(self,bufer):
        rus_word = ''
        angl_alphabet = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ' # эта строка создана с посощью чата гпт
        rus_alphabet = 'абцдефгхийклмнопкрстуфвсызАБЦДЕФГХИЙКЛМНОПКРСТУфВСЫЗ' # эта строка создана с помощью чата гпт

        for letter in bufer:
            rus_word += rus_alphabet[angl_alphabet.index(letter)]
        return rus_word




interpreter = GenCodeInterpreter()
# code = 'C'*3
code = 'G'*68 + 'NA' + 'G'*105 + 'NA' + 'G'*111 +'NA' + 'G'*103 + 'NA' + 'G'*101 + 'NA' + 'G' * 110 + 'N'
print(interpreter.eval(code))
print(interpreter.memory)

