# Processing of amino acid sequences

*Amino acids* are the monomers that make up proteins.  
There are 23 amino acids overall, 20 of them are the most studied.
The **AK_functions** library allows you to study the primary sequence of a protein in more detail through a set of functions.

## Import

To run the `AK_functions()`, first import it as a module

```
import hw2
```
and run its main function, `AK_functions` function. This function provides interface for all 5 operations from `functions_dict` dictionary. 

- First `n` arguments - protein sequences
- last one - operation from list: "palindrom", "letter_pie", "break_all", "DNA_karta", "Needlman_Wunsch"


## Features of the function
### is_palindrom 
A palindrome is a sequence of symbols that reads the same backwards as forwards. This function takes an amino acid sequence as an argument, checks whether it is a palindrome, and returns the result.
### letter_pie 
Quantitative amino acid sequence analysis. It is possible to build a pie chart by passing `plot = True` to the function.
### break_all
Function to search for specific restriction sites in the reconstructed DNA sequence. The DNA sequence is obtained using the `DNA_karta` function.
### DNA_karta 
Triplets search function. The degeneracy of the genetic code refers to the fact that most amino acids are specified by more than one codon. The function also counts the number of possible DNA sequence variants.
### Needlman_Wunsch 
The Needleman–Wunsch algorithm is an algorithm used in bioinformatics to align proteins or nucleotide sequences. Takes a list of amino acid sequences as an argument and performs pairwise alignment. Function returns a matrix of the obtained values ​​and aligned sequences.

## Examples

**AK_functions('aararlss', 'afa', 'DNA_karta')**

```
['аминокислотная последовательность aararlss могла получиться из 497664 различных нуклеотидных последовательностей. Пример нуклеотидной последовательности: GCTGCTCGTGCTCGTCTTTCTTCT',
 'аминокислотная последовательность afa могла получиться из 32 различных нуклеотидных последовательностей. Пример нуклеотидной последовательности: GCTTTTGCT']
```
**AK_functions('aararlss', 'restriction')**

```
[{'1-й TTT': 16, '1-й GTGCTC': 7}]

```
**Needlman_Wunsch('aasdfgs', 'adasfkls')**

![изображение](https://github.com/MaslovaIrina/python_2024/assets/114800146/7142762f-39e2-4aeb-8909-5d490fc69d5e)

```
'a_asdf_gs',
'adas_fkls'
```
