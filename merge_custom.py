import sys
import glob
import pandas as pd

folder = sys.argv[1]


# Patrón de archivo
pattern = '*pathogenic.txt'
# Lista de archivos que coinciden con el patrón
files = glob.glob(f'{folder}/{pattern}')
# Filtrar la lista para excluir archivos 'non_pathogenic'
files = [file for file in files if 'non_pathogenic' not in file]

df_bd = pd.DataFrame()

for i in files:
    df_i = pd.read_csv(i, sep='\t')
    df_bd = pd.concat([df_bd, df_i])

print(df_bd)


df_bd.to_csv('StrandBiasFiltered_haplotype_merge_pathogenic.txt', sep='\t', index = False)


# Patrón de archivo
pattern = '*non_pathogenic.txt'
# Lista de archivos que coinciden con el patrón
files = glob.glob(f'{folder}/{pattern}')

df_bd = pd.DataFrame()

for i in files:
    df_i = pd.read_csv(i, sep='\t')
    df_bd = pd.concat([df_bd, df_i])

print(df_bd)


df_bd.to_csv('StrandBiasFiltered_haplotype_merge_non_pathogenic.txt', sep='\t', index = False)
