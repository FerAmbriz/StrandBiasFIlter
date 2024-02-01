import pandas as pd
import io
from math import factorial
import warnings
import sys
import numpy as np
import concurrent.futures
from multiprocessing import Pool, cpu_count

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

#======================================Inputs==============================#

input_caller = sys.argv[1]
caller = sys.argv[2]
output = sys.argv[3]
out_ID = sys.argv[4]

#=====================================Functions=============================#
def read_annot_txt(path):
    with open(path) as file:
        x = file.read(); x = x.split('\n')
    lst = []
    for line in x:
        line = line.split('\t');lst.append(line)
    df = pd.DataFrame(lst); df.columns = df.iloc[0]; df = df[1:]

    # Obtener el número de columnas que tienen nombre None
    num_none_cols = df.columns.isnull().sum()

    # Crear una lista con los nuevos nombres de las columnas
    new_columns = df.columns.tolist()
    other_info_index = 1
    for i, col in enumerate(new_columns):
        if col is None:
            new_columns[i] = f'Otherinfo.{other_info_index}'
            other_info_index += 1

    # Asignar los nuevos nombres de las columnas al DataFrame
    df.columns = new_columns; df = df.dropna(subset=['Start'])
    return df

def merge_pos(chrom, start):
    return chrom + ':' + str(start)

def exact_fisher(df):
    lst_fs = []
    for i in range(0, len(df)):
        Cm = df['Rwt'][i]
        Cp = df['Fwt'][i]
        Vp = df['Fmut'][i]
        Vm = df['Rmut'][i]

        FS = (max(Vp*Cm, Vm*Cp))/((Vp*Cm)+(Vm*Cp))
        lst_fs.append(FS)
    df['FS_manual'] = lst_fs
    return df

def SB_calc(df):
    lst_sb = []
    for i in range(0, len(df)):
        a = df['Fwt'][i]
        b = df['Fmut'][i]
        c = df['Rwt'][i]
        d = df['Rmut'][i]

        SB = (abs((a/(a+b))-(d/(c+d))))/((b+d)/(a+b+c+d))
        lst_sb.append(SB)
    df['SB_manual'] = lst_sb
    return df

def gatk_sb(df):
    lst_gatk = []
    for i in range(0, len(df)):
        a = df['Fwt'][i]
        b = df['Fmut'][i]
        c = df['Rwt'][i]
        d = df['Rmut'][i]

        gatk_sb = max((((b/(a+b))*(c/(c+d)))/((a+c)/(a+b+c+d))),(((d/(c+d))*(a/(a+b)))/((a+c)/(a+b+c+d))))
        lst_gatk.append(gatk_sb)
    df['gatkSB_manual'] = lst_gatk
    return df

def pFS(df):
    lst_pFS = []
    df = df.fillna({'Fwt': 0, 'Fmut': 0, 'Rmut':0, 'Rwt':0})
 
    for i in range(0, len(df)):
        a = int(df['Fwt'][i])
        b = int(df['Fmut'][i])
        c = int(df['Rwt'][i])
        d = int(df['Rmut'][i])
        n = a+b+c+d

        pFS = ((factorial(a+b)) * (factorial(c+d)) * (factorial(a+c)) * (factorial(b+d))) / ((factorial(a))*(factorial(b))*(factorial(c))*(factorial(d))*(factorial(n)))
        pFS = 1-pFS
        lst_pFS.append(pFS)
    df['pFS_manual'] = lst_pFS
    return df

def all_strandbias(df):
    df = exact_fisher(df)
    df = SB_calc(df)
    df = gatk_sb(df)
    df = pFS(df)
    return df

def process_row_strelka(row):
	print(row)
	# Create a dictionary to store the values for the new row
	new_row = {}
	# Iterate over the values in column 'a' and 'b'
	for col, val in zip(row['Otherinfo.11'], row['Otherinfo.12']):
		# Add the value to the new row
		new_row[col] = val
	row_dict = row.to_dict()

	new_rowglobal = {}
	new_rowglobal.update(new_row); new_rowglobal.update(row_dict)
	return new_rowglobal

def process_row_haplotype(row):
    print(row['Otherinfo11'])
    # Create a dictionary to store the values for the new row
    col = []; val = []

    for i in row['Otherinfo11']:
        print(i)
        [x,y] = i.split('=')
        col.append(x); val.append(y)
    new_row = dict(zip(col, val))
    row_dict = row.to_dict()

    new_rowglobal = {}
    new_rowglobal.update(new_row); new_rowglobal.update(row_dict)
    return new_rowglobal

# Define a custom function to convert values to float
def convert_to_float(x):
	# Split the value on '.'
	parts = x.split('.')
	# If there are more than two parts
	if len(parts) > 2:
		# Only keep the first two parts
		x = f'{parts[0]}.{parts[1]}'
	# Convert the value to float
	return float(x)

def extract_strand(df, caller):
	df['merge_pos'] = list(map(merge_pos, df['Chr'], df['Start']))
	if caller == 'strelka':
		df[['Otherinfo.11', 'Otherinfo.12']] = df[['Otherinfo.11', 'Otherinfo.12']].apply(lambda x: x.str.split(':'))
		columns = pd.unique(df['Otherinfo.11'].explode())
		df_new = pd.DataFrame(columns=columns)
		# Create a pool of workers
		pool = Pool(cpu_count())
		results = pool.map(process_row_strelka, [row for _, row in df.iterrows()])
		df_bd  = df_new.append(results, ignore_index=True)

		df_bd['SB'] = list(map(float, df_bd['SB']))
		df_bd[['Fwt', 'Fmut']] = df_bd['ADF'].str.split(',', 1, expand=True)
		df_bd[['Rwt', 'Rmut']] = df_bd['ADR'].str.split(',', 1, expand=True)
		df_bd['Fwt'] = df_bd['Fwt'].astype(float)
		df_bd['Fmut'] = df_bd['Fmut'].str.replace(',', '.')
		#df_bd['Fmut'] = df_bd['Fmut'].apply(convert_to_float)
		df_bd['Fmut'] = df_bd['Fmut'].astype(float)
		df_bd['Rwt'] = df_bd['Rwt'].astype(float)
		df_bd['Rmut'] = df_bd['Rmut'].str.replace(',', '.')
		#df_bd['Rmut'] = df_bd['Rmut'].apply(convert_to_float)
		df_bd['Rmut'] = df_bd['Rmut'].astype(float)

	elif caller == 'haplotype':
		df[['Otherinfo11']] = df[['Otherinfo11']].apply(lambda x: x.str.split(';'))
		# Create a pool of workers
		pool = Pool(cpu_count())
		results = pool.map(process_row_haplotype, [row for _, row in df.iterrows()])
		df_bd = pd.DataFrame(results)

		df_bd['FS'] = list(map(float, df_bd['FS']))
		df_bd['DP'] = list(map(float, df_bd['DP']))
	elif caller == 'varscan':
		df[['Otherinfo.11', 'Otherinfo.12']] = df[['Otherinfo.11', 'Otherinfo.12']].apply(lambda x: x.str.split(':'))
		columns = pd.unique(df['Otherinfo.11'].explode())
		df_new = pd.DataFrame(columns=columns)
		# Create a pool of workers
		pool = Pool(cpu_count())
		results = pool.map(process_row_strelka, [row for _, row in df.iterrows()])
		df_bd  = df_new.append(results, ignore_index=True)
		df_bd.rename(columns={'RDF': 'Fwt', 'RDR':'Rwt', 'ADF': 'Fmut', 'ADR':'Rmut'}, inplace=True)
		df_bd['Fwt'] = df_bd['Fwt'].astype(float)
		df_bd['Rwt'] = df_bd['Rwt'].astype(float)
		df_bd['Fmut'] = df_bd['Fmut'].astype(float)
		df_bd['Rmut'] = df_bd['Rmut'].astype(float)

	return df_bd

def filtered_df(df, caller, output, out_ID):
    if caller == 'strelka':
        df = extract_strand(df, caller); print(df)
        df = all_strandbias(df)
        #df = df[df['FT']=='PASS']
        df = df[df['pFS_manual']<0.99]
        df = df[df['FS_manual']<0.99]
	# modificar a float si es numero y dejar como string '.'
        df['DP'] = df['DP'].fillna(0).astype(int)
        df['DPI'] = df['DPI'].fillna(0).astype(int)
        df['AllDepths'] = df['DP'] + df['DPI']
        df = df[df['AllDepths'] > 30]

    elif caller == 'haplotype':
        df = extract_strand(df, caller); print(df)
        df = df[df['FS']<60]
        df = df[df['DP']>30]
    elif caller == 'varscan':
        df = extract_strand(df, caller); print(df)
        df = all_strandbias(df)	
        df = df[df['pFS_manual']<0.99]
        df = df[df['FS_manual']<0.99]
	
    df['esp6500siv2_all'] = [float(x) if isinstance(x, str) and x.replace('.', '', 1).isdigit() else x for x in df['esp6500siv2_all']]
    # modificar la condición para verificar primero si el elemento es una instancia de float antes de intentar compararlo con 0.001
    df = df[(df['esp6500siv2_all'].apply(lambda x: isinstance(x, float) and x <= 0.001)) | (df['esp6500siv2_all'] == '.')]
    df = df[(df['Func.refGene']=='exonic')|(df['Func.refGene']=='splicing')]
    df = df[df['ExonicFunc.refGene']!='synonymous SNV']

    df_pathogenic = df[df['CLNSIG'].str.contains('pathogenic', case=False, na=False)]
    print('Pathogenic: ', df_pathogenic)
    df_non_pathogenic = df[~df['merge_pos'].isin(list(df_pathogenic['merge_pos']))]
    print('NonPathogenic: ', df_non_pathogenic)

    df_pathogenic.to_csv(output + '/StrandBiasFiltered_'+ caller + out_ID +'_pathogenic.txt', index = False, sep = '\t')
    df_non_pathogenic.to_csv(output + '/StrandBiasFiltered_' + caller + out_ID +'_non_pathogenic.txt', index = False, sep = '\t')
    return df


#======================================Run==================================#
if caller == 'strelka':
	df_strelka = read_annot_txt(input_caller); print(df_strelka)
	df_strelka_filtered = filtered_df(df_strelka, 'strelka', output, out_ID)
	print(df_strelka_filtered)
elif caller == 'haplotype':
	df_haplotype = read_annot_txt(input_caller); print(df_haplotype)
	df_haplotype_filtered = filtered_df(df_haplotype, 'haplotype', output, out_ID)
elif caller == 'varscan':
	df_varscan = read_annot_txt(input_caller); print(df_varscan)
	df_varscan_filtered = filtered_df(df_varscan, 'varscan', output, out_ID)

