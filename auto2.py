import status_check
from utils import *

import data_generation_2
import error_check

# rename_annotation_result()
# status_check.status_update()

df = pd.read_csv(status_table)
print('--------------------sEEPP--------------------')
for i in range(len(df)):
    if df['sEEPP dataset available'].iloc[i] != 1:
        if df['sEEPP'].iloc[i] and df['genome download done'].iloc[i] and df['Cufflink result download done'].iloc[i]:
            if df['sEEPP dataset generation'].iloc[i] != 1:
                set_seed()
                work = data_generation_2.GenerateEEPPData2(df['genome_id'].iloc[i])
                work.generate_and_save_data()
                del work
                df.loc[i, 'sEEPP dataset generation'] = 1
            if df['sEEPP dataset generation'].iloc[i] == 1 and df['sEEPP error'].iloc[i] != '0':
                df.loc[i, 'sEEPP error'] = error_check.sEEPP_error_check(df['genome_id'].iloc[i])
            if df['sEEPP dataset generation'].iloc[i] == 1 and int(df['sEEPP error'].iloc[i]) == 0:
                df.loc[i, 'sEEPP dataset available'] = 1
            df.to_csv(status_table, index=False)
print('---------------------ENC---------------------')
for i in range(len(df)):
    if df['ENC dataset available'].iloc[i] != 1:
        if df['ENC'].iloc[i] and df['genome download done'].iloc[i] and df['BRENDA result download done'].iloc[i]:
            if df['ENC dataset generation'].iloc[i] != 1:
                set_seed()
                work = data_generation_2.GenerateENCData2(df['genome_id'].iloc[i])
                work.generate_and_save_data()
                del work
                df.loc[i, 'ENC dataset generation'] = 1
            if df['ENC dataset generation'].iloc[i] == 1 and df['ENC error'].iloc[i] != '0':
                df.loc[i, 'ENC error'] = error_check.ENC_error_check(df['genome_id'].iloc[i])
            if df['ENC dataset generation'].iloc[i] == 1 and df['ENC error'].iloc[i] == 0:
                df.loc[i, 'ENC dataset available'] = 1
            df.to_csv(status_table, index=False)


# status_check.status_update()
