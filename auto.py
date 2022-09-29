import seepp_input
import seepp_label
import enc_input
import enc_label
import status_check
from utils import *

import sEEPP_data_generation
import error_check

# rename_annotation_result()
status_check.status_update()

df = pd.read_csv(status_table)
for i in range(len(df)):
    if df['sEEPP'].iloc[i] and df['genome download done'].iloc[i] and df['Cufflink result download done'].iloc[i]:
        if df['sEEPP dataset generation'].iloc[i] != 1:
            sEEPP_data_generation.set_seed()
            work = sEEPP_data_generation.GenerateEEPPData(df['genome_id'].iloc[i])
            work.generate_and_save_data()
            del work
            df.loc[i, 'sEEPP dataset generation'] = 1
        if df['sEEPP dataset generation'].iloc[i] == 1 and df['sEEPP error'].iloc[i] != '0':
            df.loc[i, 'sEEPP error'] = error_check.error_check(df['genome_id'].iloc[i], 'sEEPP')
        if df['sEEPP dataset generation'].iloc[i] == 1 and int(df['sEEPP error'].iloc[i]) == 0:
            df.loc[i, 'sEEPP dataset available'] = 1

df.to_csv(status_table, index=False)

    # if df['error'].iloc[i] != 1 and df['sEEPP'].iloc[i] == 1 and df['genome download done'].iloc[i] == 1:
    #     if df['sEEPP input generation done'].iloc[i] == 0:
    #         seepp_input.generate_seepp_input(df['genome_id'].iloc[i])
    #     if df['Cufflink result download done'].iloc[i] == 1 and df['sEEPP label generation done'].iloc[i] == 0:
    #         seepp_label.generate_seepp_label(df['genome_id'].iloc[i])
    #
    # # if df['error'].iloc[i] != 1 and df['ENC'].iloc[i] == 1 and df['genome download done'].iloc[i] == 1:
    # if df['error'].iloc[i] != 1 and df['genome download done'].iloc[i] == 1:
    #     if df['ENC input generation done'].iloc[i] == 0:
    #         enc_input.generate_enc_input(df['genome_id'].iloc[i])
    #     if df['BRENDA result download done'].iloc[i] == 1 and df['ENC label generation done'].iloc[i] == 0:
    #         enc_label.generate_enc_label(df['genome_id'].iloc[i])

status_check.status_update()
