import seepp_input
import seepp_label
import enc_input
import enc_label
import status_check
from utils import *

rename_annotation_result()
status_check.status_update()

df = pd.read_csv(status_table)
for i in range(len(df)):
    if df['error'].iloc[i] != 1 and df['sEEPP'].iloc[i] == 1 and df['genome download done'].iloc[i] == 1:
        if df['sEEPP input generation done'].iloc[i] == 0:
            seepp_input.generate_seepp_input(df['genome_id'].iloc[i])
        if df['Cufflink result download done'].iloc[i] == 1 and df['sEEPP label generation done'].iloc[i] == 0:
            seepp_label.generate_seepp_label(df['genome_id'].iloc[i])

    # if df['error'].iloc[i] != 1 and df['ENC'].iloc[i] == 1 and df['genome download done'].iloc[i] == 1:
    if df['error'].iloc[i] != 1 and df['genome download done'].iloc[i] == 1:
        if df['ENC input generation done'].iloc[i] == 0:
            enc_input.generate_enc_input(df['genome_id'].iloc[i])
        if df['BRENDA result download done'].iloc[i] == 1 and df['ENC label generation done'].iloc[i] == 0:
            enc_label.generate_enc_label(df['genome_id'].iloc[i])

status_check.status_update()
