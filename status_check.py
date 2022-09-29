from utils import *


# def data_check(gid, genome_name, type, model):
#     if model == 'sEEPP':


def status_update():
    print('updating data status..')
    df = pd.read_csv(status_table)

    for i in range(len(df)):
        genome_name = df.genome_name.iloc[i]
        gid = df.genome_id.iloc[i]

        df.loc[i, 'genome download done'] = int(os.path.exists(os.path.join(root, '전처리 이전 데이터', 'gene_input_raw', f'Merge.{gid}.csv')))
        df.loc[i, 'BRENDA result download done'] = int(os.path.exists(os.path.join(root, '전처리 이전 데이터', 'ENC_label_raw', f'{genome_name}.csv')))
        df.loc[i, 'Cufflink result download done'] = int(os.path.exists(os.path.join(root, '전처리 이전 데이터', 'sEEPP_label_raw', f'Total_{gid}.csv')))

        # df.loc[i, 'ENC input generation done'] = int(os.path.exists(os.path.join(root, 'data', 'ENC', genome_name)))
        # df.loc[i, 'ENC label generation done'] = int(os.path.exists(os.path.join(root, 'label', 'ENC', genome_name)))
        # df.loc[i, 'ENC dataset available'] = int(df['ENC label generation done'].iloc[i] and df['ENC input generation done'].iloc[i])

        # df.loc[i, 'sEEPP input generation done'] = int(os.path.exists(os.path.join(root, 'data', 'sEEPP', genome_name)))
        # df.loc[i, 'sEEPP label generation done'] = int(os.path.exists(os.path.join(root, 'label', 'sEEPP', genome_name)))
        # df.loc[i, 'sEEPP dataset available'] = int(df['sEEPP label generation done'].iloc[i] and df['sEEPP input generation done'].iloc[i])

    df.to_csv(status_table, index=False)


if __name__ == '__main__':
    status_update()