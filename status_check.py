from utils import *


def status_update():
    print('updating data status..')
    df = pd.read_csv(status_table)

    for i in range(len(df)):
        genome_name = df.genome_name.iloc[i]
        gid = df.genome_id.iloc[i]

        df.loc[i, 'genome download done'] = int(os.path.exists(os.path.join(root, '전처리 이전 데이터', 'gene_input_raw', f'Merge.{gid}.csv')))
        df.loc[i, 'BRENDA result download done'] = int(os.path.exists(os.path.join(root, '전처리 이전 데이터', 'ENC_label_raw', f'{genome_name}.csv')))
        df.loc[i, 'Cufflink result download done'] = int(os.path.exists(os.path.join(root, '전처리 이전 데이터', 'sEEPP_label_raw', f'Total_{gid}.csv')))

    df.to_csv(status_table, index=False)


if __name__ == '__main__':
    status_update()