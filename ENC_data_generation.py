import pandas as pd
from skmultilearn.model_selection import iterative_train_test_split
from sklearn.model_selection import train_test_split

import gc

from utils import *


class GenerateENCData(FilePathFinder):
    def __init__(self, genome_id):
        super().__init__(genome_id, 'ENC')
        self.read_path = self.read_path()

        self.data_dict = {}
        self.total_df = pd.DataFrame()
        self.data_cols = ['gene_id', 'amino_acid']
        self.label_cols = ['gene_id', 'enzyme', 'ec_class', 'ec_subclass', 'ec_subsubclass', 'ec_serial', 'score']

        all_ec_4d = pd.read_csv(class_ec_4digit)
        all_ec_3d = pd.read_csv(class_ec_3digit)
        self.all_ec_4d = list(all_ec_4d['EC_number_4d'])
        self.all_ec_3d = list(all_ec_3d['EC_number_3d'])

    def load_data(self):
        print('  loading data... ')
        t = time.time()

        # Merge.XXXXX 파일 읽고 PK->gene_id 변환
        data_df = pd.read_csv(self.read_path['data'])
        data_df = data_df[['PK', ㅜ'amino_acid']]
        data_df['PK'] = data_df['PK'].map(lambda x: x[2:-1])
        data_df.rename(columns={'PK': 'gene_id'}, inplace=True)

        # 아미노산 서열 길이 10 이하, 21개 문자 이외 포함하는 서열 제거
        # 단, 마지막에 ] 들어있는 것은 추가 파싱후 사용, 추후 버그 수정시 해당 기능도 삭제
        print('    checking data error... ', end='')
        drop_idx = []
        for i in range(len(data_df)):
            seq = data_df['amino_acid'].iloc[i]
            if len(str(seq)) < 10:
                drop_idx.append(i)
            elif len(set(seq) - set('ACDEFGHIKLMNPQRSTVWXY')) > 0:
                if seq.count(']') >= 1:
                    if seq.count(']') == 1 and seq[-1] == ']':
                        # print(data_df['amino_acid'].iloc[i])
                        data_df.loc[i, 'amino_acid'] = seq[:-1]
                        # print(data_df['amino_acid'].iloc[i])
                    else:
                        # print('found you')
                        drop_idx.append(i)
                else:
                    drop_idx.append(i)
        data_df = data_df.drop(drop_idx)
        print(f"{len(drop_idx)} rows removed")
        print('    sequence padding & cropping...')
        data_df['amino_acid'] = data_df['amino_acid'].map(
            lambda x: 'N'*maxlen['amino_acid'] if type(x) != str
            else (x[:maxlen['amino_acid']] if len(x) > maxlen['amino_acid'] else x+'N'*(maxlen['amino_acid']-len(x))))

        data_cols = data_df.columns

        label_df = pd.read_csv(self.read_path['label'])
        label_df['Protein name'] = label_df['Protein name'].map(lambda x: x[2:-1])
        label_df.rename(columns={'Protein name': 'gene_id'}, inplace=True)

        print('    checking label error... ', end='')
        drop_idx = []
        for i in range(len(label_df)):
            if label_df['EC Number'].iloc[i] not in self.all_ec_4d:
                drop_idx.append(i)
        label_df = label_df.drop(drop_idx)
        label_cols = label_df.columns

        label_df = label_df.sort_values(by=['gene_id', 'Score'], ascending=[True, False])
        label_df = label_df.drop_duplicates(subset=['gene_id'], keep='first')
        total = pd.merge(data_df, label_df, left_on='gene_id', right_on='gene_id', how='left')
        total['EC Number'] = total['EC Number'].replace(np.nan, 'EC 0.0.0.0')
        total['Score'] = total['Score'].replace(np.nan, 0)
        print(f"{len(drop_idx)} rows removed, {len(total)} left")

        self.total_df = total
        del data_df, label_df
        print('  time spent: %.2f(s)' % (time.time()-t))

    def split_train_val_test(self, ratio):
        print('  splitting data... ')
        t = time.time()
        valtest_ratio = sum(ratio[1:]) / sum(ratio)
        test_ratio = ratio[2] / sum(ratio[1:])

        train, valtest = train_test_split_advanced(self.total_df, valtest_ratio)
        val, test = train_test_split_advanced(valtest, test_ratio)

        self.data_dict['train'] = train
        self.data_dict['validation'] = val
        self.data_dict['test'] = test

        print(f"    train: {len(self.data_dict['train'])}, validation: {len(self.data_dict['validation'])}, test: {len(self.data_dict['test'])}")
        print('  time spent: %.2f(s)' % (time.time()-t))

    def save_x(self, save_path, df):
        print('    inputs: ', end='')
        save_result(df[['gene_id', 'amino_acid']], save_path['data_folder'], '', 'amino_acid', 'csv')
        print('amino_acid')

        print('    labels: ', end='')
        df['enzyme'] = df['EC Number'].map(lambda x: 0 if x == 'EC 0.0.0.0' else 1)
        df['ec_class'] = df['EC Number'].map(lambda x: int(x[3:].split('.')[0]))
        df['ec_subclass'] = df['EC Number'].map(lambda x: int(x[3:].split('.')[1]))
        df['ec_subsubclass'] = df['EC Number'].map(lambda x: int(x[3:].split('.')[2]))
        df['ec_serial'] = df['EC Number'].map(lambda x: int(x[3:].split('.')[3]))
        df.rename(columns={'Score': 'score'}, inplace=True)
        label = df[self.label_cols]
        label = label.sort_values(by=['gene_id'], ascending=True)

        dict = label.to_dict('records')
        enc_label = {}
        enc_label['plant_name'] = self.plant_name
        enc_label['csv_file'] = '/data/ENC/' + self.plant_name + '/amino_acid.csv'
        enc_label['label_count'] = len(dict)
        enc_label['label'] = dict
        save_result(enc_label, save_path['label_folder'], '', self.plant_name, 'json')
        print(self.plant_name)

    def save_file(self):
        print('  generating train set...')
        save_path = self.save_path('train')
        t = time.time()
        self.save_x(save_path, self.data_dict['train'])
        print('  time spent: %.2f(s)' % (time.time() - t))

        print('  generating validation set...')
        save_path = self.save_path('validation')
        t = time.time()
        self.save_x(save_path, self.data_dict['validation'])
        print('  time spent: %.2f(s)' % (time.time() - t))

        print('  generating test set...')
        save_path = self.save_path('test')
        t = time.time()
        self.save_x(save_path, self.data_dict['test'])
        print('  time spent: %.2f(s)' % (time.time() - t))

    def generate_and_save_data(self):
        t = time.time()
        print("<Generate ENC Data>")
        print(f"GenomeName: {self.genome_name}, GenomeID: {self.genome_id}")
        self.load_data()
        self.split_train_val_test([8, 1, 1])
        self.save_file()
        print('total time spent: %.2f(s)' % (time.time() - t))


def train_test_split_advanced(data, test_size):
    ec_cnt = data['EC Number'].value_counts()
    data['count'] = data['EC Number'].map(lambda x: ec_cnt[x])

    ec_single = data[data['count'] == 1].drop(['count'], axis=1)
    ec_single.sample()
    ec_single_train = ec_single.iloc[:int(len(ec_single) * (1 - test_size))]
    ec_single_test = ec_single.iloc[int(len(ec_single) * (1 - test_size)):]

    ec_multi = data[data['count'] > 1].drop(['count'], axis=1)

    ec_multi_train, ec_multi_test = train_test_split(
        ec_multi, test_size=test_size, stratify=ec_multi['EC Number']
    )

    train = pd.concat([ec_single_train, ec_multi_train], axis=0).sort_values(by=['gene_id'])
    test = pd.concat([ec_single_test, ec_multi_test], axis=0).sort_values(by=['gene_id'])
    return train, test


if __name__ == '__main__':
    set_seed()
    target_genome_id = 14205
    t = time.time()
    work = GenerateENCData(target_genome_id)
    work.generate_and_save_data()
    del work


