from sklearn.model_selection import train_test_split

import gc

from utils import *


# ENC data, label 생성
class GenerateENCData(FilePathFinder):
    def __init__(self, genome_id):
        super().__init__(genome_id, 'ENC')
        self.read_path = self.read_path()

        self.data_dict = {}
        self.total_df = pd.DataFrame()
        self.data_cols = ['gene_id', 'amino_acid']
        self.label_cols = ['gene_id', 'enzyme', 'ec_class', 'ec_subclass', 'ec_subsubclass', 'ec_serial', 'score']

        all_ec_4d = pd.read_csv(class_ec_4digit)
        self.all_ec_4d = list(all_ec_4d['EC_number_4d'])

    # 원데이터 읽기, 오류 사전 제거, 전처리
    def load_data(self):
        print('  loading data... ')
        t = time.time()

        # Merge.XXXXX 파일 읽고 PK->gene_id 변환
        data_df = pd.read_csv(self.read_path['data'])
        data_df = data_df[['PK', 'amino_acid']]
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
                        data_df.loc[i, 'amino_acid'] = seq[:-1]
                    else:
                        drop_idx.append(i)
                else:
                    drop_idx.append(i)
        data_df = data_df.drop(drop_idx)
        print(f"{len(drop_idx)} rows removed")

        # 사전에 정의된 길이로 고정
        print('    sequence padding & cropping...')
        data_df['amino_acid'] = data_df['amino_acid'].map(
            lambda x: 'N'*maxlen['amino_acid'] if type(x) != str
            else (x[:maxlen['amino_acid']] if len(x) > maxlen['amino_acid'] else x+'N'*(maxlen['amino_acid']-len(x))))

        # BRENDA annotation result 파일 읽고 Protein name->gene_id 변환
        label_df = pd.read_csv(self.read_path['label'])
        label_df['Protein name'] = label_df['Protein name'].map(lambda x: x[2:-1])
        label_df.rename(columns={'Protein name': 'gene_id'}, inplace=True)

        # 클래스로 정의되지 않는 데이터 제거 (1~3자리 EC, 또는 B 포함하는 EC)\
        # score 기준치 이하 데이터 제거 (non enzyme 으로 취급)
        print('    checking label error... ', end='')
        drop_idx = []
        for i in range(len(label_df)):
            if label_df['EC Number'].iloc[i] not in self.all_ec_4d:
                drop_idx.append(i)
            elif label_df['Score'].iloc[i] < enc_score_threshold:
                drop_idx.append(i)
        label_df = label_df.drop(drop_idx)

        # gene_id 하나당 score 가장 높은 EC number 하나만 남김
        label_df = label_df.sort_values(by=['gene_id', 'Score'], ascending=[True, False])
        label_df = label_df.drop_duplicates(subset=['gene_id'], keep='first')

        # data, label 을 left join, EC number 가 없다면 non enzyme 취급 (EC.0.0.0.0)
        total = pd.merge(data_df, label_df, left_on='gene_id', right_on='gene_id', how='left')
        total['EC Number'] = total['EC Number'].replace(np.nan, 'EC 0.0.0.0')
        total['Score'] = total['Score'].replace(np.nan, 0)
        print(f"{len(drop_idx)} rows removed, {len(total)} left")

        self.total_df = total
        del data_df, label_df, total
        print('  time spent: %.2f(s)' % (time.time()-t))

    # train, validation, test set 분할
    def split_train_val_test(self, ratio):
        print('  splitting data... ')
        t = time.time()
        valtest_ratio = sum(ratio[1:]) / sum(ratio)
        test_ratio = ratio[2] / sum(ratio[1:])

        # element 개수가 하나뿐인 class 가 있다면 sklearn 제공 함수로는 불가능
        train, valtest = train_test_split_advanced(self.total_df, valtest_ratio)
        val, test = train_test_split_advanced(valtest, test_ratio)

        self.data_dict['train'] = train
        self.data_dict['validation'] = val
        self.data_dict['test'] = test
        del train, val, test

        print(f"    train: {len(self.data_dict['train'])}, validation: {len(self.data_dict['validation'])}, test: {len(self.data_dict['test'])}")
        print('  time spent: %.2f(s)' % (time.time()-t))

    # 데이터 저장
    def save_file(self, save_path, df):
        # data 저장
        print('    inputs: ', end='')
        save_result(df[self.data_cols], save_path['data_folder'], '', 'amino_acid', 'csv')
        print('amino_acid')

        # label 사전 정의된 포멧에 맞춰 저장
        print('    labels: ', end='')
        df['enzyme'] = df['EC Number'].map(lambda x: 0 if x == 'EC 0.0.0.0' else 1)
        df['ec_class'] = df['EC Number'].map(lambda x: int(x[3:].split('.')[0]))
        df['ec_subclass'] = df['EC Number'].map(lambda x: int(x[3:].split('.')[1]))
        df['ec_subsubclass'] = df['EC Number'].map(lambda x: int(x[3:].split('.')[2]))
        df['ec_serial'] = df['EC Number'].map(lambda x: int(x[3:].split('.')[3]))
        df['score'] = df['Score'].map(lambda x: enc_max_score if x > enc_max_score else x)
        # df.rename(columns={'Score': 'score'}, inplace=True)
        label = df[self.label_cols]
        label = label.sort_values(by=['gene_id'], ascending=True)

        dict = label.to_dict('records')
        enc_label = {'plant_name': self.plant_name,
                     'csv_file': '/data/ENC/' + self.plant_name + '/amino_acid.csv',
                     'label_count': len(dict), 'label': dict}
        save_result(enc_label, save_path['label_folder'], '', self.plant_name, 'json')
        print(self.plant_name)

        del label, dict, enc_label

    # 모든 데이터 저장
    def save_all_file(self):
        for set_name in ['train', 'validation', 'test']:
            print(f'  generating {set_name} set...')
            t = time.time()
            self.save_file(self.save_path(set_name), self.data_dict[set_name])
            print('  time spent: %.2f(s)' % (time.time() - t))

    # 전체 프로세스 실행
    def generate_and_save_data(self):
        t = time.time()
        print("<Generate ENC Data>")
        print(f"GenomeName: {self.genome_name}, GenomeID: {self.genome_id}")
        self.load_data()
        self.split_train_val_test([8, 1, 1])
        self.save_all_file()
        print('total time spent: %.2f(s)' % (time.time() - t))


# unique element 도 고려하여 분할
def train_test_split_advanced(data, test_size):
    ec_cnt = data['EC Number'].value_counts()
    data['count'] = data['EC Number'].map(lambda x: ec_cnt[x])

    # unique element 들은 random split 하여 나눔
    ec_single = data[data['count'] == 1].drop(['count'], axis=1)
    ec_single.sample()
    ec_single_train = ec_single.iloc[:int(len(ec_single) * (1 - test_size))]
    ec_single_test = ec_single.iloc[int(len(ec_single) * (1 - test_size)):]

    # 2게 이상의 element 들은 sklearn train_test_split 사용
    ec_multi = data[data['count'] > 1].drop(['count'], axis=1)
    ec_multi_train, ec_multi_test = train_test_split(
        ec_multi, test_size=test_size, stratify=ec_multi['EC Number']
    )

    # 두 방법으로 나눈 데이터 각각 병합
    train = pd.concat([ec_single_train, ec_multi_train], axis=0).sort_values(by=['gene_id'])
    test = pd.concat([ec_single_test, ec_multi_test], axis=0).sort_values(by=['gene_id'])
    return train, test


if __name__ == '__main__':
    set_seed()
    target_genome_id = 14225
    t = time.time()
    work = GenerateENCData(target_genome_id)
    work.generate_and_save_data()
    del work


