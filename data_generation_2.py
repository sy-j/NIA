from skmultilearn.model_selection import iterative_train_test_split
from sklearn.model_selection import train_test_split

import gc

from utils import *


# setting with copy warning 무시
pd.set_option('mode.chained_assignment', None)


# sEEPP data, label 생성
class GenerateEEPPData2(FilePathFinder):
    def __init__(self, genome_id):
        super().__init__(genome_id, 'sEEPP')
        self.data_dict = {}
        self.total_df = pd.DataFrame()
        self.codon_cols = []
        self.data_cols = []
        self.label_cols = []

    # 원데이터 읽기, 오류 사전 제거, 전처리
    def load_data(self):
        print('  loading data... ')
        t = time.time()

        read_path = os.path.join(root, '개수 조정된 원데이터', 'sEEPP_ENC', f'{self.genome_id}_crop.csv')
        data_df = pd.read_csv(read_path)
        data_df.dropna(axis=1, how='all', inplace=True)
        self.total_df = data_df
        self.data_cols = data_df.columns[:71]
        self.codon_cols = data_df.columns[7:71]

        label_df = pd.read_csv(self.read_path()['label'])
        label_df.dropna(axis=1, how='all', inplace=True)
        self.label_cols = label_df.columns
        del data_df
        print('  time spent: %.2f(s)' % (time.time()-t))

    # train, validation, test set 분할
    def split_train_val_test(self, ratio):
        print('  splitting data... ')
        t = time.time()
        valtest_ratio = sum(ratio[1:]) / sum(ratio)
        test_ratio = ratio[2] / sum(ratio[1:])

        data_df = self.total_df[self.data_cols]  # 추후 total_df로 코드 정리 필요
        label_df = self.total_df[self.label_cols]

        if len(self.label_cols) == 2:
            # 라벨 컬럼 한개인 경우 sklearn.model_selection.train_test_split 이용
            train_x, valtest_x, train_y, valtest_y = train_test_split(
                data_df, label_df,
                test_size=valtest_ratio,
                stratify=label_df[self.label_cols[-1]]
            )
            self.data_dict['train_x'] = train_x
            self.data_dict['train_y'] = train_y
            val_x, test_x, val_y, test_y = train_test_split(
                valtest_x, valtest_y,
                test_size=test_ratio,
                stratify=valtest_y[self.label_cols[-1]]
            )
            self.data_dict['validation_x'] = val_x
            self.data_dict['validation_y'] = val_y
            self.data_dict['test_x'] = test_x
            self.data_dict['test_y'] = test_y
            del train_x, valtest_x, train_y, valtest_y, val_x, test_x, val_y, test_y, data_df, label_df
        else:
            # multilabel인 경우 skmultilearn.model_selection 이용
            label_df = label_df.drop(columns=['gene_id'])
            temp_label_cols = label_df.columns
            temp_data_cols = data_df.columns
            out_arr = np.array(label_df)
            in_arr = np.array(data_df)
            train_x, train_y, valtest_x, valtest_y = iterative_train_test_split(in_arr, out_arr,
                                                                                test_size=valtest_ratio)
            self.data_dict['train_x'] = pd.DataFrame(train_x, columns=temp_data_cols)
            temp_y_id = pd.DataFrame({'gene_id': train_x[:, 0]})
            temp_y_label = pd.DataFrame(train_y, columns=temp_label_cols)
            self.data_dict['train_y'] = pd.concat([temp_y_id, temp_y_label], axis=1)
            val_x, val_y, test_x, test_y = iterative_train_test_split(valtest_x, valtest_y, test_size=test_ratio)
            self.data_dict['validation_x'] = pd.DataFrame(val_x, columns=temp_data_cols)
            temp_y_id = pd.DataFrame({'gene_id': val_x[:, 0]})
            temp_y_label = pd.DataFrame(val_y, columns=temp_label_cols)
            self.data_dict['validation_y'] = pd.concat([temp_y_id, temp_y_label], axis=1)
            self.data_dict['test_x'] = pd.DataFrame(test_x, columns=temp_data_cols)
            temp_y_id = pd.DataFrame({'gene_id': test_x[:, 0]})
            temp_y_label = pd.DataFrame(test_y, columns=temp_label_cols)
            self.data_dict['test_y'] = pd.concat([temp_y_id, temp_y_label], axis=1)
            del in_arr, out_arr, train_x, train_y, valtest_x, valtest_y, temp_y_id, temp_y_label, data_df, label_df
        print(f"    train: {len(self.data_dict['train_x'])}, validation: {len(self.data_dict['validation_x'])}, test: {len(self.data_dict['test_x'])}")
        print('  time spent: %.2f(s)' % (time.time()-t))

    # data 저장
    def save_x(self, dir, df, codon_cols):
        print('    inputs: ', end='')

        # 염기서열들 저장
        for col in nucleotide_seq_cols:
            save_result(df[['gene_id', col]], dir, '', col, 'csv')
            print(col, end=' ')

        # codon usage 저장
        gene_id_list = list(df['gene_id'])
        temp_arr = np.array(df[codon_cols])
        temp_arr = temp_arr.flatten()
        ratio = []
        gene_id = []
        codon = list(range(1, 65)) * len(df)
        for i in range(len(df)):
            # print(temp_arr[64 * i:64 * (i + 1)])
            # print(len(temp_arr[64 * i:64 * (i + 1)]))
            # print(temp_arr[0])
            # print(len(temp_arr))
            ratio.extend(temp_arr[64 * i:64 * (i + 1)] / sum(temp_arr[64 * i:64 * (i + 1)]))
            gene_id += [gene_id_list[i]] * 64

        codon_usage = pd.DataFrame({
            'gene_id': gene_id,
            'codon': codon,
            'count': temp_arr,
            'ratio': ratio
        })
        save_result(codon_usage, dir, '', 'codon_usage', 'csv')
        print('codon_usage')

    # label 저장
    def save_y(self, dir, df, label_cols):
        print('    labels: ', end='')
        save_cols = list(label_cols[1:])
        for col in save_cols:
            temp_df = df[['gene_id', col]]
            temp_df.loc[:, 'organ'] = [col] * len(df)
            temp_df.columns = ['gene_id', 'expression_valid', 'organ']
            temp_df = temp_df[['gene_id', 'organ', 'expression_valid']]
            temp_dict = temp_df.to_dict(orient='records')
            seepp_label = {'plant_name': self.plant_name, 'csv_file': [
                f"data/sEEPP/{self.plant_name}/cds.csv",
                f"data/sEEPP/{self.plant_name}/codon_usage.csv",
                f"data/sEEPP/{self.plant_name}/promoter.csv",
                f"data/sEEPP/{self.plant_name}/terminator.csv",
                f"data/sEEPP/{self.plant_name}/utr3.csv",
                f"data/sEEPP/{self.plant_name}/utr5.csv"
            ], 'organ_list': save_cols, 'label_count': len(temp_df), 'label': temp_dict}
            fname = f"{self.plant_name}_{col}"
            save_result(seepp_label, dir, '', fname, 'json')
            print(col, end=' ')
        print('')

    # 모든 데이터 저장
    def save_file(self):
        print('  generating train set...')
        save_path = self.save_path('train')
        t = time.time()
        self.save_x(save_path['data_folder'], self.data_dict['train_x'], self.codon_cols)
        self.save_y(save_path['label_folder'], self.data_dict['train_y'], self.label_cols)
        print('  time spent: %.2f(s)' % (time.time() - t))

        print('  generating validation set...')
        save_path = self.save_path('validation')
        t = time.time()
        self.save_x(save_path['data_folder'], self.data_dict['validation_x'], self.codon_cols)
        self.save_y(save_path['label_folder'], self.data_dict['validation_y'], self.label_cols)
        print('  time spent: %.2f(s)' % (time.time() - t))

        print('  generating test set...')
        save_path = self.save_path('test')
        t = time.time()
        self.save_x(save_path['data_folder'], self.data_dict['test_x'], self.codon_cols)
        self.save_y(save_path['label_folder'], self.data_dict['test_y'], self.label_cols)
        print('  time spent: %.2f(s)' % (time.time() - t))

    # 전체 프로세스 실행
    def generate_and_save_data(self):
        t = time.time()
        print("<Generate sEEPP Data>")
        print(f"GenomeName: {self.genome_name}, GenomeID: {self.genome_id}")
        self.load_data()
        self.split_train_val_test([8, 1, 1])
        self.save_file()
        print('total time spent: %.2f(s)' % (time.time() - t))


# ENC data, label 생성
class GenerateENCData2(FilePathFinder):
    def __init__(self, genome_id):
        super().__init__(genome_id, 'ENC')
        self.data_dict = {}
        self.total_df = pd.DataFrame()
        self.data_cols = ['gene_id', 'amino_acid']
        self.label_cols = ['gene_id', 'enzyme', 'ec_class', 'ec_subclass', 'ec_subsubclass', 'ec_serial', 'score']

    # 원데이터 읽기, 오류 사전 제거, 전처리
    def load_data(self):
        print('  loading data... ')
        t = time.time()

        if os.path.exists(os.path.join(root, '개수 조정된 원데이터', 'ENC', f'{self.genome_id}_crop.csv')):
            read_path = os.path.join(root, '개수 조정된 원데이터', 'ENC', f'{self.genome_id}_crop.csv')
        elif os.path.exists(os.path.join(root, '개수 조정된 원데이터', 'sEEPP_ENC', f'{self.genome_id}_crop.csv')):
            read_path = os.path.join(root, '개수 조정된 원데이터', 'sEEPP_ENC', f'{self.genome_id}_crop.csv')
        data_df = pd.read_csv(read_path)

        self.total_df = data_df
        del data_df
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
    ec_single.sample(frac=1)
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
    target_genome_id = 14560  #14573, 14569, 14580, 14560
    t = time.time()
    work = GenerateENCData2(target_genome_id)
    work.generate_and_save_data()
    del work

    work = GenerateEEPPData2(target_genome_id)
    work.generate_and_save_data()
    del work


