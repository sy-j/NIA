import random
import numpy as np
from skmultilearn.model_selection import iterative_train_test_split
from sklearn.model_selection import train_test_split
from sklearn.utils import check_random_state
import time

import gc

from utils import *

pd.set_option('mode.chained_assignment', None)  # settingwithcopywarning 무시


def set_seed(seed=911):
    random.seed(seed)
    np.random.seed(seed)
    check_random_state(seed)


class GenerateEEPPData(FilePathFinder):
    def __init__(self, genome_id):
        super().__init__(genome_id, 'sEEPP')
        self.read_path = self.read_path()

        self.data_dict = {}
        self.data_df = pd.DataFrame()
        self.label_df = pd.DataFrame()
        self.codon_cols = []
        self.label_cols = []

    def load_data(self):
        print('loading data... ', end='')
        t = time.time()
        data_df = pd.read_csv(self.read_path['data'])
        data_df['PK'] = data_df['PK'].map(lambda x: x[2:-1])
        data_df.rename(columns={'PK': 'gene_id'}, inplace=True)
        self.data_df = data_df[data_df.columns[:71]]
        self.codon_cols = data_df.columns[7:71]

        label_df = pd.read_csv(self.read_path['label'])
        label_df.dropna(axis=1, how='all', inplace=True)
        label_df['gene_id'] = label_df['gene_id'].map(lambda x: x[2:-1])
        self.label_df = label_df
        self.label_cols = list(label_df.columns)
        del data_df, label_df
        print('time spent: %.2f(s)' % (time.time()-t))

    def split_train_val_test(self, ratio):
        print('splitting data... ', end='')
        t = time.time()
        '''
        ratio(string) : minor group ratio, split by ','
        ex) '0.2,0.5' 8:2 로 나누고 이후 2 를 1:1로 나눔
        self.data_dict 에 키 값으로 저장
        '''
        valtest_ratio = sum(ratio[1:]) / sum(ratio)
        test_ratio = ratio[2] / sum(ratio[1:])
        if len(self.label_cols) == 2:
            # 라벨 컬럼 한개인 경우 sklearn.model_selection.train_test_split 이용
            train_x, valtest_x, train_y, valtest_y = train_test_split(
                self.data_df, self.label_df,
                test_size=valtest_ratio,
                stratify=self.label_df[self.label_cols[-1]]
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
        else:
            # multilabel인 경우 skmultilearn.model_selection 이용
            label_df = self.label_df.drop(columns=['gene_id'])
            temp_label_cols = label_df.columns
            temp_data_cols = self.data_df.columns
            out_arr = np.array(label_df)
            in_arr = np.array(self.data_df)
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
        print('time spent: %.2f(s)' % (time.time()-t))

    def save_x(self, dir, df, codon_cols):
        print('inputs: ', end='')
        for col in ['promoter', 'terminator', 'utr5', 'utr3', 'cds']:
            save_result(df[['gene_id', col]], dir, self.genome_name, col, 'csv')
            # print(f"process done: {self.genome_name}: {col}.csv ")
            print(col, end=' ')

        gene_id_list = list(df['gene_id'])
        temp_arr = np.array(df[codon_cols])
        temp_arr = temp_arr.flatten()
        ratio = []
        gene_id = []
        codon = list(range(1, 65)) * len(df)
        for i in range(len(df)):
            ratio.extend(temp_arr[64 * i:64 * (i + 1)] / sum(temp_arr[64 * i:64 * (i + 1)]))
            gene_id += [gene_id_list[i]] * 64

        codon_usage = pd.DataFrame({
            'gene_id': gene_id,
            'codon': codon,
            'count': temp_arr,
            'ratio': ratio
        })
        save_result(codon_usage, dir, self.genome_name, 'codon_usage', 'csv')
        print('codon_usage')

    def save_y(self, dir, df, label_cols):
        print('labels: ', end='')
        save_cols = label_cols[1:]
        for col in save_cols:
            temp_df = df[['gene_id', col]]
            temp_df.loc[:, 'organ'] = [col] * len(df)
            temp_df.columns = ['gene_id', 'expression_valid', 'organ']
            temp_df = temp_df[['gene_id', 'organ', 'expression_valid']]
            temp_dict = temp_df.to_dict(orient='records')
            seepp_label = {}
            seepp_label['plant_name'] = self.plant_name
            seepp_label['csv_file'] = [
                f"data/sEEPP/{self.plant_name}/cds.csv",
                f"data/sEEPP/{self.plant_name}/codon_usage.csv",
                f"data/sEEPP/{self.plant_name}/promoter.csv",
                f"data/sEEPP/{self.plant_name}/terminator.csv",
                f"data/sEEPP/{self.plant_name}/utr3.csv",
                f"data/sEEPP/{self.plant_name}/utr5.csv"
            ]
            seepp_label['organ_list'] = label_cols
            seepp_label['label_count'] = len(temp_df)
            seepp_label['label'] = temp_dict
            fname = f"{self.plant_name}_{col}"
            save_result(seepp_label, dir, self.genome_name, fname, 'json')
            print(col, end=' ')
        print('')

    def save_file(self):
        print('generating train set...')
        save_path = self.save_path('train')
        t = time.time()
        self.save_x(save_path['data_folder'], self.data_dict['train_x'], self.codon_cols)
        self.save_y(save_path['label_folder'], self.data_dict['train_y'], self.label_cols)
        print('time spent: %.2f(s)' % (time.time() - t))

        print('generating validation set...')
        save_path = self.save_path('validation')
        t = time.time()
        self.save_x(save_path['data_folder'], self.data_dict['validation_x'], self.codon_cols)
        self.save_y(save_path['label_folder'], self.data_dict['validation_y'], self.label_cols)
        print('time spent: %.2f(s)' % (time.time() - t))

        print('generating test set...')
        save_path = self.save_path('test')
        t = time.time()
        self.save_x(save_path['data_folder'], self.data_dict['test_x'], self.codon_cols)
        self.save_y(save_path['label_folder'], self.data_dict['test_y'], self.label_cols)
        print('time spent: %.2f(s)' % (time.time() - t))

    def generate_and_save_data(self):
        t = time.time()
        print(f"GenomeName: {self.genome_name}, GenomeID: {self.genome_id}")
        self.load_data()
        self.split_train_val_test([8, 1, 1])
        self.save_file()
        print('total time spent: %.2f(s)' % (time.time() - t))


if __name__ == '__main__':
    set_seed()
    target_genome_id = 14205
    t = time.time()
    work = GenerateEEPPData(target_genome_id)
    work.generate_and_save_data()
    # work.load_data()
    # work.split_train_val_test([8, 1, 1])
    # work.save_file()
    del work
