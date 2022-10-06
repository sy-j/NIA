from skmultilearn.model_selection import iterative_train_test_split
from sklearn.model_selection import train_test_split

import gc

from utils import *

# setting with copy warning 무시
pd.set_option('mode.chained_assignment', None)


# sEEPP data, label 생성
class GenerateEEPPData(FilePathFinder):
    def __init__(self, genome_id):
        super().__init__(genome_id, 'sEEPP')
        self.read_path = self.read_path()

        self.data_dict = {}
        self.total_df = pd.DataFrame()
        self.codon_cols = []
        self.data_cols = []
        self.label_cols = []

    # 원데이터 읽기, 오류 사전 제거, 전처리
    def load_data(self):
        print('  loading data... ')
        t = time.time()

        # Merge.XXXXX 파일 읽고 PK->gene_id 변환
        data_df = pd.read_csv(self.read_path['data'])
        data_df['PK'] = data_df['PK'].map(lambda x: x[2:-1])
        data_df.rename(columns={'PK': 'gene_id'}, inplace=True)

        # ATGCN 이외 문자를 갖는 서열 제거
        print('    checking error... ', end='')
        drop_idx = []
        for i in range(len(data_df)):
            seq = data_df.iloc[i]
            for col in nucleotide_seq_cols:
                if str(seq[col]) == 'nan':
                    continue
                elif len(set(seq[col]) - set('ATGCN')) > 0:
                    drop_idx.append(i)
                    break
        data_df = data_df.drop(drop_idx)
        print(f"{len(drop_idx)} rows removed, {len(data_df)} left")

        # 사전에 정의된 길이로 고정
        print('    sequence padding & cropping...')
        for col in ['cds', 'terminator', 'utr3']:
            data_df[col] = data_df[col].map(
                lambda x: 'N'*maxlen[col] if type(x) != str
                else (x[:maxlen[col]] if len(x) > maxlen[col] else x+'N'*(maxlen[col]-len(x))))
        for col in ['promoter', 'utr5']:
            data_df[col] = data_df[col].map(
                lambda x: 'N'*maxlen[col] if type(x) != str
                else (x[-maxlen[col]:] if len(x) > maxlen[col] else 'N'*(maxlen[col]-len(x))+x))

        self.data_cols = data_df.columns

        # Total_XXXXX 파일 읽고 gene_id 변환
        label_df = pd.read_csv(self.read_path['label'])
        label_df.dropna(axis=1, how='all', inplace=True)
        label_df['gene_id'] = label_df['gene_id'].map(lambda x: x[2:-1])

        self.label_cols = label_df.columns

        total = pd.merge(data_df, label_df, left_on='gene_id', right_on='gene_id', how='inner')

        self.total_df = total
        self.codon_cols = data_df.columns[7:71]
        del data_df, label_df, total
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


if __name__ == '__main__':
    set_seed()
    target_genome_id = 14205
    t = time.time()
    work = GenerateEEPPData(target_genome_id)
    work.generate_and_save_data()
    del work
