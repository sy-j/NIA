import os
import pandas as pd
import json
import tqdm
import time
import random
import numpy as np

root = r'D:\NIA\학습용 데이터'
status_table = os.path.join(root, 'Genome_data_total_status.csv')
status_table2 = os.path.join(root, 'Genome_data_size.csv')
status_table3 = os.path.join(root, 'final_status_check.csv')
class_ec_4digit = os.path.join(root, 'class_ec_4digit.csv')
class_ec_3digit = os.path.join(root, 'class_ec_3digit.csv')
enc_score_threshold = 100
enc_max_score = 400
maxlen = {
    'promoter': 1000,
    'terminator': 1000,
    'utr5': 2000,
    'utr3': 2000,
    'cds': 5000,
    'amino_acid': 2000
}
amino_acid_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H',
                   'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                   'R', 'S', 'T', 'V', 'W', 'Y', 'X']

nucleotide_seq_cols = ['promoter', 'terminator', 'utr5', 'utr3', 'cds']
sEEPP_input_files = ['cds', 'promoter', 'terminator', 'utr5', 'utr3', 'codon_usage']
organs = ['leaf', 'root', 'stem', 'bud', 'flower']
raw_data_path = {
    'genome': os.path.join(root, '전처리 이전 데이터', 'gene_input_raw'),
    'RNASeq': os.path.join(root, '전처리 이전 데이터', 'sEEPP_label_raw'),
    'BRENDA': os.path.join(root, '전처리 이전 데이터', 'ENC_label_raw')
}


# genome_id, 모델명에 대한 모든 파일경로 관리
class FilePathFinder:
    def __init__(self, target_genome_id, model):
        self.genome_id = target_genome_id
        self.model = model
        plant_name, genome_name = get_name(target_genome_id)
        self.plant_name = plant_name
        self.genome_name = genome_name

    def read_path(self):
        if self.model == 'sEEPP':
            dict = {
                'data': os.path.join(root, '전처리 이전 데이터', 'gene_input_raw', f'Merge.{self.genome_id}.csv'),
                'label': os.path.join(root, '전처리 이전 데이터', 'sEEPP_label_raw', f'Total_{self.genome_id}.csv')
            }
        if self.model == 'ENC':
            dict = {
                'data': os.path.join(root, '전처리 이전 데이터', 'gene_input_raw', f'Merge.{self.genome_id}.csv'),
                'label': os.path.join(root, '전처리 이전 데이터', 'ENC_label_raw', f'{self.genome_name}.csv')
            }
        return dict

    def save_path(self, set, purpose='data_all'):
        if self.model == 'sEEPP':
            dict = {
                'data_folder': os.path.join(root, purpose, set, 'data', 'sEEPP', self.plant_name),
                'label_folder': os.path.join(root, purpose, set, 'label', 'sEEPP', self.plant_name),
                'cds': os.path.join(root, purpose, set, 'data', 'sEEPP', self.plant_name, 'cds.csv'),
                'promoter': os.path.join(root, purpose, set, 'data', 'sEEPP', self.plant_name, 'promoter.csv'),
                'terminator': os.path.join(root, purpose, set, 'data', 'sEEPP', self.plant_name, 'terminator.csv'),
                'utr5': os.path.join(root, purpose, set, 'data', 'sEEPP', self.plant_name, 'utr5.csv'),
                'utr3': os.path.join(root, purpose, set, 'data', 'sEEPP', self.plant_name, 'utr3.csv'),
                'codon_usage': os.path.join(root, purpose, set, 'data', 'sEEPP', self.plant_name, 'codon_usage.csv'),
                'leaf': os.path.join(root, purpose, set, 'label', 'sEEPP', self.plant_name, f'{self.plant_name}_leaf.json'),
                'root': os.path.join(root, purpose, set, 'label', 'sEEPP', self.plant_name, f'{self.plant_name}_root.json'),
                'stem': os.path.join(root, purpose, set, 'label', 'sEEPP', self.plant_name, f'{self.plant_name}_stem.json'),
                'bud': os.path.join(root, purpose, set, 'label', 'sEEPP', self.plant_name, f'{self.plant_name}_bud.json'),
                'flower': os.path.join(root, purpose, set, 'label', 'sEEPP', self.plant_name, f'{self.plant_name}_flower.json')
            }
        if self.model == 'ENC':
            dict = {
                'data_folder': os.path.join(root, purpose, set, 'data', 'ENC', self.plant_name),
                'label_folder': os.path.join(root, purpose, set, 'label', 'ENC', self.plant_name),
                'amino_acid': os.path.join(root, purpose, set, 'data', 'ENC', self.plant_name, 'amino_acid.csv'),
                'label': os.path.join(root, purpose, set, 'label', 'ENC', self.plant_name, f'{self.plant_name}.json')
            }
        return dict


# random seed 설정
def set_seed(seed=911):
    random.seed(seed)
    np.random.seed(seed)


# 경로 생성
def create_folder(dir):
    try:
        if not os.path.exists(dir):
            os.makedirs(dir)
    except OSError:
        print('Error: Creating directory. ' + dir)


# 경로 존재여부 확인, 덮어쓰기 옵션 선택
def file_check(settings, name):
    dir = os.path.join(root, settings['type'], settings['model'], name)
    if os.path.exists(dir):
        print('%s file already exists. proceed anyway? [y/n]' % name)
        if input() == 'y':
            return 1
        else:
            return 0
    else:
        return 1


# genome id에 해당하는 종명, 유전체명 반환
def get_name(target_genome_id):
    st = pd.read_csv(status_table)
    tmp = st[st['genome_id'] == target_genome_id]
    plant_name = tmp.plant_name.iloc[0]
    genome_name = tmp.genome_name.iloc[0]
    return plant_name, genome_name


# def type_convert(o):
#     if isinstance(o, np.int64):
#         return int(o)
#     raise TypeError


# 파일 저장
def save_result(data, root_name, folder_name, file_name, file_type):
    dir = os.path.join(root_name, folder_name)
    create_folder(dir)
    dir = os.path.join(dir, file_name+'.'+file_type)
    if file_type == 'csv':
        data.to_csv(dir, index=False)
    if file_type == 'json':
        with open(dir, 'w') as f:
            json.dump(data, f, ensure_ascii=False, indent=4)


# BRENDA annotation 결과파일 이름 변경
def rename_annotation_result():
    print('renaming BRENDA annotation result files..')
    raw_label_dir = os.path.join(root, '전처리 이전 데이터', 'ENC_label_raw')
    dir_list = os.listdir(raw_label_dir)
    for dir in dir_list:
        # print(os.path.basename(dir))
        if ' ' in dir:
            os.renames(os.path.join(raw_label_dir, dir), os.path.join(raw_label_dir, os.path.basename(dir).replace(' ', '_')))


def data_count(model, name):
    input_dirs = os.listdir(os.path.join(root, 'data', model, name))
    label_dirs = os.listdir(os.path.join(root, 'label', model, name))
    input_cnt, label_cnt = 0, 0
    for dir in input_dirs:
        df = pd.read_csv(os.path.join(os.path.join(root, 'data', model, name), dir))
        input_cnt += len(df)
    for dir in label_dirs:
        # print(dir)
        with open(os.path.join(root, 'label', model, name, dir), 'r') as tmp:
            tmp = json.load(tmp)
            label_cnt += len(tmp['label'])
    return input_cnt, label_cnt


# 전체 EC number 확인 (1~3자리까지만 나온 것, B 포함된 것 제외)
def define_total_ec():
    file_list = os.listdir(raw_data_path['BRENDA'])[:-1]
    ec = []
    for file in tqdm.tqdm(file_list):
        df = pd.read_csv(os.path.join(raw_data_path['BRENDA'], file))
        df = df[df['Score'] >= enc_score_threshold]
        ec += list(df['EC Number'].unique())

    ec = list(set(ec))
    unique_ec_4d = []
    for item in tqdm.tqdm(ec):
        if item.count('.') == 3 and item.count('B') == 0 and item not in unique_ec_4d:
            unique_ec_4d.append(item)
    unique_ec_4d.sort()
    ec_4d = pd.DataFrame(unique_ec_4d)
    ec_4d.columns = ['EC_number_4d']
    ec_4d.to_csv(class_ec_4digit, index=False)

    unique_ec_3d = []
    for item in tqdm.tqdm(unique_ec_4d):
        tmp = '.'.join(item.split('.')[:-1])
        if tmp not in unique_ec_3d:
            unique_ec_3d.append(tmp)
    unique_ec_3d.sort()
    ec_3d = pd.DataFrame(unique_ec_3d)
    ec_3d.columns = ['EC_number_3d']
    ec_3d.to_csv(class_ec_3digit, index=False)


if __name__ == '__main__':
    define_total_ec()