import os
import pandas as pd
import json
import tqdm
import time

# root = 'C:\\Users\\infoboss\\Desktop\\NIA_datafactory\\Gene'
# root = 'I:\\2022 프로젝트\\2022_NIA_AI학습용데이터\\데이터'
root = r'C:\Users\infoboss\Desktop\NIA\학습용 데이터'
# status_table = 'C:\\Users\\infoboss\\Desktop\\NIA_datafactory\\Gene\\' + 'Genome_data_total_status.csv'
status_table = os.path.join(root, 'Genome_data_total_status.csv')

maxlen = {
    'promoter': 1000,
    'terminator': 1000,
    'utr5': 2000,
    'utr3': 2000,
    'cds': 5000,
}

nucleotide_seq_cols = ['promoter', 'terminator', 'utr5', 'utr3', 'cds']
sEEPP_input_files = ['cds', 'promoter', 'terminator', 'utr5', 'utr3', 'codon_usage']
organs = ['leaf', 'root', 'stem', 'bud', 'flower']


def find_file_path(gid, set='raw'):  # 이하 모든 genome_name 은 최종적으로 plant_name 으로 통일되어야함
    plant_name, genome_name = get_name(gid)

    if set == 'raw':
        dict = {
            'genome': os.path.join(root, '전처리 이전 데이터', 'gene_input_raw', f'Merge.{gid}.csv'),
            'RNASeq': os.path.join(root, '전처리 이전 데이터', 'sEEPP_label_raw', f'Total_{gid}.csv'),
            'BRENDA': os.path.join(root, '전처리 이전 데이터', 'ENC_label_raw', f'{genome_name}.csv')
        }
    else:
        dict = {
            'cds': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', genome_name, 'cds.csv'),
            'promoter': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', genome_name, 'promoter.csv'),
            'terminator': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', genome_name, 'terminator.csv'),
            'utr5': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', genome_name, 'utr5.csv'),
            'utr3': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', genome_name, 'utr3.csv'),
            'codon_usage': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', genome_name, 'codon_usage.csv'),
            'amino_acid': os.path.join(root, '미분할 데이터', set, 'data', 'ENC', genome_name, 'amino_acid.csv'),
            'sEEPP_label': os.path.join(root, '미분할 데이터', set, 'label', 'sEEPP', genome_name),
            'ENC_label': os.path.join(root, '미분할 데이터', set, 'label', 'ENC', genome_name, f'{plant_name}.json')
        }
    return dict


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

    def save_path(self, set):
        if self.model == 'sEEPP':
            dict = {
                'data_folder': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', self.genome_name),
                'label_folder': os.path.join(root, '미분할 데이터', set, 'label', 'sEEPP', self.genome_name),
                'cds': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', self.genome_name, 'cds.csv'),
                'promoter': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', self.genome_name, 'promoter.csv'),
                'terminator': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', self.genome_name, 'terminator.csv'),
                'utr5': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', self.genome_name, 'utr5.csv'),
                'utr3': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', self.genome_name, 'utr3.csv'),
                'codon_usage': os.path.join(root, '미분할 데이터', set, 'data', 'sEEPP', self.genome_name, 'codon_usage.csv'),
                'leaf': os.path.join(root, '미분할 데이터', set, 'label', 'sEEPP', self.genome_name, f'{self.plant_name}_leaf.json'),
                'root': os.path.join(root, '미분할 데이터', set, 'label', 'sEEPP', self.genome_name, f'{self.plant_name}_root.json'),
                'stem': os.path.join(root, '미분할 데이터', set, 'label', 'sEEPP', self.genome_name, f'{self.plant_name}_stem.json'),
                'bud': os.path.join(root, '미분할 데이터', set, 'label', 'sEEPP', self.genome_name, f'{self.plant_name}_bud.json'),
                'flower': os.path.join(root, '미분할 데이터', set, 'label', 'sEEPP', self.genome_name, f'{self.plant_name}_flower.json')
            }
        if self.model == 'ENC':
            dict = {
                'data_folder': os.path.join(root, '미분할 데이터', set, 'data', 'ENC'),
                'label_folder': os.path.join(root, '미분할 데이터', set, 'label', 'ENC'),
                'amino_acid': os.path.join(root, '미분할 데이터', set, 'data', 'ENC', genome_name, 'amino_acid.csv'),
                'label': os.path.join(root, '미분할 데이터', set, 'label', 'ENC', genome_name, f'{plant_name}.json')
            }
        return dict


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


# 파일 명명 규칙에 맞춰 결과 저장
# def save_result(data, settings, name, fname, ftype):
#     dir = os.path.join(root, settings['type'], settings['model'], name, fname+'.'+ftype)
#     if ftype == 'csv':
#         data.to_csv(dir, index=False)
#     if ftype == 'json':
#         with open(dir, 'w') as f:
#             json.dump(data, f, ensure_ascii=False, indent=4)
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
    raw_label_dir = os.path.join(root, 'ENC_label_raw')
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
