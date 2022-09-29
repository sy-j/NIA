import json

import pandas as pd
import numpy as np

from utils import *

# df = pd.read_csv(status_table)
# seepp_df = df[df['sEEPP dataset available'] == 1]
# enc_df = df[df['ENC dataset available'] == 1]
#
# print('Checking sEEPP data')
# for i in range(len(seepp_df)):
#     gid = seepp_df.genome_id.iloc[i]
#     genome_name = seepp_df.genome_name.iloc[i]
#     dirs = find_file_path(gid, genome_name)
#
#     mergefile = pd.read_csv(dirs['genome'])['PK'].map(lambda x: int(x[2:-1]))
#     tot_gene = list(mergefile)
#     if len(tot_gene) != len(mergefile.unique()):
#         print(f'Error: not unique gene id detected: check Merge.{gid}.csv')
#     print(f'checking {genome_name} - {gid}')
#     for col in ['cds', 'promoter', 'terminator', 'utr5', 'utr3', 'amino_acid']:
#         try:
#             tmp = list(pd.read_csv(dirs[col])['gene_id'])
#             if len(tmp) != len(tot_gene):
#                 print(f'Error: different number of genes: {col} have {len(tmp)} genes, but should be {len(tot_gene)}')
#             else:
#                 for j in range(len(tot_gene)):
#                     if tot_gene[j] != tmp[j]:
#                         print(f'Error: gene id mismatch: P{tot_gene[j]} might be missing in {col}')
#                         break
#         except FileNotFoundError as e:
#             print(f'Error: file missing: {genome_name} do not have {col}.csv file')
#
#     for labels in os.listdir(dirs['sEEPP_label']):
#         # tmp = pd.read_json(os.path.join(dirs['sEEPP_label'], labels), orient='records')
#         with open(os.path.join(dirs['sEEPP_label'], labels), 'r') as tmp:
#             tmp = json.load(tmp)
#         if len(tmp['organ_list']) != len(os.listdir(dirs['sEEPP_label'])):
#             print(f'Error: organ list mismatch: some json file might be missing')
#         if tmp['label_count'] != len(tmp['label']):
#             print(f'Error: incorrect label count: {os.path.basename(labels)} have {len(tmp["label"])} elements, not {tmp["label_count"]}')
#         elif tmp['label_count'] != len(tot_gene):
#             print(f'Error: different number of genes: {os.path.basename(labels)} have {tmp["label_count"]} genes, but should be {len(tot_gene)}')
#         else:
#             for j in range(len(tot_gene)):
#                 if tot_gene[j] != tmp['label'][j]['gene_id']:
#                     # print(f'Error: gene id mismatch: P{tot_gene[j]} might be missing in {labels}')
#                     break


def error_check(target_genome_id, model):
    fpf = FilePathFinder(target_genome_id, model)
    print("<Inspect sEEPP data>")
    print(f"GenomeName: {fpf.genome_name}, GenomeID: {fpf.genome_id}")
    print('processing... ', end='')
    t = time.time()
    for dataset in ['train', 'validation', 'test']:
        path = fpf.save_path(dataset)

        # csv파일 5개 전부 있는지 체크
        for file in sEEPP_input_files:
            if not os.path.exists(path[file]):
                print(f"Error in folder {os.path.basename(path['data_folder'])}\n"
                      f"[DataNotExists] no file such as {os.path.basename(path[file])}")
                return 'DataNotExists'

        data = {}
        for file in sEEPP_input_files:
            data[file] = pd.read_csv(path[file])

        # 길이 전부 동일한지 체크
        len_list = [len(data['cds']), len(data['promoter']), len(data['terminator']), len(data['utr5']), len(data['utr3']), int(len(data['codon_usage'])/64)]
        if max(len_list) != min(len_list):
            print(f"Error in folder {os.path.basename(path['data_folder'])}\n"
                  f"[InconsistentDataCount] different data counts {len_list}")
            return 'InconsistentDataCount'

        data_len = len_list[0]

        # gene_id 전부 동일한지 체크
        data_gene_id = np.array(data['cds']['gene_id'])
        for file in ['promoter', 'terminator', 'utr5', 'utr3']:
            g2 = np.array(data[file]['gene_id'])
            for i in range(data_len):
                if data_gene_id[i] != g2[i]:
                    print(f"Error in file {os.path.basename(path['cds'])}, {os.path.basename(path[file])}\n"
                          f"[InconsistentDataGeneID] different gene_id in in row {i}: {data_gene_id[i]} and {g2[i]}")
                    return 'InconsistentDataGeneID'
        g2 = np.array(data['codon_usage']['gene_id'])
        for i in range(data_len):
            if data_gene_id[i] != g2[i*64]:
                print(f"Error in file {os.path.basename(path['cds'])}, {os.path.basename(path[file])}\n"
                      f"[InconsistentDataGeneID] different gene_id in in row {i}: {data_gene_id[i]} and {g2[i*64]}")
                return 'InconsistentDataGeneID'

        # json파일 하나 이상 존재하는지 확인, label_count, 길이 전부 정확한지 확인
        label = {}
        for file in organs:
            if os.path.exists(path[file]):
                with open(path[file], 'r') as tmp:
                    tmp = json.load(tmp)
                    label[file] = tmp['label']
                    if tmp['label_count'] != len(tmp['label']):
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[WrongLabelCount] label_count is {tmp['label_count']}, but actually has {len(tmp['label'])} labels")
                        return 'WrongLabelCount'
                    if tmp['label_count'] != data_len:
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[InconsistentLabelCount] label_count is {tmp['label_count']}, but actually has {len(tmp['label'])} labels")
                        return 'InconsistentLabelCount'
        if len(label) == 0:
            print(f"Error in folder {os.path.basename(path['label_folder'])}\n"
                  f"[LabelFileNotExists] no label(json) files")
            return 'LabelFileNotExists'

        # gene_id 전부 동일한지 체크
        for file, l in label.items():
            for i in range(data_len):
                if data_gene_id[i] != int(l[i]['gene_id']):
                    print(f"Error in file {os.path.basename(path['cds'])}, {os.path.basename(path[file])}\n"
                          f"[GeneIDMismatch] different gene_id in in row {i}: {data_gene_id[i]} and {l[i]['gene_id']}")
                    return 'GeneIDMismatch'

        print("no error found, time spent: %2f(s)" % (time.time()-t))
        return 0


if __name__ == '__main__':
    gid = 14225
    # df = pd.read_csv(status_table)
    error_check(gid, 'sEEPP')
    # df.to_csv(status_table, index=False)