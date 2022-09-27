import json
from utils import *

df = pd.read_csv(status_table)
seepp_df = df[df['sEEPP dataset available'] == 1]
enc_df = df[df['ENC dataset available'] == 1]

print('Checking sEEPP data')
for i in range(len(seepp_df)):
    gid = seepp_df.genome_id.iloc[i]
    genome_name = seepp_df.genome_name.iloc[i]
    dirs = find_file_path(gid, genome_name)

    mergefile = pd.read_csv(dirs['genome'])['PK'].map(lambda x: int(x[2:-1]))
    tot_gene = list(mergefile)
    if len(tot_gene) != len(mergefile.unique()):
        print(f'Error: not unique gene id detected: check Merge.{gid}.csv')
    print(f'checking {genome_name} - {gid}')
    for col in ['cds', 'promoter', 'terminator', 'utr5', 'utr3', 'amino_acid']:
        try:
            tmp = list(pd.read_csv(dirs[col])['gene_id'])
            if len(tmp) != len(tot_gene):
                print(f'Error: different number of genes: {col} have {len(tmp)} genes, but should be {len(tot_gene)}')
            else:
                for j in range(len(tot_gene)):
                    if tot_gene[j] != tmp[j]:
                        print(f'Error: gene id mismatch: P{tot_gene[j]} might be missing in {col}')
                        break
        except FileNotFoundError as e:
            print(f'Error: file missing: {genome_name} do not have {col}.csv file')

    for labels in os.listdir(dirs['sEEPP_label']):
        # tmp = pd.read_json(os.path.join(dirs['sEEPP_label'], labels), orient='records')
        with open(os.path.join(dirs['sEEPP_label'], labels), 'r') as tmp:
            tmp = json.load(tmp)
        if len(tmp['organ_list']) != len(os.listdir(dirs['sEEPP_label'])):
            print(f'Error: organ list mismatch: some json file might be missing')
        if tmp['label_count'] != len(tmp['label']):
            print(f'Error: incorrect label count: {os.path.basename(labels)} have {len(tmp["label"])} elements, not {tmp["label_count"]}')
        elif tmp['label_count'] != len(tot_gene):
            print(f'Error: different number of genes: {os.path.basename(labels)} have {tmp["label_count"]} genes, but should be {len(tot_gene)}')
        else:
            for j in range(len(tot_gene)):
                if tot_gene[j] != tmp['label'][j]['gene_id']:
                    # print(f'Error: gene id mismatch: P{tot_gene[j]} might be missing in {labels}')
                    break


# for i in range(len(enc_df)):
#     gid = enc_df.genome_id.iloc[i]
#     genome_name = enc_df.genome_name.iloc[i]
#     dirs = find_file_path(gid, genome_name)
#
#     mergefile = pd.read_csv(dirs['genome'])['PK'].map(lambda x: int(x[2:-1]))
#     tot_gene = list(mergefile)
#     if len(tot_gene) != len(mergefile.unique()):
#         print(f'Error: not unique gene id detected: check Merge.{gid}.csv')
#     print(f'checking {genome_name} - {gid}')
#     for col in ['amino_acid']:
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