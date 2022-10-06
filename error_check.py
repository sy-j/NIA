from utils import *


def sEEPP_error_check(target_genome_id):
    fpf = FilePathFinder(target_genome_id, 'sEEPP')
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


def ENC_error_check(target_genome_id):
    fpf = FilePathFinder(target_genome_id, 'ENC')
    print("<Inspect ENC data>")
    print(f"GenomeName: {fpf.genome_name}, GenomeID: {fpf.genome_id}")
    print('processing... ', end='')
    t = time.time()
    for dataset in ['train', 'validation', 'test']:
        path = fpf.save_path(dataset)

        # csv파일 있는지 체크
        if not os.path.exists(path['amino_acid']):
            print(f"Error in folder {os.path.basename(path['data_folder'])}\n"
                  f"[DataNotExists] no file such as {os.path.basename(path['amino_acid'])}")
            return 'DataNotExists'

        data = pd.read_csv(path['amino_acid'])
        data_len = len(data)

        # json 파일 하나 이상 존재하는지 확인, label_count, 길이 전부 정확한지 확인
        if os.path.exists(path['label']):
            with open(path['label'], 'r') as tmp:
                tmp = json.load(tmp)
                label = tmp['label']
                if tmp['label_count'] != len(tmp['label']):
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[WrongLabelCount] label_count is {tmp['label_count']}, but actually has {len(tmp['label'])} labels")
                    return 'WrongLabelCount'
                if tmp['label_count'] != data_len:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[InconsistentLabelCount] label_count is {tmp['label_count']}, but actually has {len(tmp['label'])} labels")
                    return 'InconsistentLabelCount'
        else:
            print(f"Error in folder {os.path.basename(path['label_folder'])}\n"
                  f"[LabelFileNotExists] no label(json) files")
            return 'LabelFileNotExists'

        # gene_id 전부 동일한지 체크
        for i in range(data_len):
            if data['gene_id'].iloc[i] != int(label[i]['gene_id']):
                print(f"Error in file {os.path.basename(path['amino_acid'])}, {os.path.basename(path['label'])}\n"
                      f"[GeneIDMismatch] different gene_id in in row {i}: {data['gene_id'].iloc[i]} and {label['gene_id']}")
                return 'GeneIDMismatch'

        print("no error found, time spent: %2f(s)" % (time.time()-t))
        return 0


if __name__ == '__main__':
    gid = 14225
    # df = pd.read_csv(status_table)
    # sEEPP_error_check(gid)
    ENC_error_check(gid)
    # df.to_csv(status_table, index=False)