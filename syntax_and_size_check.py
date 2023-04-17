import os

from utils import *

df = pd.read_csv(status_table3)


def sEEPP_error_check(idx):
    fpf = FilePathFinder(df['genome_id'].iloc[idx], 'sEEPP')
    print("<Checking sEEPP data>")
    print(f"GenomeName: {fpf.genome_name}, GenomeID: {fpf.genome_id}")
    print('processing... ', end='')
    t = time.time()
    for dataset in ['train', 'validation', 'test']:
        path = fpf.save_path(dataset, purpose='data_all')

        # csv파일 5개 전부 있는지, row 수 동일한지, 이외 파일 없는지 확인
        if 1:
            for file in sEEPP_input_files:
                if not os.path.exists(path[file]):
                    print(f"Error in folder {os.path.basename(path['data_folder'])}\n"
                          f"[sEEPP:Err1] no file such as {os.path.basename(path[file])}")
                    return 'sEEPP:Err1'
            csv_list = os.listdir(path['data_folder'])
            for csv_file in csv_list:
                if csv_file.split('.')[0] not in sEEPP_input_files:
                    print(f"Error in folder {os.path.basename(path['data_folder'])}\n"
                          f"[sEEPP:Err2] unknown data file {csv_file}")
                    return 'sEEPP:Err2'

            data = {}
            for file in sEEPP_input_files:
                data[file] = pd.read_csv(path[file])

            data_len = [len(data['cds']), len(data['promoter']), len(data['terminator']), len(data['utr5']), len(data['utr3']), int(len(data['codon_usage'])/64)]
            if max(data_len) != min(data_len):
                print(f"Error in folder {os.path.basename(path['data_folder'])}\n"
                      f"[sEEPP:Err3] different data counts {data_len}")
                return 'sEEPP:Err3'
            data_len = data_len[0]

        # 원천데이터 구문정확성
        if 1:
            for file in nucleotide_seq_cols:  # codon usage 제외(전부 string 임)
                gene_id_list = list(data[file]['gene_id'])
                for gene_id in gene_id_list:
                    # if type(gene_id) != str:  # 지금 전부다 int임
                    #     print(f"Error in file amino_acid.csv\n"
                    #           f"[Syntax1] gene_id should be string type: {gene_id}, {type(gene_id)}")
                    #     return 'Syntax1'
                    if int(gene_id) < 0:
                        print(f"Error in file {file}.csv\n"
                              f"[sEEPP:Syntax2] gene_id should be positive value: {gene_id}")
                        return 'sEEPP:Syntax2'
                seq_list = list(data[file][file])
                for seq in seq_list:
                    if type(seq) != str:
                        print(f"Error in file {file}.csv\n"
                              f"[sEEPP:Syntax3] {file} should be string type: {seq}")
                        return 'sEEPP:Syntax3'
                    if len(seq) != maxlen[file]:
                        print(f"Error in file {file}.csv\n"
                              f"[sEEPP:Syntax3] invalid {file} length: {seq}")
                        return 'sEEPP:Syntax3'

        # json파일 하나 이상 존재하는지 확인, 이외 파일 없는지 확인
        if 1:
            js = {}
            label = {}
            organ_list = []
            for file in organs:
                if os.path.exists(path[file]):
                    with open(path[file], 'r') as tmp:
                        js[file] = json.load(tmp)
                        label[file] = js[file]['label']
                        organ_list.append(file)
            if len(label) == 0:
                print(f"Error in folder {os.path.basename(path['label_folder'])}\n"
                      f"[sEEPP:Err4] no label(json) files")
                return 'sEEPP:Err4'
            if len(label) != len(os.listdir(path['label_folder'])) or len(label) != len(organ_list):
                print(f"Error in folder {os.path.basename(path['label_folder'])}\n"
                      f"[sEEPP:Err5] unknown label file")
                return 'sEEPP:Err5'

        # 라벨링데이터 구문정확성
        if 1:
            for file in organ_list:
                if type(js[file]['plant_name']) != str:
                    print(f"Error in file {os.path.basename(path[file])}\n"
                          f"[sEEPP:Syntax4] plant_name should be string type: {js[file]['plant_name']}")
                    return 'sEEPP:Syntax4'
                if type(js[file]['csv_file']) != list:
                    print(f"Error in file {os.path.basename(path[file])}\n"
                          f"[sEEPP:Syntax5] csv_file should be array type: {type(js[file]['csv_file'])}")
                    return 'sEEPP:Syntax5'
                if len(js[file]['csv_file']) != 6:
                    print(f"Error in file {os.path.basename(path[file])}\n"
                          f"[sEEPP:Syntax17] csv_file should contain 6 values: {js[file]['csv_file']}")
                    return 'sEEPP:Syntax17'
                for value in js[file]['csv_file']:
                    if type(value) != str:
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[sEEPP:Syntax6] csv_file value should be string type: {value}")
                        return 'sEEPP:Syntax6'
                if type(js[file]['organ_list']) != list:
                    print(f"Error in file {os.path.basename(path[file])}\n"
                          f"[sEEPP:Syntax7] organ_list should be array type: {type(js[file]['organ_list'])}")
                    return 'sEEPP:Syntax7'
                if len(organ_list) < 1 or len(organ_list) > 5:
                    print(f"Error in file {os.path.basename(path[file])}\n"
                          f"[sEEPP:Syntax18] invalid organ length: {organ_list}")
                    return 'sEEPP:Syntax18'
                for value in js[file]['organ_list']:
                    if type(value) != str:
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[sEEPP:Syntax8] organ_list value should be string type: {value}")
                        return 'sEEPP:Syntax8'
                    if value not in organs:
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[sEEPP:Syntax9] invalid organ value: {value}")
                        return 'sEEPP:Syntax9'
                if type(js[file]['label_count']) != int:
                    print(f"Error in file {os.path.basename(path[file])}\n"
                          f"[sEEPP:Syntax10] label_count should be int type: {js[file]['label_count']}")
                    return 'sEEPP:Syntax10'
                if type(label[file]) != list:
                    print(f"Error in file {os.path.basename(path[file])}\n"
                          f"[sEEPP:Syntax11] label should be array type: {type(label[file])}")
                    return 'sEEPP:Syntax11'
                for i in range(len(label[file])):
                    gene_id = label[file][i]['gene_id']
                    organ = label[file][i]['organ']
                    expression_valid = label[file][i]['expression_valid']
                    if type(gene_id) != int:
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[sEEPP:Syntax12] gene_id should be int type: {gene_id}")
                        return 'sEEPP:Syntax12'
                    if type(organ) != str:
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[sEEPP:Syntax13] organ should be int type: {organ}")
                        return 'sEEPP:Syntax13'
                    if type(expression_valid) != int:
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[sEEPP:Syntax14] expression_valid should be int type: {expression_valid}")
                        return 'sEEPP:Syntax14'
                    if organ not in organs:
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[sEEPP:Syntax15] invalid organ value: {organ}")
                        return 'sEEPP:Syntax15'
                    if expression_valid not in [0, 1]:
                        print(f"Error in file {os.path.basename(path[file])}\n"
                              f"[sEEPP:Syntax16] invalid expression_valid value: {expression_valid}")
                        return 'sEEPP:Syntax16'

        # json 내용 정확한지 체크
        for file in organ_list:
            if js[file]['plant_name'] != fpf.plant_name:
                print(f"Error in file {os.path.basename(path[file])}\n"
                      f"[sEEPP:Err6] plant_name is {js[file]['plant_name']}, but should be {fpf.plant_name}")
                return 'sEEPP:Err6'
            if len(js[file]['organ_list']) != len(organ_list) or set(js[file]['organ_list']) != set(organ_list):
                print(f"Error in file {os.path.basename(path[file])}\n"
                      f"[sEEPP:Err7] organ_list mismatch occurred: {js[file]['organ_list']}, {organ_list}")
                return 'sEEPP:Err7'
            if js[file]['label_count'] != len(label[file]):
                print(f"Error in file {os.path.basename(path[file])}\n"
                      f"[sEEPP:Err8] label_count is {js[file]['label_count']}, but actually has {label[file]} labels")
                return 'sEEPP:Err8'
            if js[file]['label_count'] != data_len:
                print(f"Error in file {os.path.basename(path[file])}\n"
                      f"[sEEPP:Err9] label_count is {js[file]['label_count']}, but should be {data_len} labels")
                return 'sEEPP:Err9'

        # 원천간 gene_id 전부 동일한지 체크
        for file in nucleotide_seq_cols:
            g2 = list(data[file]['gene_id'])
            for i in range(data_len):
                if gene_id_list[i] != g2[i]:
                    print(f"Error in file {os.path.basename(path['cds'])}, {os.path.basename(path[file])}\n"
                          f"[sEEPP:Err10] different gene_id in in row {i}: {gene_id_list[i]} and {g2[i]}")
                    return 'sEEPP:Err10'
        g2 = list(data['codon_usage']['gene_id'])
        for i in range(data_len):
            if gene_id_list[i] != g2[i*64]:
                print(f"Error in file {os.path.basename(path['cds'])}, {os.path.basename(path['codon_usage'])}\n"
                      f"[sEEPP:Err11] different gene_id in in row {i}: {gene_id_list[i]} and {g2[i*64]}")
                return 'sEEPP:Err11'

        # 원천, 라벨간 gene_id 전부 동일한지 체크
        for file, l in label.items():
            for i in range(data_len):
                if gene_id_list[i] != int(l[i]['gene_id']):
                    print(f"Error in file {os.path.basename(path['cds'])}, {os.path.basename(path[file])}\n"
                          f"[sEEPP:Err12] different gene_id in in row {i}: {gene_id_list[i]} and {l[i]['gene_id']}")
                    return 'sEEPP:Err12'

        df.loc[idx, f'sEEPP {dataset}'] = data_len
        df.to_csv(status_table3, index=False)

    df.loc[idx, 'organs'] = len(organ_list)
    df.loc[idx, 'sEEPP input'] = df.iloc[idx]['sEEPP train'] + df.iloc[idx]['sEEPP validation'] + df.iloc[idx]['sEEPP test']
    df.loc[idx, 'sEEPP label'] = df.loc[idx, 'sEEPP input'] * df.loc[idx, 'organs']
    print("no error found, time spent: %2f(s)" % (time.time() - t))
    return 0


def ENC_error_check(idx):
    fpf = FilePathFinder(df['genome_id'].iloc[idx], 'ENC')
    print("<Check ENC data>")
    print(f"GenomeName: {fpf.genome_name}, GenomeID: {fpf.genome_id}")
    print('processing... ', end='')
    t = time.time()
    for dataset in ['train', 'validation', 'test']:
        path = fpf.save_path(dataset, purpose='data_all')

        # amino acid 파일 존재, 유일한지 확인
        if 1:
            if not os.path.exists(path['amino_acid']):
                print(f"Error in folder {os.path.basename(path['data_folder'])}\n"
                      f"[ENCErr1] no file such as amino_acid.csv")
                return 'ENCErr1'
            csv_list = os.listdir(path['data_folder'])
            if len(csv_list) > 1:
                for csv_file in csv_list:
                    if csv_file != 'amino_acid.csv':
                        print(f"Error in folder {os.path.basename(path['data_folder'])}\n"
                              f"[ENCErr2] unknown data file {csv_file}")
                        return 'ENCErr2'

            data = pd.read_csv(path['amino_acid'])
            data_len = len(data)

        # 원천데이터 구문정확성
        if 1:
            gene_id_list = list(data['gene_id'])
            for gene_id in gene_id_list:
                # if type(gene_id) != str:  # 지금 전부다 int임
                #     print(f"Error in file amino_acid.csv\n"
                #           f"[Syntax1] gene_id should be string type: {gene_id}, {type(gene_id)}")
                #     return 'Syntax1'
                if int(gene_id) < 0:
                    print(f"Error in file amino_acid.csv\n"
                          f"[Syntax2] gene_id should be positive value: {gene_id}")
                    return 'Syntax2'
            seq_list = list(data['amino_acid'])
            for seq in seq_list:
                if type(seq) != str:
                    print(f"Error in file amino_acid.csv\n"
                          f"[Syntax3] amino_acid should be string type: {seq}")
                    return 'Syntax3'
                if len(seq) < 1:
                    print(f"Error in file amino_acid.csv\n"
                          f"[Syntax3] amino_acid should longer than 0: {seq}")
                    return 'Syntax3'

        # label json 존재, 유일한지 확인
        if 1:
            if not os.path.exists(path['label']):
                print(f"Error in folder {os.path.basename(path['label_folder'])}\n"
                      f"[ENCErr3] no file such as {os.path.basename(path['label'])}")
                return 'ENCErr3'
            json_list = os.listdir(path['label_folder'])
            if len(json_list) > 1:
                for json_file in json_list:
                    if json_file != os.path.basename(path['label']):
                        print(f"Error in folder {os.path.basename(path['label_folder'])}\n"
                              f"[ENCErr4] unknown label file {json_file}")
                        return 'ENCErr4'

            with open(path['label'], 'r') as js:
                js = json.load(js)
            label = js['label']

        # 라벨링데이터 구문정확성
        if 1:
            if type(js['plant_name']) != str:
                print(f"Error in file {os.path.basename(path['label'])}\n"
                      f"[Syntax4] plant_name should be string type: {js['plant_name']}")
                return 'Syntax4'
            if type(js['csv_file']) != str:
                print(f"Error in file {os.path.basename(path['label'])}\n"
                      f"[Syntax5] csv_file should be string type: {js['csv_file']}")
                return 'Syntax5'
            if type(js['label_count']) != int:
                print(f"Error in file {os.path.basename(path['label'])}\n"
                      f"[Syntax6] label_count should be int type: {js['label_count']}")
                return 'Syntax6'
            if type(label) != list:
                print(f"Error in file {os.path.basename(path['label'])}\n"
                      f"[Syntax7] label should be array type: {type(label)}")
                return 'Syntax7'
            for i in range(len(label)):
                gene_id = label[i]['gene_id']
                enzyme = label[i]['enzyme']
                ec_class = label[i]['ec_class']
                ec_subclass = label[i]['ec_subclass']
                ec_subsubclass = label[i]['ec_subsubclass']
                ec_serial = label[i]['ec_serial']
                score = label[i]['score']
                if type(gene_id) != int:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax8] gene_id should be int type: {gene_id}")
                    return 'Syntax8'
                if type(enzyme) != int:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax9] enzyme should be int type: {enzyme}")
                    return 'Syntax9'
                if type(ec_class) != int:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax10] ec_class should be int type: {ec_class}")
                    return 'Syntax10'
                if type(ec_subclass) != int:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax11] ec_subclass should be int type: {ec_subclass}")
                    return 'Syntax11'
                if type(ec_subsubclass) != int:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax12] ec_subsubclass should be int type: {ec_subsubclass}")
                    return 'Syntax12'
                if type(ec_serial) != int:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax13] ec_serial should be int type: {ec_serial}")
                    return 'Syntax13'
                if type(score) != float:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax14] score should be int type: {score}")
                    return 'Syntax14'
                if enzyme not in [0, 1]:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax15] invalid enzyme value: {enzyme}")
                    return 'Syntax15'
                if ec_class < 0 or ec_class > 7:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax16] invalid ec_class value: {ec_class}")
                    return 'Syntax16'
                if ec_subclass < 0 or ec_subclass > 99:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax17] invalid ec_subclass value: {ec_subclass}")
                    return 'Syntax17'
                if ec_subsubclass < 0 or ec_subsubclass > 99:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax18] invalid ec_subsubclass value: {ec_subsubclass}")
                    return 'Syntax18'
                if ec_serial < 0 or ec_serial > 999:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax19] invalid ec_serial value: {ec_serial}")
                    return 'Syntax19'
                if score < 0 or score > 400:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax20] invalid score value: {score}")
                    return 'Syntax20'
                if len(str(score).split('.')[-1]) > 2:
                    print(f"Error in file {os.path.basename(path['label'])}\n"
                          f"[Syntax21] invalid score decimal point: {score}")
                    return 'Syntax21'

        # json 내용 정확한지 체크
        if js['plant_name'] != fpf.plant_name:
            print(f"Error in file {os.path.basename(path['label'])}\n"
                  f"[ENCErr5] plant_name is {js['plant_name']}, but should be {fpf.plant_name}")
            return 'ENCErr5'
        # if js['csv_file'] !=
        if js['label_count'] != len(label):
            print(f"Error in file {os.path.basename(path['label'])}\n"
                  f"[ENCErr6] label_count is {js['label_count']}, but actually has {len(label)} labels")
            return 'ENCErr6'
        if js['label_count'] != data_len:
            print(f"Error in file {os.path.basename(path['label'])}\n"
                  f"[ENCErr7] label_count is {js['label_count']}, but should be {data_len}")
            return 'ENCErr7'

        for i in range(data_len):
            if gene_id_list[i] != label[i]['gene_id']:
                print(f"Error in file {os.path.basename(path['amino_acid'])}, {os.path.basename(path['label'])}\n"
                      f"[GeneIDMismatch] different gene_id in in row {i}: {data['gene_id'].iloc[i]} and {label['gene_id']}")
                return 'GeneIDMismatch'
        # gene_id 전부 동일한지 체크
        # score 값 유효범위 이내인지 체크
        for i in range(data_len):
            if data['gene_id'].iloc[i] != int(label[i]['gene_id']):
                print(f"Error in file {os.path.basename(path['amino_acid'])}, {os.path.basename(path['label'])}\n"
                      f"[GeneIDMismatch] different gene_id in in row {i}: {data['gene_id'].iloc[i]} and {label['gene_id']}")
                return 'GeneIDMismatch'
            if 0 < float(label[i]['score']) < enc_score_threshold or float(label[i]['score']) > enc_max_score:
                print(f"Error in file {os.path.basename(path['label'])}\n"
                      f"[ScoreOutOfRange] score out of range in row {i}: {label[i]['score']}")
                return 'ScoreOutOfRange'

        df.loc[idx, f'ENC {dataset}'] = data_len
        df.to_csv(status_table3, index=False)

    df.loc[idx, 'ENC input'] = df.iloc[idx]['ENC train'] + df.iloc[idx]['ENC validation'] + df.iloc[idx]['ENC test']
    df.loc[idx, 'ENC label'] = df.loc[idx, 'ENC input']
    print("no error found, time spent: %2f(s)" % (time.time()-t))
    return 0


def check_all_data():
    print('Checking all data...')
    t = time.time()
    for i in range(len(df)):
        if df['sEEPP error'].iloc[i] != '0' and df['sEEPP error'].iloc[i] != 0 and df['sEEPP'].iloc[i] == 1:
            df.loc[i, 'sEEPP error'] = sEEPP_error_check(i)
        if df['ENC error'].iloc[i] != '0' and df['ENC error'].iloc[i] != 0:
            df.loc[i, 'ENC error'] = ENC_error_check(i)
    df.to_csv(status_table3, index=False)
    print("all process done, total time spent: %2f(s)" % (time.time() - t))


def ENC_fix(idx):
    error_code = df['ENC error'].iloc[idx]
    fpf = FilePathFinder(df['genome_id'].iloc[idx], 'ENC')
    print("<Update ENC data>")
    print(f"GenomeName: {fpf.genome_name}, GenomeID: {fpf.genome_id}, ErrorID: {error_code}")
    print('processing... ', end='')
    t = time.time()

    if error_code == 'Syntax21':
        for dataset in ['train', 'validation', 'test']:
            path = fpf.save_path(dataset, purpose='data_all')
            with open(path['label'], 'r') as js:
                js = json.load(js)
            for i in range(len(js['label'])):
                if len(str(js['label'][i]['score']).split('.')[-1]) > 2:
                    js['label'][i]['score'] = np.round(js['label'][i]['score'], 2)
            with open(path['label'], 'w') as f:
                json.dump(js, f, ensure_ascii=False, indent=4)

    # ENC_error_check(idx)

    print("done, time spent: %2f(s)" % (time.time() - t))
    return 0


def sEEPP_fix(idx):
    error_code = df['sEEPP error'].iloc[idx]
    fpf = FilePathFinder(df['genome_id'].iloc[idx], 'sEEPP')
    print("<Update sEEPP data>")
    print(f"GenomeName: {fpf.genome_name}, GenomeID: {fpf.genome_id}, ErrorID: {error_code}")
    print('processing... ', end='')
    t = time.time()

    if error_code == 'sEEPP:Err7':
        for dataset in ['train', 'validation', 'test']:
            path = fpf.save_path(dataset, purpose='data_all')

            js = {}
            organ_list = []
            for file in organs:
                if os.path.exists(path[file]):
                    with open(path[file], 'r') as tmp:
                        js[file] = json.load(tmp)
                        organ_list.append(file)

            for file in organ_list:
                js[file]['organ_list'] = organ_list
                with open(path[file], 'w') as f:
                    json.dump(js[file], f, ensure_ascii=False, indent=4)

    # sEEPP_error_check(idx)

    print("done, time spent: %2f(s)" % (time.time() - t))
    return 0


def fix_all_data():
    print('Updating all data...')
    t = time.time()
    for i in range(len(df)):
        if df['sEEPP error'].iloc[i] != '0' and df['sEEPP error'].iloc[i] != 0 and df['sEEPP'].iloc[i] == 1:
            sEEPP_fix(i)
        if df['ENC error'].iloc[i] != '0' and df['ENC error'].iloc[i] != 0:
            ENC_fix(i)
    df.to_csv(status_table3, index=False)

    print("all process done, total time spent: %2f(s)" % (time.time() - t))


if __name__ == '__main__':
    i = 0
    # df.loc[i, 'sEEPP error'] = sEEPP_error_check(i)
    # df.loc[i, 'ENC error'] = ENC_error_check(i)
    # df.to_csv(status_table3, index=False)
    #
    # if df['ENC error'].iloc[i] != 0:
    #     ENC_fix(i)

    # fix_all_data()
    # check_all_data()

    # sEEPP_fix(95)
    ENC_error_check(95)