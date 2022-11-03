from utils import *

df = pd.read_csv(status_table)
for i in range(len(df)):
    if df['sEEPP file split'].iloc[i] != 1:
        if df['sEEPP dataset available'].iloc[i] == 1:
            fpf = FilePathFinder(df['genome_id'].iloc[i], 'sEEPP')
            print("<Split big sEEPP data/label files>")
            print(f"GenomeName: {fpf.genome_name}, GenomeID: {fpf.genome_id}")
            t = time.time()

            for set_name in ['train', 'validation', 'test']:
                print(f'  {set_name} set: ', end='')
                old_path = fpf.save_path(set_name)
                new_path = fpf.save_path(set_name, purpose='100MB 분할 데이터(검사용)')

                for file in sEEPP_input_files:
                    data = pd.read_csv(old_path[file])
                    if file == 'codon_usage':
                        for j in range(int(len(data)/64/20000)+1):
                            tmp = data[j*20000*64:min((j+1)*20000*64, len(data))]
                            save_result(tmp, new_path['data_folder'], '', file+f'_{j+1}', 'csv')
                    else:
                        for j in range(int(len(data)/20000)+1):
                            tmp = data[j*20000:min((j+1)*20000, len(data))]
                            save_result(tmp, new_path['data_folder'], '', file+f'_{j+1}', 'csv')
                print('data ', end='')

                for file in os.listdir(old_path['label_folder']):
                    with open(os.path.join(old_path['label_folder'], file), 'r') as tmp:
                        data = json.load(tmp)
                    for j in range(int(len(data['label'])/20000)+1):
                        label = data['label'][j*20000:min((j+1)*20000, len(data['label']))]
                        tmp = {
                            'plant_name': data['plant_name'],
                            'csv_file': data['csv_file'],
                            'organ_list': data['organ_list'],
                            'label_count': len(label),
                            'label': label,
                        }
                        save_result(tmp, new_path['label_folder'], '', file+f'_{j+1}', 'json')
                print('label')

            print("total time spent: %2f(s)" % (time.time() - t))
            df.loc[i, 'sEEPP file split'] = 1
            df.to_csv(status_table, index=False)

    if df['ENC file split'].iloc[i] != 1:
        if df['ENC dataset available'].iloc[i] == 1:
            fpf = FilePathFinder(df['genome_id'].iloc[i], 'ENC')
            print("<Split big ENC data/label files>")
            print(f"GenomeName: {fpf.genome_name}, GenomeID: {fpf.genome_id}")
            t = time.time()

            for set_name in ['train', 'validation', 'test']:
                print(f'  {set_name} set: ', end='')
                old_path = fpf.save_path(set_name)
                new_path = fpf.save_path(set_name, purpose='100MB 분할 데이터(검사용)')

                data = pd.read_csv(old_path['amino_acid'])
                for j in range(int(len(data)/50000)+1):
                    tmp = data[j*50000:min((j+1)*50000, len(data))]
                    save_result(tmp, new_path['data_folder'], '', f'amino_acid_{j+1}', 'csv')
                print('data ', end='')

                with open(old_path['label'], 'r') as tmp:
                    data = json.load(tmp)
                for j in range(int(len(data['label'])/50000)+1):
                    label = data['label'][j*50000:min((j+1)*50000, len(data['label']))]
                    tmp = {
                        'plant_name': data['plant_name'],
                        'csv_file': data['csv_file'],
                        'label_count': len(label),
                        'label': label,
                    }
                    save_result(tmp, new_path['label_folder'], '', fpf.plant_name+f'_{j+1}', 'json')
                print('label')

            print("total time spent: %2f(s)" % (time.time() - t))
            df.loc[i, 'ENC file split'] = 1
            df.to_csv(status_table, index=False)


cnt = 0
for path, dir, files in tqdm.tqdm(os.walk(os.path.join(root, '100MB 분할 데이터(검사용)'))):
    for file in files:
        cnt += 1
        file_name = os.path.join(path, file)
        file_size = os.path.getsize(file_name)
        # print(f'{round(file_size/1024/1024, 2)}MB')
        if file_size >= 1024*1024*100:
            print(f"Error in file {file_name}\n"
                  f"[FileSizeOutOfRange] file size {round(file_size/1024/1024, 2)}MB should be smaller than 100MB")
print(f"total {cnt} files")
