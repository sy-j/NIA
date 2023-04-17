import pandas as pd

from utils import *

df = pd.read_csv(status_table)
for i in range(len(df)):
    if df['ENC dataset available'].iloc[i] == 1:
        fpf = FilePathFinder(df['genome_id'].iloc[i], 'ENC')
        print("<Preprocess ENC data/label files>")
        print(f"GenomeName: {fpf.genome_name}, GenomeID: {fpf.genome_id}")
        t = time.time()

        for set_name in ['train', 'validation', 'test']:
            print(f'  {set_name} set')
            t = time.time()
            old_path = fpf.save_path(set_name)
            new_path = fpf.save_path(set_name, purpose='실 사용 데이터')

            if not os.path.exists(os.path.join(new_path['data_folder'], f'X_{fpf.plant_name}.npz')):
                try:
                    data = pd.read_csv(old_path['amino_acid'])
                    # if len(data)>80000:
                    #     break
                    amino_onehot = []
                    for amino_acid in tqdm.tqdm(data['amino_acid']):
                        onehot = np.stack([np.eye(21)[amino_acid_list.index(amino)] for amino in amino_acid])
                        if onehot.shape[0] != maxlen['amino_acid'] or onehot.shape[1] != len(amino_acid_list):
                            print(f"Error in file {old_path['amino_acid']}\n"
                                  f"[DataShapeError] shape {onehot.shape} not compatible")
                            exit(1)
                        amino_onehot.append(onehot)
                    amino_onehot = np.stack(amino_onehot).astype(np.uint8)
                    create_folder(new_path['data_folder'])
                    np.savez_compressed(os.path.join(new_path['data_folder'], f'X_{fpf.plant_name}.npz'), data=amino_onehot)
                except:
                    print('Memory error: ', fpf.genome_name, fpf.genome_id)
            if not os.path.exists(os.path.join(new_path['label_folder'], f'Y_{fpf.plant_name}.npz')):
                with open(old_path['label'], 'r') as tmp:
                    data = json.load(tmp)
                label = data['label']
                ec_label = []
                for j in range(len(label)):
                    d = label[j]
                    ec_label.append(f"{d['ec_class']}.{d['ec_subclass']}.{d['ec_subsubclass']}.{d['ec_serial']}")
                ec_label = np.array(ec_label)
                create_folder(new_path['label_folder'])
                np.savez_compressed(os.path.join(new_path['label_folder'], f'Y_{fpf.plant_name}.npz'), data=ec_label)

            print('  time spent: %.2f(s)' % (time.time() - t))