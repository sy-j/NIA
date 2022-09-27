from utils import *

input_cnt = 0
label_cnt = 0


def main(target_genome_id):
    global input_cnt
    global label_cnt
    df = pd.read_csv(status_table)
    tmp = df[df['genome_id'] == target_genome_id]
    # mt = pd.read_csv(matching_table)
    # tmp = mt[mt['genome_id'] == target_genome_id]
    genome_name = tmp.genome_name.iloc[0]
    history_name = tmp.history_name.iloc[0]

    # print(genome_name)

    if genome_name == 'Panicum_virgatum_5.1':
        return 0

    input_dir = root + '\\data\\ENC\\' + genome_name + '\\amino_acid.csv'
    label_dir = os.listdir(root + '\\label\\ENC\\' + genome_name)

    if os.path.exists(input_dir):
        input_df = pd.read_csv(input_dir)
        print(genome_name, '|', len(input_df))
        input_cnt += len(input_df)
    for dir in label_dir:
        if os.path.exists(dir):
            label_df = pd.read_json(dir)
            print(len(label_df['label']))
            label_cnt += len(label_df['label'])


if __name__ == '__main__':
    df = pd.read_csv(status_table)
    df = df[df['10%'] == 1]
    # df = df[df['ENC'] == 1]
    # target_genome_id = [14206]
    target_genome_id = df.genome_id
    for tgi in target_genome_id:
        main(tgi)
    print('total')
    print(input_cnt, label_cnt)