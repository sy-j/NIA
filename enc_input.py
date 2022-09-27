from utils import *


def generate_enc_input(target_genome_id):
    settings = {'type': 'data', 'model': 'ENC'}

    plant_name, genome_name = get_name(target_genome_id)

    # 중복 덮어쓰기 체크
    if not file_check(settings, genome_name):
        return 0

    print('Generate ENC input :', genome_name, '- %d' % target_genome_id)
    df = pd.read_csv(root + '\\gene_input_raw\\' + 'Merge.%d.csv' % target_genome_id)
    df['PK'] = df['PK'].map(lambda x: x[2:-1])

    gene_id = []
    amino_acid = []

    for i in tqdm.tqdm(range(len(df))):
        gene_id.append(df.PK.iloc[i])

        # 아미노산에서 ']' 제거
        amino_acid.append(df.amino_acid.iloc[i].replace(']', ''))

    df['gene_id'] = gene_id
    df['amino_acid'] = amino_acid
    enc_input = df[['gene_id', 'amino_acid']]

    create_folder(settings, genome_name)
    save_result(enc_input, settings, genome_name, 'amino_acid', 'csv')


if __name__ == '__main__':
    target_genome_id = [14302]
    for tgi in target_genome_id:
        generate_enc_input(tgi)
