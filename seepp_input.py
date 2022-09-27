import numpy as np

from utils import *


def generate_seepp_input(target_genome_id):
    settings = {'type': 'data', 'model': 'sEEPP'}

    plant_name, genome_name = get_name(target_genome_id)

    if not file_check(settings, genome_name):
        return 0

    print('Generate sEEPP input :', genome_name, '- %d' % target_genome_id)
    df = pd.read_csv(root + '\\gene_input_raw\\' + 'Merge.%d.csv'%(target_genome_id))
    df['gene_id'] = df['PK'].map(lambda x: x[2:-1])

    create_folder(settings, genome_name)
    for col in ['promoter', 'terminator', 'utr5', 'utr3', 'cds']:
        save_result(df[['gene_id', col]], settings, genome_name, col, 'csv')
        print('%s done' % col)

    codon_usage_cols = df[df.columns[7:71]]
    gene_id_list = list(df['gene_id'])

    arr = np.array(codon_usage_cols)
    flat = arr.flatten()

    ratio = []
    gene_id = []
    codon = list(range(1, 65)) * len(df)
    for i in tqdm.tqdm(range(len(df))):
        ratio.extend(flat[64*i:64*(i+1)] / sum(flat[64*i:64*(i+1)]))
        gene_id += [gene_id_list[i]]*64

    codon_usage = pd.DataFrame({
        'gene_id': gene_id,
        'codon': codon,
        'count': flat,
        'ratio': ratio
    })

    save_result(codon_usage, settings, genome_name, 'codon_usage', 'csv')
    print('codon_usage done')


if __name__ == '__main__':
    target_genome_id = 14217
    generate_seepp_input(target_genome_id)


