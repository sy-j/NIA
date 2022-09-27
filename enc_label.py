from utils import *


def generate_enc_label(target_genome_id):
    settings = {'type': 'label', 'model': 'ENC'}

    plant_name, genome_name = get_name(target_genome_id)

    if not file_check(set1t
    ings, genome_name):
        return 0

    print('Generate ENC label :', genome_name, '- %d' % target_genome_id)
    df = pd.read_csv(root + '\\ENC_label_raw\\%s.csv' % genome_name)
    gene_id = []
    enzyme = []
    ec1, ec2, ec3, ec4 = [], [], [], []
    score = []

    drop_idx = []
    for i in tqdm.tqdm(range(len(df))):
        # B 포함된 결과 제거
        if 'B' in df['EC Number'].iloc[i]:
            drop_idx.append(i)
        # 4자리가 다 안나온 결과 제거
        elif len(df['EC Number'].iloc[i].split('.')) < 4:
            drop_idx.append(i)
    df = df.drop(drop_idx)
    print('%d EC removed. %d left' % (len(drop_idx), len(df)))

    input_data = pd.read_csv(root + '\\data\\ENC\\' + genome_name + '\\amino_acid.csv')
    all_gene_id_list = list(input_data.gene_id)

    df = df.sort_values(by=['Protein name', 'Score'], ascending=[True, False])
    gene_id_list = list(df['Protein name'])
    ec_list = list(df['EC Number'])
    score_list = list(df['Score'])

    max_ec = 1  # gene 하나당 최대 매칭 개수

    cnt = max_ec
    gene = gene_id_list[0]
    for i in tqdm.tqdm(range(len(df))):
        if gene_id_list[i] != gene:
            gene = gene_id_list[i]
            cnt = max_ec
        if cnt >= 1:
            gene_id.append(int(gene.replace('{P', '').replace('}', '')))
            enzyme.append(1)
            tmp = ec_list[i].split(' ')[1].split('.')
            ec1.append(int(tmp[0]))
            ec2.append(int(tmp[1]))
            ec3.append(int(tmp[2]))
            ec4.append(int(tmp[3]))
            score.append(min(round(score_list[i], 2), 400))
            cnt -= 1

    non_enzymes = list(set(all_gene_id_list) - set(gene_id))
    gene_id += non_enzymes
    enzyme += [0]*len(non_enzymes)
    ec1 += [0]*len(non_enzymes)
    ec2 += [0]*len(non_enzymes)
    ec3 += [0]*len(non_enzymes)
    ec4 += [0]*len(non_enzymes)
    score += [0]*len(non_enzymes)

    label = pd.DataFrame(gene_id)
    label.columns = ['gene_id']
    label['enzyme'] = enzyme
    label['ec_class'] = ec1
    label['ec_subclass'] = ec2
    label['ec_subsubclass'] = ec3
    label['ec_serial'] = ec4
    label['score'] = score

    label = label.sort_values(by=['gene_id', 'score'], ascending=[True, False])
    # label.to_csv(root + '\\csv\\' + genome_name + '.csv', index=False)

    if len(input_data) != len(label):
        print('Error: data and label not matched. data %d label %d' % (len(input_data), len(label)))
        return 0

    # # split 저장하기
    # for i in range((len(label)-1)//450000+1):
    #     label_split = label[i*450000:min((i+1)*450000, len(label))]
    #     dict = label_split.to_dict('records')
    #
    #     enc_label = {}
    #     enc_label['plant_name'] = plant_name
    #     enc_label['csv_file'] = '/data/ENC/' + plant_name + '/amino_acid.csv'
    #     enc_label['label_count'] = len(dict)
    #     enc_label['label'] = dict
    #
    #     create_folder(settings, genome_name)
    #     save_result(enc_label, settings, genome_name, plant_name+f'_{i+1}', 'json')

    dict = label.to_dict('records')
    enc_label = {}
    enc_label['plant_name'] = plant_name
    enc_label['csv_file'] = '/data/ENC/' + plant_name + '/amino_acid.csv'
    enc_label['label_count'] = len(dict)
    enc_label['label'] = dict
    create_folder(settings, genome_name)
    save_result(enc_label, settings, genome_name, plant_name, 'json')


if __name__ == '__main__':
    target_genome_id = [14302]
    for tgi in target_genome_id:
        generate_enc_label(tgi)