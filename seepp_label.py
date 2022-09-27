from utils import *


def generate_seepp_label(target_genome_id):
    settings = {'type': 'label', 'model': 'sEEPP'}

    plant_name, genome_name = get_name(target_genome_id)

    if not file_check(settings, genome_name):
        return 0

    print('Generate sEEPP label :', genome_name, '- %d' % target_genome_id)

    df = pd.read_csv(root + '\\sEEPP_label_raw\\Total_%d.csv' % target_genome_id)
    df.dropna(axis=1, how='all', inplace=True)
    col_list = list(df.columns)
    col_list.remove('gene_id')
    for col in col_list:
        temp_df = df.loc[:, ['gene_id', col]]
        temp_df['organ'] = [col]*len(temp_df)
        temp_df.columns = ['gene_id', 'expression_valid', 'organ']
        temp_df = temp_df[['gene_id', 'organ', 'expression_valid']]
        temp_dict = temp_df.to_dict(orient='records')
        seepp_label = {}
        seepp_label['plant_name'] = plant_name
        seepp_label['csv_file'] = [f"data/sEEPP/{plant_name}/cds.csv",
                                f"data/sEEPP/{plant_name}/codon_usage.csv",
                                f"data/sEEPP/{plant_name}/promoter.csv",
                                f"data/sEEPP/{plant_name}/terminator.csv",
                                f"data/sEEPP/{plant_name}/utr3.csv",
                                f"data/sEEPP/{plant_name}/utr5.csv"
                                ]
        seepp_label['organ_list'] = col_list
        seepp_label['label_count'] = len(temp_df)
        seepp_label['label'] = temp_dict

        create_folder(settings, genome_name)
        save_result(seepp_label, settings, genome_name, plant_name+f'_{col}', 'json')


if __name__ == '__main__':
    target_genome_id = [14223]
    for tgi in target_genome_id:
        generate_seepp_label(tgi)