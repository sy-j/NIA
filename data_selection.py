import pandas as pd

from utils import *
import status_check

p = r'D:\NIA\학습용 데이터\개수 조정된 원데이터'

# '''
for path, dir, files in os.walk(p):
    for file in files:
        if file[-8:] == 'ec10.csv' or file[-8:] == 'crop.csv' or file[-11:] == 'dropped.csv':
            continue
        if os.path.exists(os.path.join(path, file[:-4]+'_crop.csv')):
            continue
        print(file)
        data = pd.read_csv(os.path.join(path, file))

        print(f'data size: {len(data)}', end='')
        non_enzyme = data[data['EC Number'] == 'EC 0.0.0.0']
        enzyme = data[data['EC Number'] != 'EC 0.0.0.0']

        if len(non_enzyme)/(len(non_enzyme) + len(enzyme)) >= 0.9:
            if len(enzyme) + 0.4 * len(non_enzyme) > 100000:  # 기공개 30000
            # if len(enzyme) + 0.4 * len(non_enzyme) > 30000:  # 기공개 30000
                non_enzyme = non_enzyme.sample(frac=0.4, random_state=911)

        if len(non_enzyme)/(len(non_enzyme) + len(enzyme)) >= 0.8:
            if len(enzyme) + 0.5 * len(non_enzyme) > 100000:  # 기공개 30000
            # if len(enzyme) + 0.5 * len(non_enzyme) > 30000:  # 기공개 30000
                non_enzyme = non_enzyme.sample(frac=0.5, random_state=911)

        data = pd.concat([non_enzyme, enzyme])
        print(f' -> {len(data)}', end='')

        # 기공개만
        # if len(data) > 70000:  # 기공개 70000
        #     data = data.sample(n=70000, random_state=911)    # 기공개 70000
        # print(f' -> {len(data)}', end='')

        ec_count = data.value_counts(['EC Number'])
        ec_list = data['EC Number'].unique()
        drop_idx = []
        for ec in ec_list:
            if ec_count[ec] <= 10:
                drop_idx += list(data[data['EC Number'] == ec].index)
        if len(data)-len(drop_idx) > 100000:   # 기공개 30000
        # if len(data)-len(drop_idx) > 30000:   # 기공개 30000
            data = data.drop(drop_idx)

        print(f' -> {len(data)}', end='')

        if len(data) > 100000:    # 기공개 60000
        # if len(data) > 60000:    # 기공개 60000
            val_cnt = data.value_counts(['EC Number'])
            # print(val_cnt)
            ec_count = list(val_cnt)
            # print(ec_count)
            for i in range(len(ec_count)):
                if len(data) - sum(ec_count[-i:]) <= 100000:    # 기공개 60000
                # if len(data) - sum(ec_count[-i:]) <= 60000:    # 기공개 60000
                    tmp = i+1
                    break
            ec_list = [val_cnt.index[i][0] for i in range(len(val_cnt))][-tmp:]
            # print(ec_list)
            drop_idx = []
            for ec in ec_list:
                drop_idx += list(data[data['EC Number'] == ec].index)
            data = data.drop(drop_idx)

        # 신규만
        if len(data) > 100000:
            data = data.sample(n=100000, random_state=911)
        print(f' -> {len(data)}', end='')

        data = data.sort_values(by=['gene_id'])
        non_enzyme_rate = round(data.value_counts(['EC Number'])[0]/len(data)*100, 2)
        print(f' -> {len(data)}.')
        print(f"non-enzyme: {non_enzyme_rate}%, unique EC: {data['EC Number'].nunique()}")
        data.to_csv(os.path.join(path, file[:-4]+'_crop.csv'), index=False)
        df = pd.read_csv(status_table2)
        df.loc[df[df['genome_id'] == int(file[:-4])].index, 'crop'] = len(data)
        df.loc[df[df['genome_id'] == int(file[:-4])].index, 'crop non-enzyme'] = non_enzyme_rate
        df.to_csv(status_table2, index=False)


'''
# total, merge, BAresult 기본전처리후 합쳐서 저장
df = pd.read_csv(status_table)
rename_annotation_result()
for i in range(len(df)):
# for i in [142]:
    target = [14578, 14550, 14552]
    if df['genome_id'].iloc[i] not in target:
        continue
    if df['ENC'].iloc[i]:
        if df['sEEPP'].iloc[i]:
            if df['genome download done'].iloc[i] and df['Cufflink result download done'].iloc[i] and df['BRENDA result download done'].iloc[i] \
                    and not os.path.exists(os.path.join(root, '개수 조정된 원데이터', 'sEEPP_ENC', f"{df['genome_id'].iloc[i]}.csv")):
                print(f"<Preprocessing> GenomeName: {df['genome_name'].iloc[i]}, GenomeID: {df['genome_id'].iloc[i]}, target: sEEPP & ENC")
                fpf_sEEPP = FilePathFinder(df['genome_id'].iloc[i], 'sEEPP')
                fpf_ENC = FilePathFinder(df['genome_id'].iloc[i], 'ENC')

                # Merge.XXXXX 파일 읽고 PK->gene_id 변환
                print('phase 1: ', end='')
                t = time.time()
                data_df = pd.read_csv(fpf_sEEPP.read_path()['data'])
                data_df['PK'] = data_df['PK'].map(lambda x: x[2:-1])
                data_df.rename(columns={'PK': 'gene_id'}, inplace=True)
                print('done, time spent = %.2f(s)' % (time.time()-t))
                prev_len = len(data_df)
                codon_cols = list(data_df.columns[7:71])

                # ATGCN 이외 문자를 갖는 서열 제거
                print('phase 2: ', end='')
                t = time.time()
                drop_idx = []
                for j in range(len(data_df)):
                    seq = data_df.iloc[j]
                    for col in nucleotide_seq_cols:
                        if str(seq[col]) == 'nan':
                            continue
                        elif len(set(seq[col]) - set('ATGCN')) > 0:
                            # print(set(seq[col]))
                            drop_idx.append(j)
                            break
                data_df = data_df.drop(drop_idx)
                data_df = data_df.reset_index(drop=True)
                print('done, time spent = %.2f(s), %d row removed' % (time.time() - t, len(drop_idx)))

                # 아미노산 서열 길이 10 이하, 21개 문자 이외 포함하는 서열 제거
                # 단, 마지막에 ] 들어있는 것은 추가 파싱후 사용, 추후 버그 수정시 해당 기능도 삭제
                print('phase 3: ', end='')
                t = time.time()
                drop_idx = []
                amino_acid = list(data_df['amino_acid'])
                for j in range(len(amino_acid)):
                    seq = amino_acid[j]
                    if len(str(seq)) < 10:
                        drop_idx.append(j)
                    elif len(set(seq) - set('ACDEFGHIKLMNPQRSTVWXY')) > 0:
                        if len(set(seq[:-1]) - set('ACDEFGHIKLMNPQRSTVWXY')) == 0:
                            amino_acid[j] = seq[:-1]
                        else:
                            drop_idx.append(j)
                data_df['amino_acid'] = amino_acid
                data_df = data_df.drop(drop_idx)
                print('done, time spent = %.2f(s), %d row removed' % (time.time() - t, len(drop_idx)))

                # 사전에 정의된 길이로 고정
                print('phase 4: ', end='')
                t = time.time()
                for col in ['cds', 'terminator', 'utr3']:
                    data_df[col] = data_df[col].map(
                        lambda x: 'N' * maxlen[col] if type(x) != str
                        else (x[:maxlen[col]] if len(x) > maxlen[col] else x + 'N' * (maxlen[col] - len(x))))
                for col in ['promoter', 'utr5']:
                    data_df[col] = data_df[col].map(
                        lambda x: 'N' * maxlen[col] if type(x) != str
                        else (x[-maxlen[col]:] if len(x) > maxlen[col] else 'N' * (maxlen[col] - len(x)) + x))
                data_df['amino_acid'] = data_df['amino_acid'].map(
                    lambda x: 'N' * maxlen['amino_acid'] if type(x) != str
                    else (x[:maxlen['amino_acid']] if len(x) > maxlen['amino_acid'] else x + 'N' * (
                                maxlen['amino_acid'] - len(x))))
                # amino_acid 중복 제거, [promoter terminator 5utr 3utr cds] 중복 제거
                data_df = data_df.drop_duplicates(subset=['amino_acid'])
                data_df = data_df.drop_duplicates(subset=nucleotide_seq_cols)
                print('done, time spent = %.2f(s)' % (time.time() - t))

                # Total_XXXXX 파일 읽고 gene_id 변환
                print('phase 5: ', end='')
                t = time.time()
                sEEPP_label_df = pd.read_csv(fpf_sEEPP.read_path()['label'])
                # sEEPP_label_df.dropna(axis=1, how='all', inplace=True)
                sEEPP_label_df['gene_id'] = sEEPP_label_df['gene_id'].map(lambda x: x[2:-1])
                # num_conditions = len(list(sEEPP_label_df.columns))-1

                data_df = data_df[['gene_id'] + ['amino_acid'] + nucleotide_seq_cols + codon_cols]
                data_df = pd.merge(data_df, sEEPP_label_df, left_on='gene_id', right_on='gene_id', how='inner')

                del sEEPP_label_df
                print('done, time spent = %.2f(s)' % (time.time() - t))

                # BRENDA annotation result 파일 읽고 Protein name->gene_id 변환
                print('phase 6: ', end='')
                t = time.time()
                ENC_label_df = pd.read_csv(fpf_ENC.read_path()['label'])
                ENC_label_df['Protein name'] = ENC_label_df['Protein name'].map(lambda x: x[2:-1])
                ENC_label_df.rename(columns={'Protein name': 'gene_id'}, inplace=True)
                print('done, time spent = %.2f(s)' % (time.time() - t))

                # 클래스로 정의되지 않는 데이터 제거 (1~3자리 EC, 또는 B 포함하는 EC)\
                # score 기준치 이하 데이터 제거 (non enzyme 으로 취급)
                print('phase 7: ', end='')
                t = time.time()
                drop_idx = []
                all_ec_4d = pd.read_csv(class_ec_4digit)
                all_ec_4d = list(all_ec_4d['EC_number_4d'])
                for j in range(len(ENC_label_df)):
                    if ENC_label_df['EC Number'].iloc[j] not in all_ec_4d:
                        drop_idx.append(j)
                    elif ENC_label_df['Score'].iloc[j] < enc_score_threshold:
                        drop_idx.append(j)
                ENC_label_df = ENC_label_df.drop(drop_idx)
                # gene_id 하나당 score 가장 높은 EC number 하나만 남김
                ENC_label_df = ENC_label_df.sort_values(by=['gene_id', 'Score'], ascending=[True, False])
                ENC_label_df = ENC_label_df.drop_duplicates(subset=['gene_id'], keep='first')
                ENC_label_df = ENC_label_df[['gene_id', 'EC Number', 'Score']]
                # data, label 을 left join, EC number 가 없다면 non enzyme 취급 (EC.0.0.0.0)
                total = pd.merge(data_df, ENC_label_df, left_on='gene_id', right_on='gene_id', how='left')
                total['EC Number'] = total['EC Number'].replace(np.nan, 'EC 0.0.0.0')
                total['Score'] = total['Score'].replace(np.nan, 0)

                del ENC_label_df, data_df
                print('done, time spent = %.2f(s)' % (time.time() - t))

                print('phase 8: ', end='')
                t = time.time()
                total.to_csv(os.path.join(root, '개수 조정된 원데이터', 'sEEPP_ENC', f"{df['genome_id'].iloc[i]}.csv"), index=False)
                df2 = pd.read_csv(status_table2)
                df2.loc[i, 'sEEPP total rows'] = len(total)
                # df2.loc[i, 'sEEPP num of conditions'] = num_conditions
                df2.loc[i, 'ENC total rows'] = len(total)
                df2.to_csv(status_table2, index=False)
                print('done, time spent = %.2f(s)' % (time.time() - t))
                print(f'before {prev_len} rows, {prev_len-len(total)} removed, after {len(total)} rows')
                del total

        else:
            if df['genome download done'].iloc[i] and df['BRENDA result download done'].iloc[i] \
                    and not os.path.exists(os.path.join(root, '개수 조정된 원데이터', 'ENC', f"{df['genome_id'].iloc[i]}.csv")):
                print(f"<Preprocessing> GenomeName: {df['genome_name'].iloc[i]}, GenomeID: {df['genome_id'].iloc[i]}, target: ENC")
                fpf_ENC = FilePathFinder(df['genome_id'].iloc[i], 'ENC')

                # Merge.XXXXX 파일 읽고 PK->gene_id 변환
                print('phase 1: ', end='')
                t = time.time()
                data_df = pd.read_csv(fpf_ENC.read_path()['data'])
                data_df['PK'] = data_df['PK'].map(lambda x: x[2:-1])
                data_df.rename(columns={'PK': 'gene_id'}, inplace=True)
                print('done, time spent = %.2f(s)' % (time.time() - t))
                prev_len = len(data_df)
                codon_cols = list(data_df.columns[7:71])

                # ATGCN 이외 문자를 갖는 서열 제거
                print('phase 2: ', end='')
                t = time.time()
                drop_idx = []
                for j in range(len(data_df)):
                    seq = data_df.iloc[j]
                    for col in nucleotide_seq_cols:
                        if str(seq[col]) == 'nan':
                            continue
                        elif len(set(seq[col]) - set('ATGCN')) > 0:
                            drop_idx.append(j)
                            break
                data_df = data_df.drop(drop_idx)
                print('done, time spent = %.2f(s), %d row removed' % (time.time() - t, len(drop_idx)))

                # 아미노산 서열 길이 10 이하, 21개 문자 이외 포함하는 서열 제거
                # 단, 마지막에 ] 들어있는 것은 추가 파싱후 사용, 추후 버그 수정시 해당 기능도 삭제
                print('phase 3: ', end='')
                t = time.time()
                drop_idx = []
                amino_acid = list(data_df['amino_acid'])
                for j in range(len(amino_acid)):
                    seq = amino_acid[j]
                    if len(str(seq)) < 10:
                        drop_idx.append(j)
                    elif len(set(seq) - set('ACDEFGHIKLMNPQRSTVWXY')) > 0:
                        if len(set(seq[:-1]) - set('ACDEFGHIKLMNPQRSTVWXY')) == 0:
                            amino_acid[j] = seq[:-1]
                        else:
                            drop_idx.append(j)
                data_df['amino_acid'] = amino_acid
                data_df = data_df.drop(drop_idx)
                print('done, time spent = %.2f(s), %d row removed' % (time.time() - t, len(drop_idx)))

                # 사전에 정의된 길이로 고정
                print('phase 4: ', end='')
                t = time.time()
                for col in ['cds', 'terminator', 'utr3']:
                    data_df[col] = data_df[col].map(
                        lambda x: 'N' * maxlen[col] if type(x) != str
                        else (x[:maxlen[col]] if len(x) > maxlen[col] else x + 'N' * (maxlen[col] - len(x))))
                for col in ['promoter', 'utr5']:
                    data_df[col] = data_df[col].map(
                        lambda x: 'N' * maxlen[col] if type(x) != str
                        else (x[-maxlen[col]:] if len(x) > maxlen[col] else 'N' * (maxlen[col] - len(x)) + x))
                data_df['amino_acid'] = data_df['amino_acid'].map(
                    lambda x: 'N' * maxlen['amino_acid'] if type(x) != str
                    else (x[:maxlen['amino_acid']] if len(x) > maxlen['amino_acid'] else x + 'N' * (
                            maxlen['amino_acid'] - len(x))))
                # amino_acid 중복 제거, [promoter terminator 5utr 3utr cds] 중복 제거
                data_df = data_df.drop_duplicates(subset=['amino_acid'])
                data_df = data_df.drop_duplicates(subset=nucleotide_seq_cols)
                print('done, time spent = %.2f(s)' % (time.time() - t))

                # Total_XXXXX 파일 읽고 gene_id 변환
                print('phase 5: not required')

                # BRENDA annotation result 파일 읽고 Protein name->gene_id 변환
                print('phase 6: ', end='')
                t = time.time()
                ENC_label_df = pd.read_csv(fpf_ENC.read_path()['label'])
                ENC_label_df['Protein name'] = ENC_label_df['Protein name'].map(lambda x: x[2:-1])
                ENC_label_df.rename(columns={'Protein name': 'gene_id'}, inplace=True)
                print('done, time spent = %.2f(s)' % (time.time() - t))

                # 클래스로 정의되지 않는 데이터 제거 (1~3자리 EC, 또는 B 포함하는 EC)\
                # score 기준치 이하 데이터 제거 (non enzyme 으로 취급)
                print('phase 7: ', end='')
                t = time.time()
                drop_idx = []
                all_ec_4d = pd.read_csv(class_ec_4digit)
                all_ec_4d = list(all_ec_4d['EC_number_4d'])
                for j in range(len(ENC_label_df)):
                    if ENC_label_df['EC Number'].iloc[j] not in all_ec_4d:
                        drop_idx.append(j)
                    elif ENC_label_df['Score'].iloc[j] < enc_score_threshold:
                        drop_idx.append(j)
                ENC_label_df = ENC_label_df.drop(drop_idx)
                # gene_id 하나당 score 가장 높은 EC number 하나만 남김
                ENC_label_df = ENC_label_df.sort_values(by=['gene_id', 'Score'], ascending=[True, False])
                ENC_label_df = ENC_label_df.drop_duplicates(subset=['gene_id'], keep='first')
                ENC_label_df = ENC_label_df[['gene_id', 'EC Number', 'Score']]
                # data, label 을 left join, EC number 가 없다면 non enzyme 취급 (EC.0.0.0.0)

                data_df = data_df[['gene_id'] + ['amino_acid'] + nucleotide_seq_cols + codon_cols]
                total = pd.merge(data_df, ENC_label_df, left_on='gene_id', right_on='gene_id', how='left')
                total['EC Number'] = total['EC Number'].replace(np.nan, 'EC 0.0.0.0')
                total['Score'] = total['Score'].replace(np.nan, 0)

                del ENC_label_df, data_df
                print('done, time spent = %.2f(s)' % (time.time() - t))

                print('phase 8: ', end='')
                t = time.time()
                total.to_csv(os.path.join(root, '개수 조정된 원데이터', 'ENC', f"{df['genome_id'].iloc[i]}.csv"), index=False)
                df2 = pd.read_csv(status_table2)
                df2.loc[i, 'ENC total rows'] = len(total)
                df2.to_csv(status_table2, index=False)
                print('done, time spent = %.2f(s)' % (time.time() - t))
                print(f'before {prev_len} rows, {prev_len - len(total)} removed, after {len(total)} rows')
                del total
# '''