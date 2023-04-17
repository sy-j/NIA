import pandas as pd

from utils import *


p = r'D:\NIA\학습용 데이터\s3_업로드\2.라벨링데이터\ENC'
ec_list = []
for path, dir, files in os.walk(p):
    for file in files:
        with open(os.path.join(path, file), 'r') as js:
            js = json.load(js)
        label = js['label']
        for i in range(len(label)):
            ec_num = f"EC {label[i]['ec_class']}.{label[i]['ec_subclass']}.{label[i]['ec_subsubclass']}.{label[i]['ec_serial']}"
            ec_list.append(ec_num)
        # print(f'{file}, EC number ({len(set(ec_list))}개)')
ec_list = pd.DataFrame(ec_list)
ec_list.value_counts().to_csv(r'D:\NIA\학습용 데이터\EC number count.csv')

#
#
# fpf = FilePathFinder(df)
# p = fpf.save_path()
#
# def fix_all_data():
#     print('Updating all data...')
#     t = time.time()
#     for i in range(len(df)):
#         if df['sEEPP error'].iloc[i] != '0' and df['sEEPP error'].iloc[i] != 0 and df['sEEPP'].iloc[i] == 1:
#             sEEPP_fix(i)
#         if df['ENC error'].iloc[i] != '0' and df['ENC error'].iloc[i] != 0:
#             ENC_fix(i)
#     df.to_csv(status_table3, index=False)
#
#     print("all process done, total time spent: %2f(s)" % (time.time() - t))