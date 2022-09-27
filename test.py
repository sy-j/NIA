from utils import *


import os

import pandas as pd

dir = os.getcwd()

net = r'I:\2022 프로젝트\2022_NIA_AI학습용데이터\초기데이터 1-Cycle\제출용\초기데이터10퍼센트\label\ENC'
# os.chdir(net)

print(os.getcwd())
# df = pd.read_csv(status_table)
# df.to_csv(os.path.join(net, 'test.csv'))


class Mum:
    def __init__(self, val):
        self.glob = val

    def superfunc(self):
        self.glob = -1


class Test(Mum):
    def __init__(self, val):
        super().__init__(val)
        self.val = 999

    def func(self):
        val = 1
        print(val, self.val)
        self.val = val
        print(val, self.val)


t = Test(10)
print(t.glob)
t.superfunc()
print(t.glob)

