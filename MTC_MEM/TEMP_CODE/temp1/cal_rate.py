import pandas as pd
import os

pd_save = pd.DataFrame(columns=['收益率', '收益率_去除万科'])

# your path to store data
path = r"C:\Users\13323\Documents\cal_test"
save_path = r"C:\Users\13323\Documents\cal_test\res1.xlsx"
data_name = "资金"
files_path = os.listdir(path)


def cal_rate(rate_path):
    own_value = pd.read_excel(rate_path, index_col=0, header=0)
    infoall = pd.DataFrame(own_value.loc['非金融企业（公司）债券': '剩余期限在1年以上的政府债券、准政府债券'])
    after_filter = list(filter(lambda x: str(x) == '非金融企业（公司）债券'
                                         or str(x)[0:8] == '149358SZ'  # 21万科02
                                         or str(x)[0:8] == '149297SZ'  # 20万科08
                               , infoall.index))
    info = infoall.loc[after_filter]
    rate = info.loc['非金融企业（公司）债券', '综合收益率']  # 第一个要保留的值

    diff = info.iloc[0] - info.iloc[1] - info.iloc[2]
    rate2 = diff['综合收益'] / diff['平均资金占用']  # 第2个要保留的值
    return [rate, rate2]


# walk the data
# i = 0
for file in files_path:
    file_path = os.path.join(path, file)

    try:
        for temp in os.listdir(file_path):
            if data_name in temp:
                data_path = os.path.join(file_path, temp)
                res = cal_rate(data_path)
                pd_save.loc[file] = res
    except NotADirectoryError:
        pass

pd_save.to_excel(save_path)
