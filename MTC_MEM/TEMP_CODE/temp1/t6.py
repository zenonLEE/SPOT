import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap


def load_df(file, sheet_name=None, lang='en'):
    """
    load Dataframe from file with sheet of sheet_name
    if not specified, load all
    """
    out_df_data = {}
    sheet_name = pd.ExcelFile(file).sheet_names if sheet_name is None else sheet_name
    for sheet in sheet_name:
        temp = pd.read_excel(file, sheet_name=sheet)
        skip_row = temp.index[temp.iloc[:, 0] == 'Date'].tolist()[0]
        if lang == 'en':
            df_data = pd.read_excel(file, skiprows=skip_row + 1, index_col='Date', sheet_name=sheet, na_values=['Wind'])
        elif lang == 'cn':
            skip_row = [0, 1, 3] if skip_row == 2 else [0, 2]
            df_data = pd.read_excel(file, skiprows=skip_row, index_col='日期', sheet_name=sheet, na_values=['Wind'])

        df_data.index = pd.to_datetime(df_data.index).strftime('%Y-%m-%d')
        out_df_data[sheet] = df_data
    return out_df_data


def rank(df, size):
    """
    calculate the quantile
    :param df:
    :param size:
    :return:
    """
    df_rank = pd.Series(index=df.index)
    for i in np.arange(size - 1, len(df)):
        temp = df.iloc[i - size + 1:i + 1].copy()
        df_rank.iloc[i] = temp.rank(pct=True, na_option='keep').iloc[-1]
    return df_rank


def annul_interest(df, size):
    """
    calculate the annual interest within window of size
    :param df: close of a share (pd.Series)
    :param size: window size
    """

    vol_1 = df.pct_change(periods=1)
    vol_size = vol_1.rolling(window=size).std()
    interest = (vol_size * 250 ** 0.5).dropna()
    length = 250 * 10 if len(interest) >= 2500 else len(interest) - 1
    quant = interest.rolling(length).rank(pct=True, numeric_only=True)  # rank(interest, 250 * 10)  #
    return pd.concat([interest, quant], axis=1)


def daily_interest(df, size):
    """
    calculate the daily interest within window of size
    :param df: close of a share (pd.Series)
    :param size: window size
    """

    interest = df.pct_change(periods=size).dropna()
    length = 250 * 10 if len(interest) >= 2500 else len(interest) - 1
    quant = interest.rolling(length).rank(pct=True, numeric_only=True)  # rank(interest, 250 * 10) #
    return pd.concat([interest, quant], axis=1)


def cal_index(input_data):
    """
    calculate the performance for the input share
    :param input_data: dict
    :return:
    """
    df_sheet_name, df_data = list(input_data.keys())[0], list(input_data.values())[0]
    out_data = {df_sheet_name: pd.DataFrame(index=df_data.index)}

    if df_sheet_name == 'close':
        for share in df_data.columns:
            # sheet_names = [share + i for i in ['_AI_60', '_AI_60_Q10', '_AI_250', '_AI_250_Q10',
            #                                    '_DI_5', '_DI_5_Q10', '_DI_20', '_DI_20_Q10']]
            sheet_names = [share + i for i in ['_AI_60', '_60日波动率', '_AI_250', '_250日波动率',
                                               '_DI_5', '_5日涨跌幅', '_DI_20', '_20日涨跌幅']]
            share_df = pd.concat([annul_interest(df_data[share], 60),
                                  annul_interest(df_data[share], 250),
                                  daily_interest(df_data[share], 5),
                                  daily_interest(df_data[share], 20)], axis=1)
            share_df.columns = sheet_names
            out_data[df_sheet_name] = pd.concat([out_data[df_sheet_name], share_df], axis=1)
    # elif df_sheet_name == 'chg_week' or df_sheet_name == '周涨跌幅':
    #     share_df = df_data.apply(lambda x: x.dropna().rank(pct=True, na_option='keep'), axis=0)
    #     share_df.columns = [i + '_' + '周涨跌幅' for i in df_data.columns]  # + '_分位数'
    #     out_data[df_sheet_name] = share_df
    else:

        length = 500 if df_sheet_name == 'chg_week' else 250
        share_df = df_data.apply(lambda x: x.dropna().rolling(length).rank(pct=True) if len(x.dropna()) >= length
        else x.dropna().rank(pct=True))  # rank(x, 250)
        out_data[df_sheet_name] = share_df
        index_names = {'chg': '周涨跌幅', 'amt': '成交量', 'turn': '换手率'}
        index_name = index_names[df_sheet_name.split('_')[0]]
        out_data[df_sheet_name].columns = [i + '_' + index_name for i in df_data.columns]  # + '_分位数'
    return out_data


def share_index(shares, date, all_data):
    """
    find the index of share at date
    :param all_data: dict
    :param shares: share ID (list)
    :param date: date in 'y-m-d'
    :return:
    """
    date = pd.to_datetime(date).strftime('%Y-%m-%d')
    share_df = pd.DataFrame(index=shares)
    for data_name, data in all_data.items():
        try:
            select_data = data.loc[date]
            for index in select_data.index:
                keys = index.split("_", 1)
                for share in shares:
                    if share in keys:
                        share_df.loc[share, keys[1]] = select_data[keys[0] + '_' + keys[1]]
        except KeyError:
            print(f'data are not avail at {date} in {data_name}!')
            continue
    return share_df


def draw_res(share_res, only_quant=True):
    """
    visualization of share performance
    """
    if only_quant:
        # only consider quantile
        drop_index = ['AI_60', 'AI_250', 'DI_5', 'DI_20']
        share_res = share_res.drop(drop_index, axis=1)
    share_names = share_res.index.tolist()
    share_indexes = share_res.columns.tolist()
    share_values = np.round(share_res.values, 2)

    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.figure(figsize=(10, 6))  # size of figure
    blues = ['#54BAC2', '#BEFCFF']
    grays = ['#FFFFFF'] * 6  # ['#808080', '#A0A0A0', '#C0C0C0', '#E0E0E0', '#C0C0C0', '#A0A0A0'] #FDFBFB #EDE9D0
    reds = ['#F8B1A2', '#F25945']
    cmap = ListedColormap(blues + grays + reds)

    ax = sns.heatmap(share_values, cmap=cmap, annot=True, cbar=False, fmt='.0%',
                     yticklabels=share_names, vmin=0, vmax=1)  # 'coolwarm'
    ax.set_xticklabels(share_indexes, rotation=0)
    ax.xaxis.set_ticks_position('top')
    cbar = plt.colorbar(ax.collections[0], ticks=[0, 0.2, 0.8, 1])
    cbar.ax.set_yticklabels(['0%', '20%', '80%', '100%'])
    plt.show()


def save_res(save_data):
    """
    save data to excel
    """
    if isinstance(save_data, dict):
        with pd.ExcelWriter(out_file_data, engine='openpyxl') as writer:
            for sheet_df, data_df in save_data.items():
                data_df.to_excel(writer, sheet_name=sheet_df)
    else:
        save_data.to_excel(out_file_index, engine='openpyxl')


if __name__ == "__main__":

    # give your file path
    in_file = r"D:\document\02Work\temp\多资产收盘价data.xlsx"  # input file path
    in_path, in_file_name = os.path.split(in_file)
    out_file_data = os.path.join(in_path, 'data_' + in_file_name)  # output file path
    out_file_index = os.path.join(in_path, 'index_' + in_file_name)  # output file path
    language = 'cn'  # or 'en'

    # load data for specific sheet in Excel
    target_sheets = pd.ExcelFile(in_file).sheet_names  # ['chg_week']  #
    load_data = load_df(in_file, target_sheets, lang=language)

    # calculate performance of all share at all dates
    out_dfs_data = {}
    for target_sheet in target_sheets:
        temp_dict = {target_sheet: load_data[target_sheet]}
        out_dfs_data.update(cal_index(temp_dict))
    # print(out_dfs_data)
    # save_res(out_dfs_data)  # save data to excel

    # extract performance for target share
    share_IDs = load_data[target_sheets[0]].columns  # the shares you need to extract
    res = share_index(share_IDs, '2023/10/20', out_dfs_data)
    # save_res(res)
    draw_res(res)
