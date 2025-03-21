import pandas as pd

def load_and_preprocess_data(data_path, qmh1_option, th1, input_cols, target_cols):
    try:
        df = pd.read_excel(data_path, engine='openpyxl')
    except FileNotFoundError:
        print(f"Data file not found at path: {data_path}")
        return None
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

    # 数据筛选
    df['qmh1'] = df['qmh1'].astype(int)
    df_filtered = df[df['qmh1'] == qmh1_option]
    df_filtered = df_filtered[df_filtered['eta'] >= th1]
    df_filtered = df_filtered.dropna(subset=input_cols + target_cols)

    X = df_filtered[input_cols].values
    Y = df_filtered[target_cols].values.ravel()
    return X, Y, df_filtered
