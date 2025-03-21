import time
import pandas as pd
from read import ReadData
from mem_simulator import Simulation
import json
import os

# Constants
MODEL = 'BU'
RED = '\033[91m'
ENDC = '\033[0m'


def format_pressure(p):
    """Format pressure to 2 decimal places"""
    return round(p, 2)

def update_membrane_json(params, json_path):
    """
    根据 params 参数更新 JSON 文件中膜数据部分
    params 字典中应包含：
      - per_H2O: 新的 per_H2O 数值
      - per_S: 新的 per_S 数值
      - fp1: Fp 中第一个数字的新值
      - fp2: Fp 中第二个数字的新值
    """
    # 读取原始的 JSON 数据
    with open(json_path, "r") as f:
        data = json.load(f)

    # 更新 per_H2O 和 per_S（注意保持嵌套列表的格式）
    data["membrane"]["per_H2O"] = [[params["per_H2O1"]]]
    data["membrane"]["per_S"] = [[params["per_S1"]]]

    # 更新 Fp 部分，Fp 是一个以制表符分隔的字符串
    fp_values = data["membrane"]["Fp"].split("\t")
    # 修改第一个和第二个数字
    fp_values[0] = str(params["Fp1"])
    fp_values[1] = str(params["Fp2"])
    # 将修改后的列表重新拼接成字符串
    data["membrane"]["Fp"] = "\t".join(fp_values)

    # 将更新后的数据写回文件
    with open(json_path, "w") as f:
        json.dump(data, f, indent=4)


def setup_simulation_parameters():
    """Initialize simulation parameters"""
    return {
        'T': 513.15,
        'P': 50.0,  # Format P to 2 decimal places
        'Fp1': 0.5,
        'Fp2': 0.00000000,
        'per_H2O1': 0.000002,
        'per_S1': 200,
        'model': MODEL
    }


def initialize_data(params):
    """Initialize data structures for simulation"""
    json_path = os.path.join(os.path.dirname(__file__), "in_mem.json")
    update_membrane_json(params, json_path)
    # 创建 ReadData 实例，传入指定模型
    in_data = ReadData(kn_model=params['model'])

    # 设置温度和压力条件
    in_data.feed_para['condition']['T'][0] = params['T']
    in_data.feed_para['condition']['P'][0] = format_pressure(params['P'])

    # 设置流量参数
    in_data.mem_data.Fp[0] = params['Fp1']
    in_data.mem_data.Fp[1] = params['Fp2']

    # 直接赋值更新渗透性参数（使用 JSON 文件中的键名）
    in_data.mem_para['per_H2O1'] = [[params['per_H2O1']]]
    in_data.mem_para['per_S1'] = [[params['per_S1']]]
    # 验证赋值是否成功
    # print("After assignment:")
    # print(in_data.mem_data[['per_H2O1', 'per_S1']])

    return in_data


def run_simulation(feed_data, reactor_data, mem_data, chem_data):
    """Run the simulation with nested loops and time tracking"""
    n = 0
    overall_start = time.time()

    print(f"Starting simulation with:")
    print(f"{RED}feed_data shape: {feed_data.shape[0]}{ENDC}")
    print(f"{RED}reactor_data shape: {reactor_data.shape[0]}{ENDC}")
    print(f"{RED}mem_data shape: {mem_data.shape[0]}{ENDC}")

    for i in range(feed_data.shape[0]):
        loop_i_start = time.time()

        for j in range(reactor_data.shape[0]):
            loop_j_start = time.time()

            for k in range(mem_data.shape[0]):
                loop_k_start = time.time()

                # Initialize and run simulation
                sim = Simulation(
                    reactor_data.iloc[j],
                    chem_data,
                    feed_data.iloc[i],
                    mem_data.iloc[k],
                    eos=1,
                    drop=1
                )

                sim.sim(
                    save_profile=0,
                    loop='indirect',
                    rtol=0.001,
                    r_target=None
                )

                n += 1

                # Log iteration time
                loop_k_end = time.time()
                print(f"{RED}Iteration (i={i}, j={j}, k={k}) took {loop_k_end - loop_k_start:.2f} seconds.{ENDC}")

            loop_j_end = time.time()
            print(f"{RED}Inner loop j={j} completed in {loop_j_end - loop_j_start:.2f} seconds.{ENDC}")

        loop_i_end = time.time()
        print(f"{RED}Outer loop i={i} completed in {loop_i_end - loop_i_start:.2f} seconds.{ENDC}")

    overall_end = time.time()
    print(f"{RED}Total time for all iterations: {overall_end - overall_start:.2f} seconds.{ENDC}")
    return n


def main():
    # Set up pandas display options
    pd.set_option('display.max_columns', None)

    # Initialize parameters
    params = setup_simulation_parameters()

    # Initialize data structures
    in_data = initialize_data(params)

    # Extract necessary data for simulation
    reactor_data = in_data.reactor_data
    feed_data = in_data.feed_data
    chem_data = in_data.chem
    mem_data = in_data.mem_data

    # Run simulation
    total_iterations = run_simulation(feed_data, reactor_data, mem_data, chem_data)
    print(f"{RED}Completed {total_iterations} total iterations{ENDC}")


if __name__ == "__main__":
    main()