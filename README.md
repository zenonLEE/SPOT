# SPOT: StepwisePSO-based Package for Optimization of Targeted Process Simulation

**Optimization of membrane-assisted methanol synthesis using SPOT: A machine learning-based closed-loop framework**  
📄 *Authors: Yuanming Li, Zhenyu Du, Shuai Deng, Xiaonan Wang, Hao Wang, Shuangjun Li*  
🏛 *Affiliations: Korea University, Peking University, Tianjin University, Tsinghua University*  

**SPOT** (StepwisePSO-based Package for Optimization of Targeted process simulation) is a novel machine learning-based closed-loop framework designed for optimizing membrane-assisted methanol synthesis. The framework integrates **predictive modeling, feature importance analysis, and a stepwise particle swarm optimization (stepwisePSO) algorithm** to refine process parameters efficiently.

---

## 🔥 Key Features
- **Membrane-Assisted Methanol Synthesis Optimization**
- **Stepwise Particle Swarm Optimization (StepwisePSO) Algorithm**
- **Machine Learning-enhanced Process Parameter Tuning**
- **Data-Efficient Closed-Loop Optimization Framework**

---

## 📌 System Requirements
**Before running the code, you must install AspenONE v14.0.**
- Make sure AspenONE v14.0 is correctly installed and licensed on your system.
- The simulation requires Aspen properties for process calculations.

---

## 🚀 Installation & Usage

### **1️⃣ Clone the Repository**
```bash
git clone https://github.com/zenonLEE/SPOT.git
cd SPOT
```

### **2️⃣ Set Up a Virtual Environment Using Anaconda (Recommended)**
```bash
conda create --name spot_env python=3.8 -y
conda activate spot_env
```

### **3️⃣ Install Dependencies**
```bash
pip install -r requirements.txt
```

### **4️⃣ Run Optimization**
```bash
python pso_and_simulation.py
```

---

## ⚙️ Configuration Parameters
The `config/config.yaml` file contains key hyperparameters and settings. Modify the following parameters as needed:

- **Variable Parameters**: Adjust parameters such as `per_H2O1`, `per_S1`, `Fp1`, `T0`, and `P0` to optimize performance.
- **Model Settings**: Change `model_name`, `input_columns`, and `output_columns` as required.
- **Data Path**: Update `data_path` to specify the dataset location.
- **Optimization Settings**: Modify iteration counts and optimization strategies.
- **PSO Parameters**: Adjust `num_particles`, `max_iter`, `mutation_probability`, and other hyperparameters for tuning.
- **Random Seed**: Set `seed` for reproducibility.

Refer to `config/config.yaml` for detailed configurations.

---

## 📜 Citation
If you use **SPOT** in your research, please cite:
```bibtex
@article{Li2025SPOT,
  author    = {Yuanming Li, Zhenyu Du, Shuai Deng, Xiaonan Wang, Hao Wang, Shuangjun Li},
  title     = {Optimization of membrane-assisted methanol synthesis using SPOT: A machine learning-based closed-loop framework},
  journal   = {},
  year      = {2025},
  doi       = {}
}
```

---

## 📩 Contact
If you have any questions, please feel free to contact us:
📧 Email: lym7499500@gmail.com  

