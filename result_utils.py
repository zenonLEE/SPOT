import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
import warnings
import os
import logging
import json

warnings.filterwarnings("ignore")



def load_config(config_path='./config/config.yaml'):
    """
    Load the configuration from a YAML file with UTF-8 encoding.
    Ensure that all numerical parameters are correctly typed.
    """
    with open(config_path, 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)

    # Provide default values if keys are missing
    fixed_params = config.get('fixed_parameters', {})
    variable_params = config.get('variable_parameters', {})

    # Ensure fixed_parameters and variable_parameters are dictionaries
    if not isinstance(fixed_params, dict):
        fixed_params = {}
    if not isinstance(variable_params, dict):
        variable_params = {}

    # Convert bounds and parameters to appropriate types
    for param, bounds in variable_params.items():
        try:
            bounds['lower'] = float(bounds.get('lower', 0))
            bounds['upper'] = float(bounds.get('upper', 1))
        except (ValueError, TypeError) as e:
            logging.error(f"Error in variable parameter bounds for {param}: {e}")

    for param, value in fixed_params.items():
        try:
            fixed_params[param] = float(value)
        except ValueError as e:
            logging.error(f"Error converting fixed parameter {param} to float: {e}")

    config['fixed_parameters'] = fixed_params
    config['variable_parameters'] = variable_params

    return config

def update_membrane_json(params, json_path):
    """
    Update the membrane data section in the JSON file based on the given params.
    The params dictionary should contain:
      - per_H2O: New value for per_H2O
      - per_S: New value for per_S
      - fp1: New value for the first number in Fp
      - fp2: New value for the second number in Fp
    """
    # Read the original JSON data
    with open(json_path, "r") as f:
        data = json.load(f)

    # Update per_H2O and per_S (ensuring the nested list format is maintained)
    data["membrane"]["per_H2O"] = [[params["per_H2O1"]]]
    data["membrane"]["per_S"] = [[params["per_S1"]]]

    # Update the Fp section, which is a tab-separated string
    fp_values = data["membrane"]["Fp"].split("\t")
    # Modify the first and second values
    fp_values[0] = str(params["Fp1"])
    fp_values[1] = str(params["Fp2"])
    # Reassemble the modified list into a string
    data["membrane"]["Fp"] = "\t".join(fp_values)

    # Write the updated data back to the file
    with open(json_path, "w") as f:
        json.dump(data, f, indent=4)

def over_sample_high_eta(df_data, eta_threshold=0.9, multiple=5):

    df_high = df_data[df_data['eta_mod_heat'] > eta_threshold]
    if len(df_high) == 0:
        return df_data
    df_oversampled_high = pd.concat([df_high]*multiple, ignore_index=True)
    df_augmented = pd.concat([df_data, df_oversampled_high], ignore_index=True)
    return df_augmented

def plot_feature_importance(model, input_cols, run_folder):
    """
    Plot and save feature importance for tree-based models.
    """
    try:
        if hasattr(model, 'feature_importances_'):
            importances = model.feature_importances_
            indices = np.argsort(importances)[::-1]
            plt.figure(figsize=(10, 6))
            plt.title("Feature Importances")
            plt.bar(range(len(input_cols)), importances[indices], align='center')
            plt.xticks(range(len(input_cols)), [input_cols[i] for i in indices], rotation=90)
            plt.xlabel('Feature')
            plt.ylabel('Importance')
            plt.tight_layout()
            plt.savefig(os.path.join(run_folder, "feature_importances.png"))
            plt.close()
    except Exception as e:
        print(f"Error plotting feature importances: {e}")


def create_results_directory(base_dir='results'):
    """
    Create a base results directory. If it exists, create a timestamped subdirectory.
    """
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    # Create a timestamped folder to avoid overwriting
    timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
    run_folder = os.path.join(base_dir, f"run_{timestamp}")
    os.makedirs(run_folder)
    return run_folder


def plot_model_evaluation(y_test_original, y_pred_original, model_name, save_dir):
    """
    Plot and save model evaluation visualizations including scatter plot and residuals distribution.

    Args:
        y_test_original: Original scale test values
        y_pred_original: Original scale predicted values
        model_name: Name of the model being evaluated
        save_dir: Directory to save the plots
    """
    # True vs Predicted scatter plot
    plt.figure(figsize=(6, 6))
    plt.scatter(y_test_original, y_pred_original, alpha=0.5)
    plt.xlabel('True Values')
    plt.ylabel('Predicted Values')
    plt.title(f'{model_name.capitalize()} True vs Predicted')
    plt.plot(
        [y_test_original.min(), y_test_original.max()],
        [y_test_original.min(), y_test_original.max()],
        'r--'
    )
    plt.grid(True)
    scatter_plot_path = os.path.join(save_dir, f"{model_name}_true_vs_predicted.png")
    plt.savefig(scatter_plot_path)
    plt.close()

    # Residuals distribution histogram
    residuals = y_test_original - y_pred_original
    plt.figure(figsize=(6, 4))
    plt.hist(residuals, bins=50, color='blue', edgecolor='k', alpha=0.7)
    plt.xlabel('Residuals')
    plt.ylabel('Frequency')
    plt.title(f'{model_name.capitalize()} Residuals Distribution')
    plt.grid(True)
    residuals_plot_path = os.path.join(save_dir, f"{model_name}_residuals_distribution.png")
    plt.savefig(residuals_plot_path)
    plt.close()

    return scatter_plot_path, residuals_plot_path


def plot_optimization_history(best_result, run_folder):
    """
    Plot and save optimization history visualizations including cost, eta, and diversity history.

    Args:
        best_result: Dictionary containing optimization history data
        run_folder: Directory to save the plots

    Returns:
        tuple: Paths to the saved plots
    """
    # Cost history plot
    plt.figure(figsize=(6, 4))
    plt.plot(best_result['cost_history'], marker='o')
    plt.title('Cost History')
    plt.xlabel('Iteration')
    plt.ylabel('Minimum Cost')
    plt.grid(True)
    cost_history_path = os.path.join(run_folder, "best_run_cost_history.png")
    plt.savefig(cost_history_path)
    plt.close()

    # Eta value history plot
    plt.figure(figsize=(6, 4))
    plt.plot(best_result['eta_history'], marker='o', color='green')
    plt.title('Eta Value History')
    plt.xlabel('Iteration')
    plt.ylabel('Eta')
    plt.grid(True)
    eta_history_path = os.path.join(run_folder, "best_run_eta_history.png")
    plt.savefig(eta_history_path)
    plt.close()

    # Particle diversity history plot
    plt.figure(figsize=(6, 4))
    plt.plot(best_result['diversity_history'], marker='o', color='red')
    plt.title('Particle Diversity History')
    plt.xlabel('Iteration')
    plt.ylabel('Diversity')
    plt.grid(True)
    diversity_history_path = os.path.join(run_folder, "best_run_diversity_history.png")
    plt.savefig(diversity_history_path)
    plt.close()

    return cost_history_path, eta_history_path, diversity_history_path

