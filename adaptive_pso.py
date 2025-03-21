import numpy as np
import pandas as pd
import logging
import os
import matplotlib.pyplot as plt
from scipy.stats import qmc
from multiprocessing import Pool
from scipy.optimize import minimize
from tqdm import tqdm



class StepwisePSO:
    def __init__(
            self,
            model,
            scaler_X,
            scaler_Y,
            input_cols,
            bounds_lower,
            bounds_upper,
            num_particles=30,
            max_iter=100,
            early_stopping_rounds=10,
            early_stopping_tolerance=1e-5,
            w=0.7,    # Inertia weight
            c1=1.5,   # Cognitive factor
            c2=1.5,   # Social factor
            objective_function=None
    ):
        """
        Basic PSO implementation with staged search and local search enhancement.

        Args:
            model: Trained model for predicting y (e.g., eta).
            scaler_X, scaler_Y: Scalers for input and output normalization.
            input_cols: List of input column names.
            bounds_lower, bounds_upper: Lower and upper bounds for each dimension.
            num_particles: Number of particles.
            max_iter: Maximum iterations.
            early_stopping_rounds: Rounds without improvement to trigger early stopping.
            early_stopping_tolerance: Tolerance for "no improvement".
            w, c1, c2: PSO parameters (inertia, cognitive, social).
            objective_function: Optional custom objective function; defaults to -eta if None.
        """
        self.model = model
        self.scaler_X = scaler_X
        self.scaler_Y = scaler_Y
        self.input_cols = input_cols

        self.bounds_lower = np.array(bounds_lower, dtype=float)
        self.bounds_upper = np.array(bounds_upper, dtype=float)

        self.num_particles = int(num_particles)
        self.max_iter = int(max_iter)
        self.early_stopping_rounds = int(early_stopping_rounds)
        self.early_stopping_tolerance = float(early_stopping_tolerance)

        self.w = float(w)
        self.c1 = float(c1)
        self.c2 = float(c2)

        self.objective_function = objective_function

        self.prediction_cache = {}
        self.max_cache_size = 10000

        self.dim = len(self.bounds_lower)

    def _check_cache_size(self):
        """Control the cache size to avoid unlimited growth."""
        if len(self.prediction_cache) > self.max_cache_size:
            num_to_delete = int(0.2 * self.max_cache_size)
            for _ in range(num_to_delete):
                self.prediction_cache.pop(next(iter(self.prediction_cache)))

    def _initialize_particles(self):
        """Initialize particle positions randomly within bounds and set velocities to zero."""
        particles = []
        for _ in range(self.num_particles):
            rand_unit = np.random.rand(self.dim)
            pos = self.bounds_lower + rand_unit * (self.bounds_upper - self.bounds_lower)
            particles.append(pos)
        particles = np.array(particles, dtype=float)
        velocities = np.zeros_like(particles, dtype=float)
        return particles, velocities

    def get_predictions(self, particle):
        """Return the predicted value (eta) for a given particle using caching."""
        key = tuple(particle)
        if key in self.prediction_cache:
            return self.prediction_cache[key]

        df_input = pd.DataFrame([particle], columns=self.input_cols)
        X_scaled = self.scaler_X.transform(df_input)
        y_scaled = self.model.predict(X_scaled)

        if isinstance(y_scaled, (list, np.ndarray)):
            y_pred = self.scaler_Y.inverse_transform(y_scaled.reshape(-1, 1))[0, 0]
        else:
            y_pred = np.inf
            logging.error(f"Model prediction not array-like: {y_scaled}")

        try:
            y_pred = float(y_pred)
        except Exception as e:
            logging.error(f"Conversion error: {e}")
            y_pred = np.inf

        self.prediction_cache[key] = y_pred
        self._check_cache_size()
        return y_pred

    def get_predictions_with_uncertainty(self, particle):
        """Return the predicted value and uncertainty (sigma) for a given particle using caching."""
        key = tuple(particle)
        if key in self.prediction_cache:

            return self.prediction_cache[key]

        df_input = pd.DataFrame([particle], columns=self.input_cols)
        X_scaled = self.scaler_X.transform(df_input)
        try:

            y_scaled, sigma = self.model.predict(X_scaled, return_std=True)
        except TypeError:

            y_scaled = self.model.predict(X_scaled)
            sigma = 0.0

        if isinstance(y_scaled, (list, np.ndarray)):
            y_pred = self.scaler_Y.inverse_transform(y_scaled.reshape(-1, 1))[0, 0]
        else:
            y_pred = np.inf
            logging.error(f"Model prediction not array-like: {y_scaled}")

        try:
            y_pred = float(y_pred)
            sigma = float(sigma)
        except Exception as e:
            logging.error(f"Conversion error: {e}")
            y_pred = np.inf
            sigma = 0.0

        self.prediction_cache[key] = (y_pred, sigma)
        self._check_cache_size()
        return y_pred, sigma


    def evaluate_particle(self, position):
        """Evaluate the fitness (cost) of a particle (to be minimized)."""
        if self.objective_function is not None:
            cost = self.objective_function(position)
        else:
            cost = self._internal_objective_function(position)

        if not isinstance(cost, float):
            logging.error(f"Cost is not float: {cost}, forcing to inf.")
            cost = np.inf
        return cost

    def _internal_objective_function(self, particle):
        """
        Default objective: maximize predicted eta by minimizing - (eta + k * sigma),
        where sigma is the prediction uncertainty and k is a trade-off parameter.
        """
        eta_pred, sigma = self.get_predictions_with_uncertainty(particle)
        if np.isnan(eta_pred):
            return np.inf

        exploration_factor = 0.15
        return - (eta_pred + exploration_factor * sigma)

    def _update_velocity_position(self, particles, velocities, personal_best_pos, global_best_pos):
        """
        Standard PSO velocity and position update.
        After updating, apply reflection: if a particle goes out of bounds,
        set it to the boundary and reverse its velocity for that dimension.
        """
        for i in range(self.num_particles):
            r1 = np.random.rand(self.dim)
            r2 = np.random.rand(self.dim)
            velocities[i] = (self.w * velocities[i] +
                             self.c1 * r1 * (personal_best_pos[i] - particles[i]) +
                             self.c2 * r2 * (global_best_pos - particles[i]))
            particles[i] += velocities[i]

            for d in range(self.dim):
                if particles[i][d] < self.bounds_lower[d]:
                    particles[i][d] = self.bounds_lower[d]
                    velocities[i][d] = -velocities[i][d]
                elif particles[i][d] > self.bounds_upper[d]:
                    particles[i][d] = self.bounds_upper[d]
                    velocities[i][d] = -velocities[i][d]
        return particles, velocities

    def pso_optimize(self, run_folder="pso_run"):
        """
        Full-space PSO search.
        """
        os.makedirs(run_folder, exist_ok=True)
        particles, velocities = self._initialize_particles()
        personal_best_pos = particles.copy()
        personal_best_cost = [self.evaluate_particle(p) for p in particles]

        global_best_index = np.argmin(personal_best_cost)
        global_best_pos = personal_best_pos[global_best_index].copy()
        global_best_cost = personal_best_cost[global_best_index]

        cost_history = []
        eta_history = []
        best_costs = []
        no_improvement_count = 0

        logging.info("Starting full-space PSO optimization...")
        for iteration in tqdm(range(self.max_iter), desc="PSO Iterations", unit="iter"):
            particles, velocities = self._update_velocity_position(
                particles, velocities, personal_best_pos, global_best_pos
            )

            for i in range(self.num_particles):
                cost = self.evaluate_particle(particles[i])
                if cost < personal_best_cost[i]:
                    personal_best_cost[i] = cost
                    personal_best_pos[i] = particles[i].copy()
                    if cost < global_best_cost:
                        global_best_cost = cost
                        global_best_pos = particles[i].copy()

            cost_history.append(global_best_cost)
            best_eta = -global_best_cost  # Since cost = -eta
            eta_history.append(best_eta)
            best_costs.append(global_best_cost)

            if len(best_costs) > self.early_stopping_rounds:
                recent_improvement = best_costs[-self.early_stopping_rounds] - best_costs[-1]
                if recent_improvement < self.early_stopping_tolerance:
                    no_improvement_count += 1
                else:
                    no_improvement_count = 0
                if no_improvement_count >= self.early_stopping_rounds:
                    logging.info(f"Early stopping triggered at iteration {iteration + 1}.")
                    break

            logging.info(f"[Full PSO] Iter {iteration + 1:3d} | Best cost = {global_best_cost:.6f} | Eta = {best_eta:.4f}")

        final_eta = -global_best_cost
        self._visualize_optimization(run_folder, cost_history, eta_history)
        logging.info("Full-space PSO optimization finished.")
        return {
            'best_params': global_best_pos,
            'best_cost': global_best_cost,
            'best_eta': final_eta,
            'cost_history': cost_history,
            'eta_history': eta_history,
            'iterations': iteration + 1
        }

    def pso_optimize_subset(self, search_dims, fixed_values, run_folder="pso_subset", max_iter=None, stage_name="Subspace PSO"):
        """
        PSO optimization in a subspace defined by search_dims, with other dimensions fixed.
        The 'stage_name' parameter is used to output more descriptive logging messages.

        Args:
            search_dims: List/array of indices to be searched (e.g., [4,5] for TP optimization).
            fixed_values: Dictionary {dim_idx: value} for dimensions not in search_dims.
            run_folder (str): Folder to save results.
            max_iter (int): Overrides self.max_iter for this stage, if provided.
            stage_name (str): Descriptive name for this stage (e.g., "Stage 1: TP Optimization").

        Returns:
            Dictionary with the best solution and history information.
        """
        if max_iter is None:
            max_iter = self.max_iter

        os.makedirs(run_folder, exist_ok=True)
        search_dims = list(search_dims)
        sub_dim = len(search_dims)

        sub_bounds_lower = self.bounds_lower[search_dims]
        sub_bounds_upper = self.bounds_upper[search_dims]

        sub_particles = []
        sub_velocities = []
        for _ in range(self.num_particles):
            rand_unit = np.random.rand(sub_dim)
            pos_sub = sub_bounds_lower + rand_unit * (sub_bounds_upper - sub_bounds_lower)
            sub_particles.append(pos_sub)
            sub_velocities.append(np.zeros_like(pos_sub))
        sub_particles = np.array(sub_particles, dtype=float)
        sub_velocities = np.array(sub_velocities, dtype=float)

        def get_full_position(subpos):
            full = np.zeros_like(self.bounds_lower)
            for d, val in fixed_values.items():
                full[d] = val
            for i, dim in enumerate(search_dims):
                full[dim] = subpos[i]
            return full

        personal_best_subpos = sub_particles.copy()
        personal_best_cost = []
        for sp in sub_particles:
            full_pos = get_full_position(sp)
            cost = self.evaluate_particle(full_pos)
            personal_best_cost.append(cost)
        personal_best_cost = np.array(personal_best_cost)

        gbest_idx = np.argmin(personal_best_cost)
        global_best_subpos = personal_best_subpos[gbest_idx].copy()
        global_best_cost = personal_best_cost[gbest_idx]

        cost_history = []
        eta_history = []
        best_costs = []
        no_improvement_count = 0

        logging.info(f"Starting {stage_name} optimization in the subspace...")
        for iteration in tqdm(range(max_iter), desc=f"{stage_name} Iterations", unit="iter"):
            for i in range(self.num_particles):
                r1 = np.random.rand(sub_dim)
                r2 = np.random.rand(sub_dim)
                sub_velocities[i] = (self.w * sub_velocities[i] +
                                     self.c1 * r1 * (personal_best_subpos[i] - sub_particles[i]) +
                                     self.c2 * r2 * (global_best_subpos - sub_particles[i]))
                sub_particles[i] += sub_velocities[i]
                for j in range(sub_dim):
                    if sub_particles[i][j] < sub_bounds_lower[j]:
                        sub_particles[i][j] = sub_bounds_lower[j]
                        sub_velocities[i][j] = -sub_velocities[i][j]
                    elif sub_particles[i][j] > sub_bounds_upper[j]:
                        sub_particles[i][j] = sub_bounds_upper[j]
                        sub_velocities[i][j] = -sub_velocities[i][j]

            for i in range(self.num_particles):
                full_pos = get_full_position(sub_particles[i])
                cost = self.evaluate_particle(full_pos)
                if cost < personal_best_cost[i]:
                    personal_best_cost[i] = cost
                    personal_best_subpos[i] = sub_particles[i].copy()
                    if cost < global_best_cost:
                        global_best_cost = cost
                        global_best_subpos = sub_particles[i].copy()

            cost_history.append(global_best_cost)
            best_eta = -global_best_cost
            eta_history.append(best_eta)
            best_costs.append(global_best_cost)

            if len(best_costs) > self.early_stopping_rounds:
                recent_improvement = best_costs[-self.early_stopping_rounds] - best_costs[-1]
                if recent_improvement < self.early_stopping_tolerance:
                    no_improvement_count += 1
                else:
                    no_improvement_count = 0
                if no_improvement_count >= self.early_stopping_rounds:
                    logging.info(f"[{stage_name}] Early stopping triggered at iteration {iteration + 1}.")
                    break

            logging.info(f"[{stage_name}] Iter {iteration + 1:3d} | Best cost = {global_best_cost:.6f} | Eta = {best_eta:.4f}")

        final_eta = -global_best_cost
        self._visualize_optimization(run_folder, cost_history, eta_history)
        best_subpos = global_best_subpos
        best_fullpos = get_full_position(best_subpos)
        logging.info(f"{stage_name} optimization finished.")
        return {
            'best_subpos': best_subpos,
            'best_params': best_fullpos,
            'best_cost': global_best_cost,
            'best_eta': final_eta,
            'cost_history': cost_history,
            'eta_history': eta_history,
            'iterations': iteration + 1
        }

    def _local_search(self, position, maxiter=50):
        """
        Perform local search near the given position using Nelder-Mead.
        Returns improved position and its cost if successful; otherwise returns the original.
        In addition, the returned position is clipped to the specified bounds.
        """

        def objective(x):
            return self.evaluate_particle(x)

        result = minimize(objective, position, method='Nelder-Mead', options={'maxiter': maxiter})

        # 对结果进行投影，确保每个维度在预设的边界内
        improved_position = np.clip(result.x, self.bounds_lower, self.bounds_upper)
        improved_cost = self.evaluate_particle(improved_position)

        if result.success and improved_cost < self.evaluate_particle(position):
            logging.info("Local search converged and found an improved solution (after clipping).")
            return improved_position, improved_cost
        else:
            logging.warning("Local search did not converge, returning original solution.")
            return position, self.evaluate_particle(position)

    def two_stage_optimize(self,
                           stage1_dims,
                           stage1_fixed,
                           stage2_dims,
                           stage2_fixed,
                           stage1_iter=100,
                           stage2_iter=100,
                           folder_stage1="stage1",
                           folder_stage2="stage2",
                           local_search_maxiter=50):
        """
        Two-stage optimization:
          - Stage 1: Optimize TP parameters (e.g., T0 and P0) using PSO.
          - Stage 2: Optimize the remaining parameters with TP fixed (or updated).
          - Finally, refine Stage 2 result using local search.
        """
        logging.info("=== Starting Two-Stage Optimization: Stage 1 (TP Optimization) ===")
        res1 = self.pso_optimize_subset(
            search_dims=stage1_dims,
            fixed_values=stage1_fixed,
            run_folder=folder_stage1,
            max_iter=stage1_iter,
            stage_name="Stage 1: TP Optimization"
        )
        best_params_stage1 = res1['best_params']
        best_eta_stage1 = res1['best_eta']
        logging.info(f"[Two-Stage][Stage 1] Best eta = {best_eta_stage1:.4f}, Best parameters = {best_params_stage1}")

        # Update fixed values for Stage 2 with Stage 1 results
        for d in stage1_dims:
            stage2_fixed[d] = best_params_stage1[d]

        logging.info("=== Starting Two-Stage Optimization: Stage 2 (Other Parameters Optimization) ===")
        res2 = self.pso_optimize_subset(
            search_dims=stage2_dims,
            fixed_values=stage2_fixed,
            run_folder=folder_stage2,
            max_iter=stage2_iter,
            stage_name="Stage 2: Other Parameters Optimization"
        )
        best_params_stage2 = res2['best_params']
        best_eta_stage2 = res2['best_eta']
        best_cost_stage2 = res2['best_cost']
        logging.info(f"[Two-Stage][Stage 2] Preliminary best eta = {best_eta_stage2:.4f}, Best parameters = {best_params_stage2}")

        logging.info("=== Starting Local Search to refine Stage 2 result ===")
        loc_pos, loc_cost = self._local_search(best_params_stage2, maxiter=local_search_maxiter)
        loc_eta = -loc_cost

        if loc_cost < best_cost_stage2:
            logging.info(f"[Two-Stage] Local search improved solution: cost reduced from {best_cost_stage2:.6f} to {loc_cost:.6f}, eta improved to {loc_eta:.4f}")
            res2['best_params'] = loc_pos
            res2['best_cost'] = loc_cost
            res2['best_eta'] = loc_eta
            best_params_stage2 = loc_pos
            best_eta_stage2 = loc_eta
        else:
            logging.info("Local search did not find a better solution; keeping Stage 2 result.")

        return {
            'stage1_result': res1,
            'stage2_result': res2,
            'final_best_params': best_params_stage2,
            'final_best_eta': best_eta_stage2,
            'final_best_cost': best_cost_stage2
        }

    def _visualize_optimization(self, run_folder, cost_history, eta_history):
        """Visualize the optimization process: plot cost and eta versus iterations."""
        plt.figure(figsize=(12, 4))

        plt.subplot(1, 2, 1)
        plt.plot(cost_history, color='blue')
        plt.xlabel('Iteration')
        plt.ylabel('Best Cost')
        plt.title('Cost over Iterations')

        plt.subplot(1, 2, 2)
        plt.plot(eta_history, color='green')
        plt.xlabel('Iteration')
        plt.ylabel('Best Eta')
        plt.title('Eta over Iterations')

        plt.tight_layout()
        fig_path = os.path.join(run_folder, 'pso_basic_progress.png')
        plt.savefig(fig_path, dpi=120)
        plt.close()
        logging.info(f"Optimization progress figure saved to {fig_path}")


