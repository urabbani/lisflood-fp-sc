"""
Visualization utilities for calibration results.
Generates hydrograph comparisons, inundation maps, and calibration progress plots.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Any

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10


def plot_calibration_progress(history: List[Dict[str, Any]], 
                           output_path: str):
    """
    Plot calibration progress over iterations.
    
    Args:
        history: List of calibration iteration results
        output_path: Path to save plot
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    iterations = [r["iteration"] for r in history]
    nse_scores = [r["metrics"]["nse"] for r in history]
    kge_scores = [r["metrics"]["kge"] for r in history]
    iou_scores = [r["metrics"]["iou"] for r in history]
    composite_scores = [r["score"] for r in history]
    
    # NSE progress
    axes[0, 0].plot(iterations, nse_scores, 'b-', linewidth=1.5)
    axes[0, 0].set_xlabel('Iteration')
    axes[0, 0].set_ylabel('NSE')
    axes[0, 0].set_title('Nash-Sutcliffe Efficiency')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].axhline(y=0.85, color='r', linestyle='--', alpha=0.5, label='Target (0.85)')
    axes[0, 0].legend()
    
    # KGE progress
    axes[0, 1].plot(iterations, kge_scores, 'g-', linewidth=1.5)
    axes[0, 1].set_xlabel('Iteration')
    axes[0, 1].set_ylabel('KGE')
    axes[0, 1].set_title('Kling-Gupta Efficiency')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].axhline(y=0.80, color='r', linestyle='--', alpha=0.5, label='Target (0.80)')
    axes[0, 1].legend()
    
    # IoU progress
    axes[1, 0].plot(iterations, iou_scores, 'm-', linewidth=1.5)
    axes[1, 0].set_xlabel('Iteration')
    axes[1, 0].set_ylabel('IoU')
    axes[1, 0].set_title('Intersection over Union')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].axhline(y=0.70, color='r', linestyle='--', alpha=0.5, label='Target (0.70)')
    axes[1, 0].legend()
    
    # Composite score progress
    axes[1, 1].plot(iterations, composite_scores, 'k-', linewidth=1.5)
    axes[1, 1].set_xlabel('Iteration')
    axes[1, 1].set_ylabel('Composite Score')
    axes[1, 1].set_title('Weighted Composite Score')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   Saved calibration progress: {output_path}")


def plot_pareto_front(pareto_front: List[Dict[str, Any]],
                     output_path: str):
    """
    Plot Pareto front of non-dominated solutions.
    
    Args:
        pareto_front: List of non-dominated solutions
        output_path: Path to save plot
    """
    if len(pareto_front) == 0:
        print("   Warning: No Pareto front solutions to plot")
        return
    
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Extract metrics
    nse = [r["metrics"]["nse"] for r in pareto_front]
    kge = [r["metrics"]["kge"] for r in pareto_front]
    iou = [r["metrics"]["iou"] for r in pareto_front]
    scores = [r["score"] for r in pareto_front]
    
    # Color by composite score
    scatter = ax.scatter(nse, kge, iou, c=scores, cmap='viridis', s=100, alpha=0.7)
    
    # Add labels for best solution
    best_idx = np.argmax(scores)
    ax.scatter(nse[best_idx], kge[best_idx], iou[best_idx], 
              c='red', s=200, marker='*', label='Best Solution')
    
    ax.set_xlabel('NSE', fontsize=12)
    ax.set_ylabel('KGE', fontsize=12)
    ax.set_zlabel('IoU', fontsize=12)
    ax.set_title('Pareto Front - Multi-Objective Calibration', fontsize=14, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label('Composite Score', rotation=270, labelpad=20)
    
    # Set axis limits
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   Saved Pareto front: {output_path}")


def plot_hydrograph_comparison(gauge_data: Dict[str, np.ndarray],
                             simulated_data: Dict[str, np.ndarray],
                             output_path: str):
    """
    Plot observed vs simulated discharge hydrograph.
    
    Args:
        gauge_data: Dict with "timestamps" and "discharge"
        simulated_data: Dict with "timestamps" and "values"
        output_path: Path to save plot
    """
    fig, ax = plt.subplots(figsize=(14, 6))
    
    obs_time = gauge_data["timestamps"]
    obs_q = gauge_data["discharge"]
    sim_time = simulated_data["timestamps"]
    sim_q = simulated_data["values"]
    
    # Align lengths (crop to shorter)
    min_len = min(len(obs_q), len(sim_q))
    obs_q = obs_q[:min_len]
    sim_q = sim_q[:min_len]
    obs_time = obs_time[:min_len]
    sim_time = sim_time[:min_len]
    
    # Plot
    ax.plot(obs_time, obs_q, 'b-', linewidth=2, label='Observed', alpha=0.7)
    ax.plot(sim_time, sim_q, 'r--', linewidth=2, label='Simulated', alpha=0.7)
    
    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.set_ylabel('Discharge (m³/s)', fontsize=12)
    ax.set_title('Hydrograph Comparison', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Calculate metrics
    from calibration.metrics import CalibrationMetrics
    metrics = CalibrationMetrics()
    nse = metrics.nash_sutcliffe(obs_q, sim_q)
    kge = metrics.kling_gupta(obs_q, sim_q)
    
    # Add text box with metrics
    textstr = f'NSE = {nse:.3f}\nKGE = {kge:.3f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   Saved hydrograph comparison: {output_path}")


def plot_inundation_comparison(observed_mask: np.ndarray,
                             simulated_mask: np.ndarray,
                             output_path: str):
    """
    Plot observed vs simulated inundation extent.
    
    Args:
        observed_mask: Observed flood mask (2D array)
        simulated_mask: Simulated flood mask (2D array)
        output_path: Path to save plot
    """
    # Ensure same shape
    min_rows = min(observed_mask.shape[0], simulated_mask.shape[0])
    min_cols = min(observed_mask.shape[1], simulated_mask.shape[1])
    observed_mask = observed_mask[:min_rows, :min_cols]
    simulated_mask = simulated_mask[:min_rows, :min_cols]
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Observed
    im1 = axes[0].imshow(observed_mask, cmap='Blues', vmin=0, vmax=1)
    axes[0].set_title('Observed Inundation', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('X (pixels)')
    axes[0].set_ylabel('Y (pixels)')
    plt.colorbar(im1, ax=axes[0], label='Flooded')
    
    # Simulated
    im2 = axes[1].imshow(simulated_mask, cmap='Blues', vmin=0, vmax=1)
    axes[1].set_title('Simulated Inundation', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('X (pixels)')
    axes[1].set_ylabel('Y (pixels)')
    plt.colorbar(im2, ax=axes[1], label='Flooded')
    
    # Difference
    diff_mask = simulated_mask.astype(int) - observed_mask.astype(int)
    # -1: False positive (sim=1, obs=0)
    # 0: Correct (both 0 or both 1)
    # 1: False negative (sim=0, obs=1)
    
    diff_cmap = plt.cm.RdYlBu
    im3 = axes[2].imshow(diff_mask, cmap=diff_cmap, vmin=-1, vmax=1)
    axes[2].set_title('Difference (Sim - Obs)', fontsize=12, fontweight='bold')
    axes[2].set_xlabel('X (pixels)')
    axes[2].set_ylabel('Y (pixels)')
    cbar = plt.colorbar(im3, ax=axes[2])
    cbar.set_label('Difference')
    cbar.set_ticks([-1, 0, 1])
    cbar.set_ticklabels(['False Positive', 'Correct', 'False Negative'])
    
    # Calculate IoU
    from calibration.metrics import CalibrationMetrics
    metrics = CalibrationMetrics()
    iou = metrics.iou(observed_mask, simulated_mask)
    f1 = metrics.f1_score(observed_mask, simulated_mask)
    
    # Add text with metrics
    textstr = f'IoU = {iou:.3f}\nF1 = {f1:.3f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axes[2].text(0.05, 0.95, textstr, transform=axes[2].transAxes,
                fontsize=11, verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   Saved inundation comparison: {output_path}")


def plot_parameter_evolution(history: List[Dict[str, Any]],
                           param_names: List[str],
                           output_path: str):
    """
    Plot parameter evolution over iterations.
    
    Args:
        history: List of calibration results
        param_names: List of parameter names to plot
        output_path: Path to save plot
    """
    n_params = len(param_names)
    fig, axes = plt.subplots(n_params, 1, figsize=(14, 3 * n_params))
    
    if n_params == 1:
        axes = [axes]
    
    iterations = [r["iteration"] for r in history]
    
    for i, param_name in enumerate(param_names):
        values = [r["params"].get(param_name, np.nan) for r in history]
        axes[i].plot(iterations, values, 'b-', linewidth=1.5)
        axes[i].set_xlabel('Iteration')
        axes[i].set_ylabel(param_name)
        axes[i].set_title(f'{param_name} Evolution', fontsize=11, fontweight='bold')
        axes[i].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   Saved parameter evolution: {output_path}")
