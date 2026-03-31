"""
Self-Calibrating LISFLOOD-FP
Main entry point for automatic flood model calibration.
"""

import argparse
import yaml
import os
import json
import logging
from pathlib import Path

from models.lisflood.adapter import LISFLOODAdapter
from models.lisflood82.adapter import LISFLOOD82Adapter
from calibration.loop import AutoCalibrationLoop
from data.preprocessing import load_observations
from satellite.preprocessor import SatelliteInundation

logger = logging.getLogger(__name__)


def load_config(config_path: str) -> dict:
    """Load YAML configuration file"""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def setup_directories(config: dict):
    """Create necessary directories"""
    directories = [
        "data/outputs",
        "data/outputs/lisflood",
        "data/observations/gauges",
        "data/observations/satellite",
        "visualization",
        "logs"
    ]
    
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
    
    print("✅ Directory structure created")


def main():
    parser = argparse.ArgumentParser(
        description="Self-Calibrating LISFLOOD-FP",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--config",
        type=str,
        default="config/project.yaml",
        help="Path to configuration file"
    )
    
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=None,
        help="Override max iterations from config"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed progress"
    )
    
    parser.add_argument(
        "--save-plots",
        action="store_true",
        help="Generate calibration visualization plots"
    )
    
    args = parser.parse_args()
    
    # Load configuration
    print(f"\n🔧 Loading configuration from: {args.config}")
    config = load_config(args.config)

    # Configure logging
    log_config = config.get("logging", {})
    log_level = getattr(logging, log_config.get("level", "INFO").upper(), logging.INFO)
    log_file = log_config.get("file")
    log_console = log_config.get("console", True)

    handlers = []
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    if log_console:
        handlers.append(logging.StreamHandler())

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers if handlers else None
    )

    # Setup directories
    setup_directories(config)
    
    # Initialize model adapter
    print(f"\n🏗️  Initializing {config['model']['type']} adapter...")
    model_config = config["model"]
    
    if config["model"]["type"] == "lisflood-fp":
        model = LISFLOODAdapter(model_config)
    elif config["model"]["type"] == "lisflood82":
        model = LISFLOOD82Adapter(model_config)
    else:
        raise ValueError(f"Unsupported model: {config['model']['type']}. Use 'lisflood-fp' or 'lisflood82'")
    
    # Load observations
    print("\n📊 Loading observation data...")
    
    # Load gauge data
    gauge_config = config["observations"]["gauges"][0]
    gauge_data = load_observations(gauge_config["path"], "gauge")
    
    # Load satellite inundation
    satellite_config = config["observations"]["satellite"][0]
    
    # Auto-extract via GEE if enabled
    if "auto_extract" in satellite_config and satellite_config["auto_extract"]["enabled"]:
        print("   Extracting flood mask from Google Earth Engine...")
        # (GEE integration placeholder)
        # inundation_mask = extract_from_gee(satellite_config["auto_extract"])
        print("   ⚠️  GEE extraction not yet implemented")
        print(f"   Using existing file: {satellite_config['path']}")
    
    # Load satellite data
    inundation_mask = load_observations(satellite_config["path"], "satellite")
    
    # Prepare storm event
    storm_event = {
        "start_date": config["input_data"]["rainfall"]["start_date"],
        "end_date": config["input_data"]["rainfall"]["end_date"],
        "rainfall_data": config["input_data"]["rainfall"]["path"],
        "boundary_discharge": config["input_data"]["boundary_discharge"],
        "observations": {
            "gauge": gauge_data,
            "satellite": {"mask": inundation_mask}
        }
    }
    
    # Initialize calibration loop
    print("\n🔄 Initializing calibration loop...")
    loop = AutoCalibrationLoop(
        model_adapter=model,
        storm_event=storm_event,
        target_metrics=config["calibration"]["target_metrics"],
        weights=config["calibration"]["objective_weights"]
    )
    
    # Run calibration
    max_iterations = args.max_iterations or config["calibration"]["max_iterations"]
    
    print(f"\n🚀 Starting calibration (max {max_iterations} iterations)...")
    results = loop.run(max_iterations=max_iterations, verbose=args.verbose)
    
    # Save results
    output_dir = config["model"]["output_dir"]
    results_path = os.path.join(output_dir, "calibration_results.json")
    
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n📁 Calibration results saved to: {results_path}")
    
    # Save best parameters to .par file
    best_params_path = os.path.join(output_dir, "calibrated.par")
    model.save_config(best_params_path)
    
    # Generate plots if requested
    if args.save_plots:
        print("\n📈 Generating visualization plots...")
        from visualization.plots import (
            plot_calibration_progress,
            plot_pareto_front,
            plot_hydrograph_comparison,
            plot_inundation_comparison
        )
        
        plot_calibration_progress(
            results["history"],
            os.path.join("visualization", "calibration_progress.png")
        )
        
        plot_pareto_front(
            results["pareto_front"],
            os.path.join("visualization", "pareto_front.png")
        )
        
        plot_hydrograph_comparison(
            storm_event["observations"]["gauge"],
            model.get_outputs()["discharge"],
            os.path.join("visualization", "hydrograph_comparison.png")
        )
        
        plot_inundation_comparison(
            storm_event["observations"]["satellite"]["mask"],
            model.get_outputs()["inundation"]["grid"],
            os.path.join("visualization", "inundation_comparison.png")
        )
        
        print("   Plots saved to visualization/")
    
    # Summary
    print("\n" + "="*60)
    print("📊 CALIBRATION SUMMARY")
    print("="*60)
    print(f"Iterations: {results['iterations']}")
    print(f"Converged: {results['converged']}")
    print(f"\nBest Metrics:")
    print(f"  NSE: {results['best_metrics']['nse']:.4f} (target: {config['calibration']['target_metrics']['nse']})")
    print(f"  KGE: {results['best_metrics']['kge']:.4f} (target: {config['calibration']['target_metrics']['kge']})")
    print(f"  IoU: {results['best_metrics']['iou']:.4f} (target: {config['calibration']['target_metrics']['iou']})")
    print(f"\nBest Score: {results['best_score']:.4f}")
    print(f"Pareto Front: {len(results['pareto_front'])} solutions")
    print(f"Lessons Learned: {len(results['lessons'])}")
    
    # Export lessons to MetaClaw if enabled
    if config.get("metaclaw", {}).get("enabled", False):
        print("\n🧠 Exporting lessons to MetaClaw...")
        # (MetaClaw integration placeholder)
        # export_lessons_to_metaclaw(results["lessons"], config["metaclaw"])
        print("   ⚠️  MetaClaw integration not yet implemented")
    
    print("\n✅ Calibration complete!")


if __name__ == "__main__":
    main()
