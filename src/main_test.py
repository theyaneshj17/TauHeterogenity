import os
import json
import glob
import logging
from datetime import datetime
from typing import List, Dict, Any
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("test_execution.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def find_trained_models(base_dir: str, pattern: str = "Results_train_outer_fold=*_type=0_*") -> List[Dict[str, str]]:
    """Find all trained models matching the pattern."""
    search_pattern = os.path.join(base_dir, pattern)
    model_dirs = glob.glob(search_pattern)
    model_dirs.sort()
    
    models = []
    for model_dir in model_dirs:
        model_path = os.path.join(model_dir, "model.pth.tar")
        gmm_path_mat = os.path.join(model_dir, "gm_params.mat")
        gmm_path_json = os.path.join(model_dir, "gmm_params.json")
        
        # Prefer .mat file, fall back to .json
        gmm_path = gmm_path_mat if os.path.exists(gmm_path_mat) else gmm_path_json
        
        if os.path.exists(model_path) and os.path.exists(gmm_path):
            # Extract metadata from directory name
            dir_name = os.path.basename(model_dir)
            parts = dir_name.split('_')
            
            metadata = {
                'model_dir': model_dir,
                'model_path': model_path,
                'gmm_path': gmm_path,
                'dir_name': dir_name
            }
            
            # Try to extract parameters
            for part in parts:
                if '=' in part:
                    key, value = part.split('=')
                    metadata[key] = value
            
            models.append(metadata)
            logger.info(f"Found model: {dir_name}")
        else:
            logger.warning(f"Incomplete model in {model_dir}")
    
    return models


def test_single_model(model_info: Dict[str, str], data_path: str, script_path: str, 
                     output_base_dir: str) -> Dict[str, Any]:
    """Test a single model."""
    model_name = model_info['dir_name']
    output_dir = os.path.join(output_base_dir, f"Test_{model_name}")
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Testing: {model_name}")
    logger.info(f"{'='*60}")
    
    # Construct test command
    command = (
        f"python {script_path} "
        f"--mode test "
        f"--model {model_info['model_path']} "
        f"--gmmmodel {model_info['gmm_path']} "
        f"--data ADNI "
        f"--output_dir {output_dir} "
        f"--batch 32 "
        f"--workers 0"
    )
    
    try:
        logger.info(f"Executing: {command}")
        return_code = os.system(command)
        
        if return_code != 0:
            logger.error(f"Test failed with return code {return_code}")
            return {
                'model_name': model_name,
                'status': 'failed',
                'error': f'Return code {return_code}'
            }
        
        # Load test results
        metrics_file = os.path.join(output_dir, 'test_metrics.json')
        if os.path.exists(metrics_file):
            with open(metrics_file, 'r') as f:
                metrics = json.load(f)
            
            result = {
                'model_name': model_name,
                'status': 'success',
                'output_dir': output_dir,
                **model_info,
                **metrics
            }
            
            logger.info(f"Test completed: Accuracy={metrics['accuracy']:.2f}%, Loss={metrics['loss']:.4f}")
            return result
        else:
            logger.error(f"Metrics file not found: {metrics_file}")
            return {
                'model_name': model_name,
                'status': 'failed',
                'error': 'Metrics file not found'
            }
    
    except Exception as e:
        logger.error(f"Error testing model: {e}")
        return {
            'model_name': model_name,
            'status': 'failed',
            'error': str(e)
        }


def aggregate_test_results(results: List[Dict[str, Any]], output_dir: str):
    """Aggregate and save test results."""
    logger.info("\n" + "="*60)
    logger.info("AGGREGATING RESULTS")
    logger.info("="*60)
    
    # Create summary DataFrame
    summary_data = []
    for result in results:
        if result['status'] == 'success':
            summary_data.append({
                'Model': result['model_name'],
                'Outer_Fold': result.get('outer', 'N/A'),
                'K': result.get('k', 'N/A'),
                'LR': result.get('lr', 'N/A'),
                'Reassign': result.get('reassign', 'N/A'),
                'Accuracy': result.get('accuracy', 0.0),
                'Loss': result.get('loss', 0.0),
                'Num_Samples': result.get('num_samples', 0),
                'Num_Clusters': result.get('num_clusters', 0)
            })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(output_dir, 'test_summary.csv')
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"Summary saved to: {summary_file}")
        
        # Print statistics
        logger.info("\nTest Statistics:")
        logger.info(f"Total models tested: {len(results)}")
        logger.info(f"Successful: {len(summary_data)}")
        logger.info(f"Failed: {len(results) - len(summary_data)}")
        
        if len(summary_data) > 0:
            logger.info(f"\nMean Accuracy: {summary_df['Accuracy'].mean():.2f}%")
            logger.info(f"Std Accuracy: {summary_df['Accuracy'].std():.2f}%")
            logger.info(f"Mean Loss: {summary_df['Loss'].mean():.4f}")
            logger.info(f"Std Loss: {summary_df['Loss'].std():.4f}")
    
    # Save full results as JSON
    results_file = os.path.join(output_dir, 'test_results_full.json')
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=4)
    logger.info(f"Full results saved to: {results_file}")


def batch_test_models(trained_models_dir: str, data_path: str, script_path: str, 
                     output_dir: str, pattern: str = "Results_train_outer_fold=*_type=0_*"):
    """Test all trained models in batch."""
    logger.info("="*60)
    logger.info("BATCH TESTING MODELS")
    logger.info("="*60)
    logger.info(f"Trained models directory: {trained_models_dir}")
    logger.info(f"Output directory: {output_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Save configuration
    config = {
        'trained_models_dir': trained_models_dir,
        'data_path': data_path,
        'script_path': script_path,
        'output_dir': output_dir,
        'pattern': pattern,
        'timestamp': datetime.now().isoformat()
    }
    
    with open(os.path.join(output_dir, 'test_config.json'), 'w') as f:
        json.dump(config, f, indent=4)
    
    # Find all trained models
    models = find_trained_models(trained_models_dir, pattern)
    logger.info(f"\nFound {len(models)} trained models to test")
    
    if not models:
        logger.error("No models found!")
        return
    
    # Test each model
    results = []
    for i, model_info in enumerate(models):
        logger.info(f"\nTesting model {i+1}/{len(models)}")
        result = test_single_model(model_info, data_path, script_path, output_dir)
        results.append(result)
    
    # Aggregate results
    aggregate_test_results(results, output_dir)
    
    logger.info("\n" + "="*60)
    logger.info("BATCH TESTING COMPLETED")
    logger.info("="*60)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Batch test trained models')
    parser.add_argument('--trained_models_dir', type=str, required=True,
                        help='Directory containing trained models')
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to dataset CSV')
    parser.add_argument('--script_path', type=str, required=True,
                        help='Path to main.py')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for test results')
    parser.add_argument('--pattern', type=str, 
                        default='Results_train_outer_fold=*_type=0_*',
                        help='Pattern to match model directories')
    
    args = parser.parse_args()
    
    batch_test_models(
        trained_models_dir=args.trained_models_dir,
        data_path=args.data_path,
        script_path=args.script_path,
        output_dir=args.output_dir,
        pattern=args.pattern
    )