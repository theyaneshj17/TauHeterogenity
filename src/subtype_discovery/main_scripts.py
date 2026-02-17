import os
import json
from datetime import datetime
import itertools
import numpy as np
from sklearn.model_selection import KFold
import concurrent.futures
from typing import List, Tuple, Dict, Any
from dataclasses import dataclass
import pandas as pd
import matplotlib.pyplot as plt
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("ncv_execution.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


@dataclass
class InnerFoldTask:
    """Data class for inner fold tasks."""
    params: tuple
    train_idx: np.ndarray
    val_idx: np.ndarray
    fold_idx: int
    outer_fold: int
    data_path: str
    output_dir: str
    script_path: str


@dataclass
class OuterFoldTask:
    """Data class for outer fold tasks."""
    params: tuple
    train_idx: np.ndarray
    test_idx: np.ndarray
    outer_fold: int
    data_path: str
    output_dir: str
    script_path: str


def get_actual_dataset_size(data_path: str) -> int:
    """Determine the actual number of samples in the dataset."""
    try:
        if data_path.endswith('.csv'):
            df = pd.read_csv(data_path)
            return len(df)
        else:
            logger.warning(f"Unknown file type for {data_path}, using default size")
            return 318
    except Exception as e:
        logger.error(f"Error reading dataset size: {e}")
        return 318


def find_best_hyperparams_by_k(output_dir: str, param_combinations: List[tuple],
                               metrics: List[str] = ['aic', 'bic']) -> Dict[int, Dict[str, Any]]:
    """
    Select optimal learning rate and reassignment frequency for each k based on validation loss.
    
    Following methodology: "Hyperparameters such as learning rate (lr) and iteration counts 
    between cluster assignments (re) were selected for each k based on averaged cross-entropy loss."
    
    Returns best hyperparameters for EACH k value.
    """
    results = []
    
    # Collect all inner fold results
    for file in os.listdir(output_dir):
        if file.startswith('inner_fold_results_') and file.endswith('.json'):
            try:
                with open(os.path.join(output_dir, file), 'r') as f:
                    results.append(json.load(f))
            except Exception as e:
                logger.warning(f"Could not read {file}: {e}")
    
    if not results:
        logger.warning("No inner fold results found")
        return {}
    
    # Aggregate results by parameters
    aggregated = {}
    for result in results:
        params = (
            result['parameters']['num_cluster'],
            result['parameters']['lr'],
            result['parameters']['reassign']
        )
        if params not in aggregated:
            aggregated[params] = {metric: [] for metric in metrics}
            aggregated[params]['val_loss'] = []
        
        for metric in metrics:
            if metric in result['metrics']:
                metric_value = result['metrics'][metric]
                if np.isfinite(metric_value):
                    aggregated[params][metric].append(metric_value)
        
        if 'val_loss' in result['metrics']:
            val_loss = result['metrics']['val_loss']
            if np.isfinite(val_loss):
                aggregated[params]['val_loss'].append(val_loss)
    
    # Calculate mean for each parameter combination
    mean_metrics = {}
    for params, values in aggregated.items():
        mean_metrics[params] = {}
        for metric in metrics + ['val_loss']:
            if values[metric]:
                mean_metrics[params][metric] = np.mean(values[metric])
            else:
                mean_metrics[params][metric] = float('inf')
        
        logger.info(f"Parameters {params}: val_loss={mean_metrics[params]['val_loss']:.4f}")
    
    # For each k, find best lr and reassign based on validation loss
    k_values = sorted(set(params[0] for params in mean_metrics.keys()))
    best_params_by_k = {}
    
    for k in k_values:
        k_params = {params: metrics for params, metrics in mean_metrics.items() if params[0] == k}
        
        if k_params:
            valid_params = {p: m for p, m in k_params.items() if np.isfinite(m['val_loss'])}
            
            if valid_params:
                best_param = min(valid_params.keys(), key=lambda x: valid_params[x]['val_loss'])
                best_params_by_k[k] = {
                    'params': best_param,
                    'val_loss': valid_params[best_param]['val_loss']
                }
                
                for metric in metrics:
                    if np.isfinite(valid_params[best_param].get(metric, float('inf'))):
                        best_params_by_k[k][metric] = valid_params[best_param][metric]
                
                logger.info(f"k={k}: best lr={best_param[1]}, reassign={best_param[2]} "
                           f"(val_loss={valid_params[best_param]['val_loss']:.4f})")
            else:
                logger.warning(f"No valid parameter combinations for k={k}")
    
    # Create visualizations
    try:
        k_values = sorted(set(params[0] for params in aggregated.keys()))
        
        if len(k_values) > 1:
            plot_data = {metric: {k: [] for k in k_values} for metric in metrics + ['val_loss']}
            
            for params, values in aggregated.items():
                k = params[0]
                for metric in metrics + ['val_loss']:
                    if values[metric]:
                        plot_data[metric][k].extend(values[metric])
            
            # Box plots for each metric
            for metric in metrics + ['val_loss']:
                plt.figure(figsize=(10, 6))
                data_to_plot = [plot_data[metric][k] for k in k_values if plot_data[metric][k]]
                labels = [f'k={k}' for k in k_values if plot_data[metric][k]]
                
                if data_to_plot:
                    plt.boxplot(data_to_plot, labels=labels)
                    plt.title(f'{metric.upper()} by Number of Clusters', fontsize=24, fontweight='bold')
                    plt.ylabel(metric.upper(), fontsize=20, fontweight='bold')
                    plt.xlabel('Number of Clusters (k)', fontsize=20, fontweight='bold')
                    plt.xticks(fontsize=16, fontweight='bold')
                    plt.yticks(fontsize=16, fontweight='bold')
                    plt.grid(True, linestyle='--', alpha=0.7)
                    plt.savefig(os.path.join(output_dir, f'{metric}_by_k.png'), dpi=300, bbox_inches='tight')
                    plt.close()
    except Exception as e:
        logger.warning(f"Error creating plots: {e}")
    
    return best_params_by_k


def process_inner_fold(task: InnerFoldTask) -> Tuple[float, float, float, tuple]:
    """Process a single inner fold task."""
    num_cluster, lr, reassign = task.params
    
    # Filter indices to be within dataset size
    actual_size = get_actual_dataset_size(task.data_path)
    valid_train_idx = [idx for idx in task.train_idx if idx < actual_size]
    valid_val_idx = [idx for idx in task.val_idx if idx < actual_size]
    
    if len(valid_train_idx) < len(task.train_idx) or len(valid_val_idx) < len(task.val_idx):
        logger.warning(f"Filtered indices for inner fold {task.fold_idx}: "
                      f"Train {len(task.train_idx)}->{len(valid_train_idx)}, "
                      f"Val {len(task.val_idx)}->{len(valid_val_idx)}")
    
    train_indices_str = ','.join(map(str, valid_train_idx))
    val_indices_str = ','.join(map(str, valid_val_idx))
    
    run_id = f"outer_fold_{task.outer_fold}_inner_fold_{task.fold_idx}_k{num_cluster}_lr{lr}_reassign{reassign}"
    logger.info(f"Starting inner fold: {run_id}")
    
    # Construct command for inner fold
    command = (
        f"python {task.script_path} "
        f"--data {task.data_path} "
        f"--lr {lr} "
        f"--num_cluster {num_cluster} "
        f"--reassign {reassign} "
        f"--fold_idx {task.fold_idx} "
        f"--train_indices {train_indices_str} "
        f"--val_indices {val_indices_str} "
        f"--t 1 "  # Inner fold indicator
        f"--outer_fold {task.outer_fold} "
        f"--output_dir {task.output_dir} "
        f"--epochs 2"
    )
    
    try:
        logger.info(f"Executing: {command}")
        return_code = os.system(command)
        
        if return_code != 0:
            logger.error(f"Command failed with return code {return_code}")
            raise RuntimeError(f"Command failed")
        
        # Load results
        profile = f"train_outer_fold={task.outer_fold}_inner_fold={task.fold_idx}_type=1_k={num_cluster}_lr={lr}_reassign={reassign}"
        results_dir = os.path.join(task.output_dir, f'Results_{profile}')
        results_file = os.path.join(results_dir, 'results.json')
        
        training_results = {}
        if os.path.exists(results_file):
            with open(results_file, 'r') as f:
                training_results = json.load(f)
        
        # Extract validation loss
        val_losses = training_results.get('val_loss', [float('inf')])
        val_loss = np.mean(val_losses) if val_losses else float('inf')
        
        # Load GMM parameters for AIC/BIC
        aic, bic = float('inf'), float('inf')
        gmm_params_file = os.path.join(results_dir, 'gmm_params.json')
        
        if os.path.exists(gmm_params_file):
            try:
                with open(gmm_params_file, 'r') as f:
                    gmm_params = json.load(f)
                if 'aic' in gmm_params and 'bic' in gmm_params:
                    aic, bic = gmm_params['aic'], gmm_params['bic']
                    logger.info(f"Using AIC/BIC from GMM: AIC={aic:.4f}, BIC={bic:.4f}")
            except Exception as e:
                logger.warning(f"Error loading GMM parameters: {e}")
        else:
            logger.warning(f"No GMM parameters found at {gmm_params_file}")
    
    except Exception as e:
        logger.error(f"Error processing {run_id}: {e}")
        val_loss = float('inf')
        aic = float('inf')
        bic = float('inf')
    
    # Save inner fold results
    inner_fold_results = {
        'run_id': run_id,
        'parameters': {
            'num_cluster': num_cluster,
            'lr': lr,
            'reassign': reassign,
            'outer_fold': task.outer_fold
        },
        'fold_idx': task.fold_idx,
        'metrics': {
            'val_loss': val_loss,
            'aic': float(aic),
            'bic': float(bic)
        },
        'num_samples': len(valid_val_idx),
        'timestamp': datetime.now().isoformat()
    }
    
    results_filename = os.path.join(task.output_dir, f'inner_fold_results_{run_id}.json')
    with open(results_filename, 'w') as f:
        json.dump(inner_fold_results, f, indent=4)
    
    logger.info(f"Completed {run_id}: val_loss={val_loss:.4f}, AIC={aic:.4f}, BIC={bic:.4f}")
    return val_loss, aic, bic, task.params


def process_outer_fold_task(task: OuterFoldTask) -> Dict[str, Any]:
    """Process a single outer fold task with specific hyperparameters."""
    num_cluster, lr, reassign = task.params
    
    # Filter indices
    actual_size = get_actual_dataset_size(task.data_path)
    valid_train_idx = [idx for idx in task.train_idx if idx < actual_size]
    valid_test_idx = [idx for idx in task.test_idx if idx < actual_size]
    
    if len(valid_train_idx) < len(task.train_idx) or len(valid_test_idx) < len(task.test_idx):
        logger.warning(f"Filtered indices for outer fold {task.outer_fold}: "
                      f"Train {len(task.train_idx)}->{len(valid_train_idx)}, "
                      f"Test {len(task.test_idx)}->{len(valid_test_idx)}")
    
    train_indices_str = ','.join(map(str, valid_train_idx))
    test_indices_str = ','.join(map(str, valid_test_idx))
    
    run_id = f"outer_fold_{task.outer_fold}_k{num_cluster}_lr{lr}_reassign{reassign}"
    logger.info(f"Starting outer fold evaluation: {run_id}")
    
    # Construct command for outer fold
    command = (
        f"python {task.script_path} "
        f"--data {task.data_path} "
        f"--lr {lr} "
        f"--num_cluster {num_cluster} "
        f"--reassign {reassign} "
        f"--fold_idx {task.outer_fold} "
        f"--train_indices {train_indices_str} "
        f"--val_indices {test_indices_str} "
        f"--t 0 "  # Outer fold indicator
        f"--outer_fold {task.outer_fold} "
        f"--output_dir {task.output_dir}"
    )
    
    try:
        logger.info(f"Executing: {command}")
        return_code = os.system(command)
        if return_code != 0:
            logger.error(f"Training failed with return code {return_code}")
            raise RuntimeError(f"Command failed")
        
        # Load results
        profile = f"train_outer_fold={task.outer_fold}_type=0_k={num_cluster}_lr={lr}_reassign={reassign}"
        results_dir = os.path.join(task.output_dir, f'Results_{profile}')
        results_file = os.path.join(results_dir, 'results.json')
        
        if os.path.exists(results_file):
            with open(results_file, 'r') as f:
                results = json.load(f)
        else:
            logger.error(f"Results file not found: {results_file}")
            results = {
                'train_loss': [0.0],
                'val_loss': [0.0],
                'train_accuracy': [0.0],
                'val_accuracy': [0.0],
                'train_ami': [0.0],
                'val_ami': [0.0]
            }
        
        # Load GMM parameters
        gmm_params = None
        gmm_params_file = os.path.join(results_dir, 'gmm_params.json')
        if os.path.exists(gmm_params_file):
            try:
                with open(gmm_params_file, 'r') as f:
                    gmm_params = json.load(f)
            except Exception as e:
                logger.warning(f"Error loading GMM parameters: {e}")
        
        # Extract results
        evaluation_results = {
            'params': task.params,
            'outer_fold': task.outer_fold,
            'num_cluster': num_cluster,
            'lr': lr,
            'reassign': reassign,
            'train_loss': results.get('train_loss', [0.0])[-1] if results.get('train_loss') else 0.0,
            'val_loss': results.get('val_loss', [0.0])[-1] if results.get('val_loss') else 0.0,
            'train_accuracy': results.get('train_accuracy', [0.0])[-1] if results.get('train_accuracy') else 0.0,
            'val_accuracy': results.get('val_accuracy', [0.0])[-1] if results.get('val_accuracy') else 0.0,
            'train_ami': results.get('train_ami', [0.0])[-1] if results.get('train_ami') else 0.0,
            'val_ami': results.get('val_ami', [0.0])[-1] if results.get('val_ami') else 0.0,
        }
        
        # Add AIC/BIC if available
        if gmm_params:
            if 'aic' in gmm_params:
                evaluation_results['aic'] = gmm_params['aic']
            if 'bic' in gmm_params:
                evaluation_results['bic'] = gmm_params['bic']
        
        # Save individual outer fold results
        results_filename = os.path.join(task.output_dir, f'outer_fold_results_{run_id}.json')
        with open(results_filename, 'w') as f:
            json.dump(evaluation_results, f, indent=4)
        
        logger.info(f"Completed {run_id}")
        return evaluation_results
    
    except Exception as e:
        logger.error(f"Error processing {run_id}: {e}")
        return {
            'params': task.params,
            'outer_fold': task.outer_fold,
            'num_cluster': num_cluster,
            'lr': lr,
            'reassign': reassign,
            'error': str(e),
            'train_loss': float('inf'),
            'val_loss': float('inf'),
            'train_accuracy': 0.0,
            'val_accuracy': 0.0,
            'train_ami': 0.0,
            'val_ami': 0.0
        }


def nested_cross_validation(data_path: str, script_path: str, output_dir: str = '',
                           n_outer_folds: int = 10, n_inner_folds: int = 2,
                           max_workers: int = 4) -> Dict[str, Any]:
    """
    Implement parallel nested cross-validation with hyperparameter tuning.
    
    For each outer fold:
    1. Run inner CV to find best hyperparameters for each k
    2. Evaluate each k with its best hyperparameters on outer test set
    """
    n_inner_folds = max(2, n_inner_folds)
    n_outer_folds = max(2, n_outer_folds)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Get dataset size
    try:
        if data_path.endswith('.csv'):
            df = pd.read_csv(data_path)
            df = df.dropna()
            n_samples = len(df)
        else:
            n_samples = get_actual_dataset_size(data_path)
    except Exception as e:
        logger.error(f"Error reading dataset: {e}")
        n_samples = 318
    
    logger.info(f"Dataset contains {n_samples} samples")
    
    # Hyperparameter grid
    hyperparameter_grid = {
        'num_cluster': [3],
        'lr': [ 0.05],
        'reassign': [1.0]
    }
    
    param_combinations = list(itertools.product(*hyperparameter_grid.values()))
    logger.info(f"Testing {len(param_combinations)} parameter combinations "
               f"across {n_outer_folds} outer folds with {n_inner_folds} inner folds each")
    
    # Set up cross-validation
    outer_cv = KFold(n_splits=n_outer_folds, shuffle=True, random_state=31)
    inner_cv = KFold(n_splits=n_inner_folds, shuffle=True, random_state=31)
    
    all_indices = np.arange(n_samples)
    all_outer_results = []
    
    try:
        # Outer cross-validation loop
        for outer_fold, (train_idx, test_idx) in enumerate(outer_cv.split(all_indices)):
            train_indices = all_indices[train_idx]
            test_indices = all_indices[test_idx]
            
            logger.info(f"\n{'='*60}")
            logger.info(f"Outer Fold {outer_fold + 1}/{n_outer_folds}")
            logger.info(f"Train size: {len(train_indices)}, Test size: {len(test_indices)}")
            logger.info(f"{'='*60}")
            
            # Inner cross-validation for hyperparameter tuning
            inner_fold_tasks = []
            
            for fold_idx, (inner_train_idx, inner_val_idx) in enumerate(inner_cv.split(train_indices)):
                inner_train_indices = train_indices[inner_train_idx]
                inner_val_indices = train_indices[inner_val_idx]
                
                logger.info(f"Inner fold {fold_idx}: Train={len(inner_train_indices)}, Val={len(inner_val_indices)}")
                
                for params in param_combinations:
                    task = InnerFoldTask(
                        params=params,
                        train_idx=inner_train_indices,
                        val_idx=inner_val_indices,
                        fold_idx=fold_idx,
                        outer_fold=outer_fold,
                        data_path=data_path,
                        output_dir=output_dir,
                        script_path=script_path
                    )
                    inner_fold_tasks.append(task)
            
            # Execute inner fold tasks in parallel
            with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(process_inner_fold, task) for task in inner_fold_tasks]
                
                for future in concurrent.futures.as_completed(futures):
                    try:
                        result = future.result()
                        logger.info(f"Completed inner fold: {result[3]}")
                    except Exception as e:
                        logger.error(f"Inner fold task failed: {e}")
            
            # Find best hyperparameters for each k
            best_params_by_k = find_best_hyperparams_by_k(output_dir, param_combinations)
            
            logger.info(f"Best hyperparameters for each k in outer fold {outer_fold}:")
            for k, info in best_params_by_k.items():
                logger.info(f"  k={k}: lr={info['params'][1]}, reassign={info['params'][2]} "
                           f"(val_loss={info['val_loss']:.4f})")
            
            # Create outer fold tasks for each k
            outer_fold_tasks = []
            for k, params_info in best_params_by_k.items():
                best_params = params_info['params']
                task = OuterFoldTask(
                    params=best_params,
                    train_idx=train_indices,
                    test_idx=test_indices,
                    outer_fold=outer_fold,
                    data_path=data_path,
                    output_dir=output_dir,
                    script_path=script_path
                )
                outer_fold_tasks.append(task)
            
            # Execute outer fold tasks
            with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, max_workers // 2)) as executor:
                futures = [executor.submit(process_outer_fold_task, task) for task in outer_fold_tasks]
                
                for future in concurrent.futures.as_completed(futures):
                    try:
                        result = future.result()
                        all_outer_results.append(result)
                        logger.info(f"Completed outer fold: k={result['num_cluster']}, "
                                  f"val_loss={result['val_loss']:.4f}, "
                                  f"val_accuracy={result['val_accuracy']:.4f}")
                    except Exception as e:
                        logger.error(f"Outer fold task failed: {e}")
            
            # Save progress
            with open(os.path.join(output_dir, 'nested_cv_results.json'), 'w') as f:
                json.dump({
                    'outer_results': all_outer_results,
                    'hyperparameter_grid': hyperparameter_grid,
                    'timestamp': datetime.now().isoformat(),
                }, f, indent=4)
        
        # Generate final summary
        summary = analyze_cv_results_by_k(all_outer_results, output_dir)
        
        return {
            'outer_results': all_outer_results,
            'summary': summary
        }
        
    except Exception as e:
        logger.error(f"Execution stopped due to error: {str(e)}")
        
        if all_outer_results:
            with open(os.path.join(output_dir, 'nested_cv_results_partial.json'), 'w') as f:
                json.dump({
                    'outer_results': all_outer_results,
                    'hyperparameter_grid': hyperparameter_grid,
                    'error': str(e),
                    'timestamp': datetime.now().isoformat(),
                }, f, indent=4)
        raise


def analyze_cv_results_by_k(results: List[Dict[str, Any]], output_dir: str) -> Dict[str, Any]:
    """Analyze cross-validation results and create visualizations."""
    if not results:
        return {}
    
    metrics = ['train_loss', 'val_loss', 'train_accuracy', 'val_accuracy', 'train_ami', 'val_ami']
    if results and 'aic' in results[0]:
        metrics.extend(['aic', 'bic'])
    
    # Group results by k
    results_by_k = {}
    for result in results:
        k = result['num_cluster']
        if k not in results_by_k:
            results_by_k[k] = []
        results_by_k[k].append(result)
    
    # Calculate summary statistics for each k
    summary_by_k = {}
    for k, k_results in results_by_k.items():
        data = {metric: [] for metric in metrics}
        
        for result in k_results:
            for metric in metrics:
                if metric in result:
                    data[metric].append(result[metric])
        
        summary_by_k[k] = {
            'mean': {metric: np.mean(values) for metric, values in data.items() if values},
            'std': {metric: np.std(values) for metric, values in data.items() if values},
            'min': {metric: np.min(values) for metric, values in data.items() if values},
            'max': {metric: np.max(values) for metric, values in data.items() if values},
            'count': len(k_results)
        }
    
    logger.info("Summary statistics by k:")
    for k, stats in summary_by_k.items():
        logger.info(f"\nk={k} (evaluated in {stats['count']} outer folds):")
        for metric in metrics:
            if metric in stats['mean']:
                logger.info(f"  {metric}: Mean={stats['mean'][metric]:.4f}, Std={stats['std'][metric]:.4f}")
    
    # Find best k by metric
    best_k_by_metric = {}
    for metric in metrics:
        valid_k = [k for k, stats in summary_by_k.items() if metric in stats['mean']]
        
        if valid_k:
            if metric in ['val_loss', 'aic', 'bic']:
                best_k = min(valid_k, key=lambda k: summary_by_k[k]['mean'][metric])
            else:
                best_k = max(valid_k, key=lambda k: summary_by_k[k]['mean'][metric])
            
            best_k_by_metric[metric] = {
                'k': best_k,
                'value': summary_by_k[best_k]['mean'][metric]
            }
            
            logger.info(f"Best k by {metric}: k={best_k} (mean={summary_by_k[best_k]['mean'][metric]:.4f})")
    
    # Create visualizations
    try:
        # Performance comparison across k
        plt.figure(figsize=(15, 10))
        n_metrics = len(metrics)
        n_cols = min(3, n_metrics)
        n_rows = (n_metrics + n_cols - 1) // n_cols
        
        for i, metric in enumerate(metrics):
            ax = plt.subplot(n_rows, n_cols, i + 1)
            
            k_values = []
            mean_values = []
            std_values = []
            
            for k, stats in summary_by_k.items():
                if metric in stats['mean']:
                    k_values.append(k)
                    mean_values.append(stats['mean'][metric])
                    std_values.append(stats['std'][metric])
            
            if k_values:
                sort_idx = np.argsort(k_values)
                k_values = [k_values[i] for i in sort_idx]
                mean_values = [mean_values[i] for i in sort_idx]
                std_values = [std_values[i] for i in sort_idx]
                
                ax.errorbar(k_values, mean_values, yerr=std_values, fmt='o-', capsize=5)
                
                if metric in best_k_by_metric:
                    best_k = best_k_by_metric[metric]['k']
                    best_value = best_k_by_metric[metric]['value']
                    ax.plot(best_k, best_value, 'r*', markersize=15, label=f'Best k={best_k}')
                
                ax.set_title(metric.replace('_', ' ').title(), fontsize=24, fontweight='bold')
                ax.set_xlabel('Number of Clusters (k)', fontsize=20, fontweight='bold')
                ax.set_ylabel(metric.replace('_', ' ').title(), fontsize=20, fontweight='bold')
                ax.grid(True, linestyle='--', alpha=0.7)
                ax.set_xticks(k_values)
                ax.tick_params(axis='both', which='major', labelsize=16)
                
                if metric in best_k_by_metric:
                    ax.legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'all_metrics_by_k.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Box plots
        for metric in metrics:
            data_to_plot = []
            labels = []
            
            for k in sorted(results_by_k.keys()):
                values = [result[metric] for result in results_by_k[k] if metric in result]
                if values:
                    data_to_plot.append(values)
                    labels.append(f'k={k}')
            
            if data_to_plot:
                plt.figure(figsize=(10, 6))
                plt.boxplot(data_to_plot, labels=labels)
                plt.title(f'{metric.replace("_", " ").title()} by Number of Clusters',
                         fontsize=24, fontweight='bold')
                plt.ylabel(metric.replace('_', ' ').title(), fontsize=20, fontweight='bold')
                plt.xticks(fontsize=16, fontweight='bold')
                plt.yticks(fontsize=16, fontweight='bold')
                plt.grid(True, linestyle='--', alpha=0.7)
                plt.savefig(os.path.join(output_dir, f'{metric}_boxplot_by_k.png'),
                           dpi=300, bbox_inches='tight')
                plt.close()
        
        # Summary table
        summary_data = []
        for k in sorted(summary_by_k.keys()):
            row = {'k': k}
            for metric in metrics:
                if metric in summary_by_k[k]['mean']:
                    row[f'{metric}_mean'] = summary_by_k[k]['mean'][metric]
                    row[f'{metric}_std'] = summary_by_k[k]['std'][metric]
            summary_data.append(row)
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(os.path.join(output_dir, 'summary_by_k.csv'), index=False)
        
    except Exception as e:
        logger.warning(f"Error creating visualizations: {e}")
    
    # Create report
    with open(os.path.join(output_dir, 'nested_cv_analysis_report.txt'), 'w') as f:
        f.write("NESTED CROSS-VALIDATION ANALYSIS REPORT\n")
        f.write("="*60 + "\n\n")
        
        f.write("SUMMARY BY NUMBER OF CLUSTERS\n")
        f.write("-"*60 + "\n")
        for k in sorted(summary_by_k.keys()):
            f.write(f"\nk={k} (evaluated in {summary_by_k[k]['count']} outer folds)\n\n")
            for metric in metrics:
                if metric in summary_by_k[k]['mean']:
                    f.write(f"{metric}: {summary_by_k[k]['mean'][metric]:.4f} "
                           f"± {summary_by_k[k]['std'][metric]:.4f}\n")
        
        f.write("\n\nBEST K BY METRIC\n")
        f.write("-"*60 + "\n")
        for metric in metrics:
            if metric in best_k_by_metric:
                f.write(f"{metric}: k={best_k_by_metric[metric]['k']} "
                       f"(value={best_k_by_metric[metric]['value']:.4f})\n")
    
    # Determine overall best k
    final_summary = {
        'summary_by_k': summary_by_k,
        'best_k_by_metric': best_k_by_metric
    }
    
    if 'val_ami' in best_k_by_metric:
        final_summary['best_k'] = best_k_by_metric['val_ami']['k']
        final_summary['best_k_metric'] = 'val_ami'
    elif 'val_accuracy' in best_k_by_metric:
        final_summary['best_k'] = best_k_by_metric['val_accuracy']['k']
        final_summary['best_k_metric'] = 'val_accuracy'
    elif 'bic' in best_k_by_metric:
        final_summary['best_k'] = best_k_by_metric['bic']['k']
        final_summary['best_k_metric'] = 'bic'
    
    with open(os.path.join(output_dir, 'final_summary.json'), 'w') as f:
        json.dump(final_summary, f, indent=4)
    
    return final_summary


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Nested cross-validation for clustering')
    parser.add_argument('--data_path', type=str, required=True, help='Path to dataset CSV')
    parser.add_argument('--script_path', type=str, required=True, help='Path to main.py')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory')
    parser.add_argument('--outer_folds', type=int, default=10, help='Number of outer folds')
    parser.add_argument('--inner_folds', type=int, default=2, help='Number of inner folds')
    parser.add_argument('--max_workers', type=int, default=4, help='Max parallel workers')
    
    args = parser.parse_args()
    
    logger.info("Starting nested cross-validation")
    logger.info(f"Data: {args.data_path}")
    logger.info(f"Script: {args.script_path}")
    logger.info(f"Output: {args.output_dir}")
    logger.info(f"Outer folds: {args.outer_folds}, Inner folds: {args.inner_folds}")
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Save configuration
    with open(os.path.join(args.output_dir, 'config.json'), 'w') as f:
        json.dump({
            'data_path': args.data_path,
            'script_path': args.script_path,
            'output_dir': args.output_dir,
            'outer_folds': args.outer_folds,
            'inner_folds': args.inner_folds,
            'max_workers': args.max_workers,
            'timestamp': datetime.now().isoformat()
        }, f, indent=4)
    
    try:
        results = nested_cross_validation(
            data_path=args.data_path,
            script_path=args.script_path,
            output_dir=args.output_dir,
            n_outer_folds=args.outer_folds,
            n_inner_folds=args.inner_folds,
            max_workers=args.max_workers
        )
        
        logger.info("Nested cross-validation completed successfully!")
        
        if 'summary' in results and 'best_k' in results['summary']:
            logger.info(f"Best k: {results['summary']['best_k']}")
            logger.info(f"Selected by: {results['summary']['best_k_metric']}")
        
    except Exception as e:
        logger.error(f"Failed: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())