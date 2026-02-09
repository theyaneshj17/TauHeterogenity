import os
import json
import numpy as np
import torch
from torch.utils.data.sampler import Sampler
from itertools import permutations
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
from scipy import io
from sklearn import mixture
from sklearn.mixture._gaussian_mixture import _compute_precision_cholesky

import network


def label_stab(images_lists_unmatched, images_lists_prev, label_unmatched):
    """Stabilize labels across epochs for classification."""
    num_cluster = len(images_lists_unmatched)
    a = list(range(num_cluster))
    permute = list(permutations(a))
    num_same_subject = [0] * len(permute)
    
    for i in range(len(permute)):
        images_lists_perm = [[] for _ in range(num_cluster)]
        for j in range(num_cluster):
            images_lists_perm[j] = images_lists_unmatched[permute[i][j]]
        for k in range(num_cluster):
            num_same_subject[i] += len(set(images_lists_prev[k]) & set(images_lists_perm[k]))
    
    perm_idx = int(np.argmax(np.array(num_same_subject)))

    images_lists_now = [[] for _ in range(num_cluster)]
    for j in range(num_cluster):
        images_lists_now[j] = images_lists_unmatched[permute[perm_idx][j]]

    label_now = [0] * len(label_unmatched)
    for j in range(num_cluster):
        for k in range(len(images_lists_now[j])):
            label_now[images_lists_now[j][k]] = j

    return np.array(label_now), images_lists_now


def load_model(path):
    """Load a saved model from checkpoint."""
    if not os.path.isfile(path):
        print(f"Error: No model file found at '{path}'")
        return None

    try:
        print(f"Loading model from: {path}")
        checkpoint = torch.load(path, map_location=torch.device('cpu'))
        
        # Determine input dimension
        if 'state_dict' in checkpoint:
            state_dict = checkpoint['state_dict']
            if 'features.0.weight' in state_dict:
                input_dim = state_dict['features.0.weight'].shape[1]
                print(f"Detected input_dim={input_dim} from weights")
            else:
                input_dim = 68
                print(f"Using default input_dim={input_dim}")
        else:
            input_dim = 68
            print(f"Using default input_dim={input_dim}")
        
        # Get output dimension
        if 'state_dict' in checkpoint and 'top_layer.bias' in checkpoint['state_dict']:
            output_dim = checkpoint['state_dict']['top_layer.bias'].size()[0]
            print(f"Detected output_dim={output_dim} from top_layer.bias")
        else:
            try:
                folder_name = os.path.basename(os.path.dirname(path))
                if 'k=' in folder_name:
                    k_part = folder_name.split('k=')[1].split('_')[0]
                    output_dim = int(k_part)
                    print(f"Inferred output_dim={output_dim} from folder name")
                else:
                    output_dim = 2
                    print(f"Using default output_dim={output_dim}")
            except:
                output_dim = 2
                print(f"Using default output_dim={output_dim}")
        
        # Create model
        print(f"Creating model with: input_dim={input_dim}, output_dim={output_dim}")
        model = network.improved_mlp(
            input_dim=input_dim,
            output_dim=output_dim,
            dropout_rate=0.4,
            l2_lambda=1e-2,
            decoding=True
        )
        
        # Handle DataParallel state_dict
        def rename_key(key):
            if 'module' not in key:
                return key
            return ''.join(key.split('.module'))
            
        checkpoint['state_dict'] = {rename_key(key): val for key, val in checkpoint['state_dict'].items()}
        
        # Load weights
        model.load_state_dict(checkpoint['state_dict'])
        print("Model loaded successfully")
        return model
    
    except Exception as e:
        print(f"Error loading model: {e}")
        return None


def load_gmm(path):
    """Load GMM parameters from .mat file."""
    gm_struct = io.loadmat(path)['gm_params'][0, 0]

    weights = gm_struct['weights'].flatten()
    means = gm_struct['means']
    covariances = gm_struct['covariances']

    if isinstance(covariances, np.ndarray) and covariances.dtype == np.object_:
        covariances = covariances[0, 0] if covariances.shape == (1, 1) else covariances[0]

    covariances = np.asarray(covariances)
    assert covariances.ndim == 3, f"Covariance shape should be (n_components, n_features, n_features), got {covariances.shape}"

    n_components = int(gm_struct['n_components'][0, 0])

    # Build GMM
    clf = mixture.GaussianMixture(n_components=n_components, covariance_type='full')
    clf.weights_ = weights
    clf.means_ = means
    clf.covariances_ = covariances
    clf.precisions_cholesky_ = _compute_precision_cholesky(covariances, 'full')

    return clf


class UnifLabelSampler(Sampler):
    """Samples elements uniformly across pseudolabels."""
    
    def __init__(self, N, images_lists):
        self.N = N
        self.images_lists = images_lists
        self.indexes = self.generate_indexes_epoch()

    def generate_indexes_epoch(self):
        nmb_non_empty_clusters = 0
        for i in range(len(self.images_lists)):
            if len(self.images_lists[i]) != 0:
                nmb_non_empty_clusters += 1

        size_per_pseudolabel = int(self.N / nmb_non_empty_clusters) + 1
        res = np.array([])
        
        for i in range(len(self.images_lists)):
            if len(self.images_lists[i]) == 0:
                continue
            indexes = np.random.choice(
                self.images_lists[i],
                size_per_pseudolabel,
                replace=(len(self.images_lists[i]) <= size_per_pseudolabel)
            )
            res = np.concatenate((res, indexes))
        
        np.random.shuffle(res)
        res = list(res.astype('int'))
        if len(res) >= self.N:
            return res[:self.N]
        res += res[:(self.N - len(res))]
        return res
    
    def __iter__(self):
        return iter(self.indexes)
    
    def __len__(self):
        return len(self.indexes)


class AverageMeter(object):
    """Computes and stores the average and current value."""
    
    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count


def plotandsave(results_dir, results, train_index, test_index, gmm_params, model, optimizer, epoch):
    """Plot training metrics, create visualizations, and save results."""
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    
    # Set global font properties
    plt.rcParams['font.size'] = 20
    plt.rcParams['font.weight'] = 'bold'
    
    # Save results
    with open(os.path.join(results_dir, 'results.json'), 'w') as f:
        json.dump(results, f)
    
    # Save GMM parameters
    if gmm_params is not None:
        gmm_params_json = {}
        for key, value in gmm_params.items():
            if isinstance(value, np.ndarray):
                gmm_params_json[key] = value.tolist()
            elif isinstance(value, np.generic):
                gmm_params_json[key] = value.item()
            else:
                gmm_params_json[key] = value
        
        with open(os.path.join(results_dir, 'gmm_params.json'), 'w') as f:
            json.dump(gmm_params_json, f, indent=4)
        
        # Save as .mat for compatibility
        save_mat = {key: value for key, value in gmm_params_json.items() 
                   if key in ['weights', 'means', 'covariances', 'aic', 'bic', 'n_components']}
        io.savemat(os.path.join(results_dir, 'gm_params.mat'), {'gm_params': save_mat})
    
    # TSNE Visualization
    try:
        features = np.array(results['features'][f'epoch_{epoch + 1}'])
        x_embedded = TSNE(n_components=2).fit_transform(features)
        
        figure, axesSubplot = plt.subplots(figsize=(10, 8))
        scatter = axesSubplot.scatter(x_embedded[:, 0], x_embedded[:, 1],
                                     c=results['labels'][f'epoch_{epoch + 1}'],
                                     cmap='viridis', alpha=0.8)
        
        cbar = plt.colorbar(scatter)
        cbar.set_label('Cluster', fontsize=24, fontweight='bold')
        cbar.ax.tick_params(labelsize=20)
        
        title = f'TSNE Visualization (Epoch {epoch+1})'
        if gmm_params and 'aic' in gmm_params:
            title += f"\nAIC: {gmm_params['aic']:.2f}, BIC: {gmm_params['bic']:.2f}"
        axesSubplot.set_title(title, fontsize=28, fontweight='bold')
        axesSubplot.tick_params(axis='both', which='major', labelsize=20, width=2, length=6)
        
        plt.savefig(os.path.join(results_dir, 'tsne_visualization.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        # Train/Test visualization
        figure, axesSubplot = plt.subplots(figsize=(10, 8))
        train_mask = np.zeros(len(features), dtype=bool)
        train_mask[train_index] = True
        
        axesSubplot.scatter(x_embedded[train_mask, 0], x_embedded[train_mask, 1],
                          c='blue', alpha=0.6, label='Train')
        axesSubplot.scatter(x_embedded[~train_mask, 0], x_embedded[~train_mask, 1],
                          c='red', alpha=0.6, label='Test')
        
        axesSubplot.set_title('TSNE: Train vs Test', fontsize=28, fontweight='bold')
        axesSubplot.tick_params(axis='both', which='major', labelsize=20, width=2, length=6)
        axesSubplot.legend(fontsize=24, frameon=True, edgecolor='black', borderpad=1)
        
        plt.savefig(os.path.join(results_dir, 'tsne_train_test.png'), dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Error creating TSNE visualization: {e}")
    
    # Create metric plots
    metrics = [
        ('loss', 'Loss', 'train_loss', 'val_loss'),
        ('accuracy', 'Accuracy (%)', 'train_accuracy', 'val_accuracy'),
        ('ami', 'AMI Score', 'train_ami', 'val_ami')
    ]
    
    for filename, ylabel, train_key, val_key in metrics:
        plt.figure(figsize=(10, 6))
        plt.plot(results[train_key], 'b-', linewidth=3, label=f'Train {ylabel}')
        plt.plot(results[val_key], 'r-', linewidth=3, label=f'Val {ylabel}')
        plt.xlabel('Epoch', fontsize=24, fontweight='bold')
        plt.ylabel(ylabel, fontsize=24, fontweight='bold')
        plt.title(f'Training and Validation {ylabel}', fontsize=28, fontweight='bold')
        plt.legend(fontsize=24, frameon=True, edgecolor='black', borderpad=1)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.xticks(fontsize=20, fontweight='bold')
        plt.yticks(fontsize=20, fontweight='bold')
        plt.savefig(os.path.join(results_dir, f'{filename}.png'), dpi=300, bbox_inches='tight')
        plt.close()
    
    # Save indices and model
    pd.DataFrame(train_index).to_csv(os.path.join(results_dir, 'train_idx.csv'), index=False)
    pd.DataFrame(test_index).to_csv(os.path.join(results_dir, 'test_idx.csv'), index=False)
    
    torch.save({
        'epoch': epoch,
        'state_dict': model.state_dict(),
        'optimizer': optimizer.state_dict(),
        'train_index': train_index.tolist(),
        'test_index': test_index.tolist(),
    }, os.path.join(results_dir, 'model.pth.tar'))
    
    print(f"Results saved to {results_dir}")