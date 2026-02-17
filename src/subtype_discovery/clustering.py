import time
from sklearn import mixture
import numpy as np


class ReassignedDataset:
    """Dataset with reassigned pseudo-labels."""
    
    def __init__(self, dataset, pseudolabels):
        self.data = list(zip(dataset, pseudolabels))
        self.len = len(self.data)

    def __getitem__(self, index):
        img, pseudolabel = self.data[index]
        return [img, pseudolabel]

    def __len__(self):
        return self.len


def calculate_dynamic_reg_covar(features, method='adaptive', factor=0.1, min_reg=1e-6):
    """Calculate dynamic regularization parameter based on feature characteristics."""
    feature_vars = np.var(features, axis=0)
    mean_var = np.mean(feature_vars)
    min_var = np.min(feature_vars)
    median_var = np.median(feature_vars)
    
    if min_var < 1e-10:
        print(f"Warning: Very low minimum variance detected: {min_var}")
        min_var = 1e-10
    
    var_ratio = np.max(feature_vars) / min_var
    
    if method == 'adaptive':
        if var_ratio > 100:
            reg_value = max(0.05 * mean_var, min_reg)
        elif var_ratio > 10:
            reg_value = max(0.01 * mean_var, min_reg)
        else:
            reg_value = max(0.005 * mean_var, min_reg)
    elif method == 'min_based':
        reg_value = max(factor * min_var, min_reg)
    elif method == 'median_based':
        reg_value = max(factor * median_var, min_reg)
    else:
        reg_value = max(0.01 * mean_var, min_reg)
    
    return reg_value


def run_gmm(x, train_idx, val_idx, k, dynamic_reg=True):
    """Run Gaussian Mixture Model clustering and calculate information criteria."""
    train_features = x[train_idx]
    
    # Determine regularization parameter
    if dynamic_reg:
        reg_covar = calculate_dynamic_reg_covar(train_features, method='adaptive')
        print(f"Using dynamic regularization: {reg_covar:.8f}")
    else:
        reg_covar = 1e-2
        print(f"Using fixed regularization: {reg_covar}")
    
    # Initialize GMM
    clf = mixture.GaussianMixture(
        n_components=k,
        covariance_type='full',
        reg_covar=reg_covar,
        random_state=31,
        n_init=10,
        max_iter=100
    )
    
    label_array = np.ones(len(x))
    
    # Fit GMM on training data only
    clf.fit(train_features)
    
    # Calculate log-likelihood on training and validation sets
    train_log_likelihood = clf.score(train_features) * len(train_idx)
    
    # Predict on both sets
    label_train = clf.predict(x[train_idx])
    label_test = clf.predict(x[val_idx])
    
    # Calculate validation log-likelihood
    val_features = x[val_idx]
    val_log_likelihood = clf.score(val_features) * len(val_idx)
    
    # Calculate number of free parameters
    features_dim = x.shape[1]
    n_parameters = (k * features_dim +  # Means
                   k * features_dim * (features_dim + 1) / 2 +  # Covariances
                   (k - 1))  # Weights
    
    # Calculate AIC/BIC on validation set (for model selection)
    n_samples_val = len(val_idx)
    aic_val = 2 * n_parameters - 2 * val_log_likelihood
    bic_val = np.log(n_samples_val) * n_parameters - 2 * val_log_likelihood
    
    # Calculate AIC/BIC on training set 
    n_samples_train = len(train_idx)
    aic_train = 2 * n_parameters - 2 * train_log_likelihood
    bic_train = np.log(n_samples_train) * n_parameters - 2 * train_log_likelihood
    
    # Store GMM parameters
    gmm_params = {
        'weights': clf.weights_.tolist(),
        'means': clf.means_.tolist(),
        'covariances': clf.covariances_.tolist(),
        'precisions': clf.precisions_.tolist(),
        'converged': bool(clf.converged_),
        'n_iter': int(clf.n_iter_),
        'n_components': int(k),
        'covariance_type': 'full',
        'reg_covar': float(reg_covar),
        
        # Log-likelihoods
        'train_log_likelihood': float(train_log_likelihood),
        'val_log_likelihood': float(val_log_likelihood),
        
        # Information criteria (validation set - for model selection)
        'aic': float(aic_val),
        'bic': float(bic_val),
        

        'train_aic': float(aic_train),
        'train_bic': float(bic_train),
        'n_parameters': int(n_parameters),
        'n_samples_train': int(n_samples_train),
        'n_samples_val': int(n_samples_val),
        'features_dim': int(features_dim)
    }
    
    # Combine labels
    label_array[train_idx] = label_train
    label_array[val_idx] = label_test
    
    return label_array, gmm_params


class GMM:
    """Gaussian Mixture Model clustering wrapper."""
    
    def __init__(self, k, dynamic_reg=True):
        self.k = k
        self.gmm = None
        self.aic = None
        self.bic = None
        self.images_lists = None
        self.gmm_params = None
        self.dynamic_reg = dynamic_reg
    
    def cluster(self, data, train_idx, val_idx):
        """Perform GMM-based clustering."""
        label_array, self.gmm_params = run_gmm(data, train_idx, val_idx, self.k, self.dynamic_reg)
        
        # Store AIC/BIC
        if self.gmm_params:
            self.aic = self.gmm_params['aic']
            self.bic = self.gmm_params['bic']
        
        # Create image lists
        self.images_lists = [[] for _ in range(self.k)]
        for i in range(len(data)):
            self.images_lists[int(label_array[i])].append(i)
        
        return label_array, self.gmm_params
    
    def get_information_criteria(self):
        """Get calculated information criteria."""
        if self.gmm_params is None:
            raise ValueError("Model must be fitted using cluster() before accessing information criteria")
            
        return {
            'aic': self.aic,
            'bic': self.bic,
            'gmm_params': self.gmm_params
        }