

import argparse
import time
import random
import os
import json

import network
import clustering
import petdataset
from util import label_stab, UnifLabelSampler, AverageMeter, plotandsave, load_model, load_gmm

from sklearn.metrics import adjusted_mutual_info_score
from sklearn.model_selection import KFold
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.backends.cudnn as cudnn

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")




print('Torch:', torch.__version__); 
print('CUDA:', torch.version.cuda, torch.cuda.is_available())



def parse_args():
    parser = argparse.ArgumentParser(description='PyTorch Implementation of DebaNet')
    parser.add_argument('--mode', type=str, choices=['train', 'test'], default='train',
                        help='Mode: train or test')
    
    # Common arguments
    parser.add_argument('--data', type=str, default='ADNI', help='Data type to use')
    parser.add_argument('--batch', type=int, default=32, help='Mini-batch size')
    parser.add_argument('--workers', type=int, default=4, help='Number of data loading workers (use 0 for Windows)')
    parser.add_argument('--seed', type=int, default=31, help='Random seed')
    parser.add_argument('--output_dir', type=str, default='./Results', help='Directory to store output files')
    
    # Training-specific arguments
    parser.add_argument('--num_cluster', '--k', type=int, default=3, help='Number of clusters')
    parser.add_argument('--lr', type=float, default=0.001, help='Learning rate')
    parser.add_argument('--reassign', type=float, default=1.0, help='Epochs between cluster reassignments')
    parser.add_argument('--wd', type=float, default=-4, help='Weight decay power')
    parser.add_argument('--epochs', type=int, default=1000, help='Total number of epochs')
    parser.add_argument('--fold_seed', type=int, default=31, help='Random seed for cross-validation')
    parser.add_argument('--fold_idx', type=int, default=0, help='Fold index for cross-validation')
    parser.add_argument('--train_indices', type=str, default='', help='Comma-separated training indices')
    parser.add_argument('--val_indices', type=str, default='', help='Comma-separated validation indices')
    parser.add_argument('--t', type=int, default=0, help='Inner/Outer fold indicator')
    parser.add_argument('--outer_fold', type=int, default=-1, help='Outer fold number')
    parser.add_argument('--patience', type=int, default=10, help='Patience for early stopping')
    parser.add_argument('--min_lr', type=float, default=1e-5, help='Minimum learning rate')
    
    # Test-specific arguments
    parser.add_argument('--model', type=str, default='', help='Path to trained model checkpoint')
    parser.add_argument('--gmmmodel', type=str, default='', help='Path to GMM parameters')
    
    return parser.parse_args()


def validate_indices(train_indices, val_indices=None, test_indices=None):
    """Validate no overlap between train, validation and test indices."""
    train_set = set(train_indices)
    
    if val_indices is not None:
        val_set = set(val_indices)
        if len(train_set.intersection(val_set)) > 0:
            raise ValueError("Data leakage: overlap between train and validation sets")
    
    if test_indices is not None:
        test_set = set(test_indices)
        if len(train_set.intersection(test_set)) > 0:
            raise ValueError("Data leakage: overlap between train and test sets")
        if val_indices is not None and len(val_set.intersection(test_set)) > 0:
            raise ValueError("Data leakage: overlap between validation and test sets")
    
    return True


def compute_features(dataloader, model, dataset_size, indices=None):
    """Extract features from the model for specified data samples."""
    model.eval()
    all_features = []
    
    with torch.no_grad():
        for input_tensor in dataloader:
            batch_features = model(input_tensor).cpu().numpy()
            all_features.append(batch_features)
    
    all_features = np.vstack(all_features)
    
    if indices is not None:
        valid_indices = [idx for idx in indices if idx < all_features.shape[0]]
        
        if len(valid_indices) < len(indices):
            print(f"WARNING: {len(indices) - len(valid_indices)} indices out of bounds")
        
        if not valid_indices:
            raise ValueError("No valid indices remain!")
            
        return all_features[valid_indices], np.array(valid_indices)
    
    return all_features, np.arange(all_features.shape[0])


def train_epoch(loader, model, crit, opt, epoch, total_epochs, args):
    """Train the model for one epoch."""
    batch_count = 0
    model.train()
    optimizer_tl = torch.optim.AdamW(
        model.top_layer.parameters(),
        lr=args.lr,
        weight_decay=10**args.wd,
    )
    
    total_loss = 0
    for input_tensor, target in loader:


                # Forward pass
        if model.decoding:
            classification, decoded = model.get_training_outputs(input_tensor)
            recon_loss = nn.functional.mse_loss(decoded, input_tensor)
            class_loss = crit(classification, target)
            recon_weight = max(0.1 * (1.0 - epoch/total_epochs), 0.01)
            loss = class_loss + recon_weight * recon_loss + model.get_l2_regularization()

 


        else:
            classification = model(input_tensor)
            loss = crit(classification, target) + model.get_l2_regularization()

        # Backward pass
        opt.zero_grad()
        optimizer_tl.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.3)
        opt.step()
        optimizer_tl.step()
        total_loss += loss.item()




        batch_count += 1
        
    return total_loss / len(loader)


def feedforward(loader, model, crit):
    """Evaluate model performance."""
    losses = AverageMeter()
    model.eval()
    total = 0
    correct = 0
    outputs = None
    
    with torch.no_grad():
        for i, (input_tensor, target) in enumerate(loader):
            model_output = model(input_tensor)
            
            if model.decoding:
                classification_output = model_output[0] if isinstance(model_output, tuple) else model_output
            else:
                classification_output = model_output
                
            if len(classification_output.shape) == 1:
                classification_output = classification_output.unsqueeze(0)
                
            if i == 0:
                outputs = classification_output
            else:
                outputs = torch.cat((outputs, classification_output), dim=0)
                
            if len(target.shape) == 1:
                target = target.view(-1)
                
            loss = crit(classification_output, target)
            _, predicted = torch.max(classification_output, 1)
            total += target.size(0)
            correct += (predicted == target).sum().item()
            losses.update(loss.item(), input_tensor.size(0))
        
    return losses.avg, 100*correct/total, outputs.cpu().numpy()


def train_mode(args):
    """Training mode."""
    print("\n" + "="*60)
    print("TRAINING MODE")
    print("="*60 + "\n")
    
    # Set random seeds
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(args.seed)
    np.random.seed(args.seed)
    random.seed(args.seed)

    print(f'K: {args.num_cluster}')
    print(f'Learning rate: {args.lr}')
    print(f'Reassign: {args.reassign}')

    NoofVertices = 68

    # Initialize model
    model = network.improved_mlp(
        input_dim=NoofVertices,
        output_dim=args.num_cluster,
        dropout_rate=0.4,
        l2_lambda=1e-2,
        decoding=True
    )

    #print(model)

    fd = int(model.top_layer.weight.size()[1])
    model.top_layer = None
    
    cudnn.benchmark = True

    # Setup optimizer
    optimizer = torch.optim.AdamW(
        filter(lambda x: x.requires_grad, model.parameters()),
        lr=args.lr,
        weight_decay=10**args.wd,
    )

    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=0.2,
        patience=10,
        min_lr=args.min_lr,
        verbose= True
    )

    criterion = nn.CrossEntropyLoss()

    # Load dataset
    print("Loading dataset...")
    dataset = petdataset.PetDataset(datatype=args.data)
    print(f"Total samples: {len(dataset)}")
    
    dataloader = torch.utils.data.DataLoader(
        dataset,
        batch_size=args.batch,
        num_workers=args.workers,
        pin_memory=True
    )

    # Train/validation split
    if args.train_indices and args.val_indices:
        train_index = np.array([int(x) for x in args.train_indices.split(',')])
        test_index = np.array([int(x) for x in args.val_indices.split(',')])
    else:
        skf = KFold(n_splits=10, shuffle=True, random_state=args.fold_seed)
        all_indices = np.arange(len(dataset))
        idx_list = list(skf.split(all_indices))
        train_index, test_index = idx_list[args.fold_idx]
        train_index = all_indices[train_index]
        test_index = all_indices[test_index]

    validate_indices(train_index, test_indices=test_index)

    # Initialize GMM
    GMM = clustering.GMM(args.num_cluster, dynamic_reg=True)

    # Results tracking
    results = {
        'train_loss': [],
        'train_accuracy': [],
        'train_ami': [],
        'val_loss': [],
        'val_accuracy': [],
        'val_ami': [],
        'features': {},
        'labels': {},
        'val_predictions': {},
    }

    # Training loop
    for epoch in range(args.epochs):
        end = time.time()
        
        # Feature extraction
        model.top_layer = None
        model.features = nn.Sequential(*list(model.features.children())[:-1])
        
        train_features, valid_train_indices = compute_features(dataloader, model, len(dataset), train_index)
        
        # Log feature variance on training data only
        feature_variance = np.var(train_features, axis=0)
        #print(f'Variance of Encoded Features for Epoch {epoch + 1}:', feature_variance)
        high_variance_features = np.sum(feature_variance > 0.1)
        #print(f'Number of High Variance Features (threshold 0.1): {high_variance_features}')
        
        # Now compute features for test set, getting valid indices
        test_features, valid_test_indices = compute_features(dataloader, model, len(dataset), test_index)
        
        features = np.zeros((len(dataset), train_features.shape[1]), dtype='float32')
        for i, idx in enumerate(valid_train_indices):
            features[idx] = train_features[i]
        for i, idx in enumerate(valid_test_indices):
            features[idx] = test_features[i]
        
        train_index = valid_train_indices
        test_index = valid_test_indices
        
        # GMM clustering
        label_unmatched, gmm_params = GMM.cluster(features, train_index, test_index)
        images_lists_unmatched = GMM.images_lists

        # Label stabilization
        if epoch == 0:
            label_prev = label_unmatched.astype(int)
            images_lists_prev = images_lists_unmatched
        else:
            label_prev = label_now
            images_lists_prev = images_lists_now

        label_now, images_lists_now = label_stab(images_lists_unmatched, images_lists_prev, label_unmatched)


        results['features'][f'epoch_{epoch + 1}'] = features.tolist()
        results['labels'][f'epoch_{epoch+1}'] = label_now.tolist()

        # Compute AMI
        train_ami = adjusted_mutual_info_score(label_prev[train_index], label_now[train_index])
        val_ami = adjusted_mutual_info_score(label_prev[test_index], label_now[test_index])
        results['train_ami'].append(train_ami)
        results['val_ami'].append(val_ami)

        # Create datasets
        train_dataset = clustering.ReassignedDataset(dataset.x_data.numpy()[train_index], label_now[train_index])
        test_dataset = clustering.ReassignedDataset(dataset.x_data.numpy()[test_index], label_now[test_index])

        # Create samplers
        images_lists_train_abs = [[] for _ in range(len(images_lists_now))]
        images_lists_test_abs = [[] for _ in range(len(images_lists_now))]


        train_indices_set = set(train_index)
        test_indices_set = set(test_index)
        
        # Ensure no duplicates in train/test indices but keep original ordering
        train_index_unique = np.array([idx for idx in train_index if idx in train_indices_set], dtype=int)
        test_index_unique = np.array([idx for idx in test_index if idx in test_indices_set], dtype=int)
        
        # Update train_index and test_index for later use, preserving the original order
        train_index = train_index_unique
        test_index = test_index_unique
        
        for label, image_list in enumerate(images_lists_now):
            images_lists_train_abs[label] = [i for i, idx in enumerate(train_index) if idx in image_list]
            images_lists_test_abs[label] = [i for i, idx in enumerate(test_index) if idx in image_list]





        train_sampler = UnifLabelSampler(int(args.reassign * len(train_dataset)), images_lists_train_abs)
        test_sampler = UnifLabelSampler(int(args.reassign * len(test_dataset)), images_lists_test_abs)



       
        train_dataloader = torch.utils.data.DataLoader(
            train_dataset,
            batch_size=args.batch,
            num_workers=args.workers,
            sampler=train_sampler,
            pin_memory=True
        )

        test_dataloader = torch.utils.data.DataLoader(
            test_dataset,
            batch_size=args.batch,
            num_workers=args.workers,
            pin_memory=True
        )

        # Prepare for training
        tmp = list(model.features.children())
        tmp.append(nn.ReLU())
        model.features = nn.Sequential(*tmp)
        
        if epoch == 0:
            model.top_layer = nn.Linear(fd, len(images_lists_now))
            model.top_layer.weight.data.normal_(0, 0.01)
            model.top_layer.bias.data.zero_()
        else:
            model.top_layer = top_layer_tmp

        # Train
        loss_t = train_epoch(train_dataloader, model, criterion, optimizer, epoch, args.epochs, args)

        top_layer_tmp = model.top_layer
        
        # Evaluate
        train_dataloader_wo_sampling = torch.utils.data.DataLoader(
            train_dataset,
            batch_size=args.batch,
            num_workers=args.workers,
            pin_memory=True
        )

        loss_t, acc_t, output_train = feedforward(train_dataloader_wo_sampling, model, criterion)
        loss_v, acc_v, output_test = feedforward(test_dataloader, model, criterion)

        results['train_loss'].append(loss_t)
        results['train_accuracy'].append(acc_t)
        results['val_loss'].append(loss_v)
        results['val_accuracy'].append(acc_v)
        results['val_predictions'][f'epoch_{epoch + 1}'] = output_test.tolist()
        
        print(f'Epoch [{epoch}] | Time: {time.time() - end:.2f}s | '
              f'Train loss: {loss_t:.2f} | Val loss: {loss_v:.2f} | '
              f'Train acc: {acc_t:.1f} | Val acc: {acc_v:.1f} | '
              f'Train AMI: {train_ami:.1f} | Val AMI: {val_ami:.1f} | Inclass size: {[len(x) for x in images_lists_now]}')


        scheduler.step(loss_v)
        
        # Dynamic L2
        if epoch > 0 and loss_v > results['val_loss'][-1]:
            model.l2_lambda *= 1.2
        else:
            model.l2_lambda *= 0.8
        model.l2_lambda = max(1e-6, min(1e-2, model.l2_lambda))

    # Save results
    if args.t == 1:
        profile = f"train_outer_fold={args.outer_fold}_inner_fold={args.fold_idx}_type={args.t}_k={args.num_cluster}_lr={args.lr}_reassign={args.reassign}"
    else:
        profile = f"train_outer_fold={args.outer_fold}_type={args.t}_k={args.num_cluster}_lr={args.lr}_reassign={args.reassign}"

    results_dir = os.path.join(args.output_dir, f'Results_{profile}')
    plotandsave(results_dir, results, train_index, test_index, gmm_params, model, optimizer, epoch)
    
    print(f"\nTraining complete! Results saved to: {results_dir}")


def test_mode(args):
    """Testing mode."""
    print("\n" + "="*60)
    print("TESTING MODE")
    print("="*60 + "\n")
    
    # Set random seeds
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(args.seed)
    np.random.seed(args.seed)
    random.seed(args.seed)
    
    # Load model
    print(f"Loading model from: {args.model}")
    model = load_model(args.model)
    if model is None:
        print(f"ERROR: Failed to load model")
        exit(1)
    print("Model loaded successfully!")
    
    # Load GMM
    print(f"Loading GMM from: {args.gmmmodel}")
    gmm = load_gmm(args.gmmmodel)
    if gmm is None:
        print(f"ERROR: Failed to load GMM")
        exit(1)
    print(f"GMM loaded successfully! ({gmm.n_components} components)")
    
    cudnn.benchmark = True
    
    # Load dataset
    print("\nLoading dataset...")
    dataset = petdataset.PetDataset(datatype=args.data)
    print(f"Total samples: {len(dataset)}")
    
    dataloader = torch.utils.data.DataLoader(
        dataset,
        batch_size=args.batch,
        num_workers=args.workers,
        pin_memory=False
    )
    
    # Extract features
    print("\n" + "-"*60)
    print("Extracting features...")
    print("-"*60)
    
    top_layer = model.top_layer
    model.top_layer = None
    model.features = nn.Sequential(*list(model.features.children())[:-1])
    
    features, valid_indices = compute_features(dataloader, model, len(dataset))
    print(f"Features: {features.shape}")
    
    # Apply GMM
    print("\n" + "-"*60)
    print("Applying GMM clustering...")
    print("-"*60)
    
    label = gmm.predict(features)
    unique_labels, label_counts = np.unique(label, return_counts=True)
    print(f"GMM clusters: {dict(zip(unique_labels, label_counts))}")
    
    # Get model predictions
    print("\n" + "-"*60)
    print("Getting model predictions...")
    print("-"*60)
    
    tmp = list(model.features.children())
    tmp.append(nn.ReLU())
    model.features = nn.Sequential(*tmp)
    model.top_layer = top_layer
    
    test_dataset = clustering.ReassignedDataset(dataset.x_data.numpy()[valid_indices], label)
    test_dataloader = torch.utils.data.DataLoader(
        test_dataset,
        batch_size=args.batch,
        num_workers=args.workers,
        pin_memory=False
    )
    
    criterion = nn.CrossEntropyLoss()
    _, _, outputs = feedforward(test_dataloader, model, criterion)
    predicted_labels = np.argmax(outputs, axis=1)
    
    unique_preds, pred_counts = np.unique(predicted_labels, return_counts=True)
    print(f"Model predictions: {dict(zip(unique_preds, pred_counts))}")
    
    # Label stabilization
    print("\n" + "-"*60)
    print("Performing label stabilization...")
    print("-"*60)
    
    n_clusters = len(unique_labels)
    image_list_unmatched = [[] for _ in range(n_clusters)]
    for i, lbl in enumerate(label):
        if int(lbl) < n_clusters:
            image_list_unmatched[int(lbl)].append(i)
    
    image_list_predicted = [[] for _ in range(n_clusters)]
    for i, pred_lbl in enumerate(predicted_labels):
        if int(pred_lbl) < n_clusters:
            image_list_predicted[int(pred_lbl)].append(i)
    
    stabilized_label, _ = label_stab(image_list_unmatched, image_list_predicted, label)
    unique_stable, stable_counts = np.unique(stabilized_label, return_counts=True)
    print(f"Stabilized labels: {dict(zip(unique_stable, stable_counts))}")
    
    # Final evaluation
    print("\n" + "-"*60)
    print("Final evaluation...")
    print("-"*60)
    
    test_dataset_final = clustering.ReassignedDataset(
        dataset.x_data.numpy()[valid_indices],
        stabilized_label
    )
    test_dataloader_final = torch.utils.data.DataLoader(
        test_dataset_final,
        batch_size=args.batch,
        num_workers=args.workers,
        pin_memory=False
    )
    
    loss, acc, final_outputs = feedforward(test_dataloader_final, model, criterion)
    
    print("\n" + "="*60)
    print("FINAL RESULTS")
    print("="*60)
    print(f"Loss: {loss:.4f}")
    print(f"Accuracy: {acc:.2f}%")
    print(f"Clusters: {dict(zip(unique_stable, stable_counts))}")
    
    # Save results
    print("\n" + "-"*60)
    print("Saving results...")
    print("-"*60)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    subject_ids = dataset.subject_ids.values[valid_indices]
    
    results_df = pd.DataFrame({
        'SubjectID': subject_ids,
        'GMM_Label': label,
        'Model_Prediction': predicted_labels,
        'Stabilized_Label': stabilized_label
    })
    
    results_df.to_csv(os.path.join(args.output_dir, 'predicted_labels.csv'), index=False)
    pd.DataFrame(final_outputs).to_csv(os.path.join(args.output_dir, 'model_outputs.csv'), index=False)
    
    metrics = {
        'loss': float(loss),
        'accuracy': float(acc),
        'num_samples': int(len(valid_indices)),
        'num_clusters': int(gmm.n_components),
        'cluster_distribution': {int(k): int(v) for k, v in zip(unique_stable, stable_counts)}
    }
    
    with open(os.path.join(args.output_dir, 'test_metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=4)
    
    print(f"\nResults saved to: {args.output_dir}")
    print("\n" + "="*60)
    print("TESTING COMPLETED")
    print("="*60 + "\n")


def main():
    args = parse_args()
    
    if args.mode == 'train':
        train_mode(args)
    elif args.mode == 'test':
        test_mode(args)
    else:
        raise ValueError(f"Invalid mode: {args.mode}")


if __name__ == '__main__':
    main()