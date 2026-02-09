# Tau Heterogenity

## Prerequisites

* Access to a system with **Python 3.11.4** available as a module
* CUDA 12.4–compatible GPU (for PyTorch with CUDA support)

---

## Steps Required

### 1. Load Python Module

```bash
module load python/3.11.4
```

### 2. Create a Virtual Environment

```bash
python3.11 -m venv venv_tau
```

### 3. Activate the Virtual Environment

```bash
source venv_tau/bin/activate
```

### 4. Install PyTorch Dependencies (CUDA 12.4)

```bash
pip install torchvision==0.20.0+cu124 torchaudio==2.5.0+cu124 \
  --index-url https://download.pytorch.org/whl/cu124
```

### 5. Install Remaining Requirements

```bash
pip install -r requirements.txt
```

---

## Training the Model

Run the following command to start model training:

```bash
python main_scripts.py \
  --data_path /N/slate/thjaya/Final/pet.csv \
  --script_path /N/slate/thjaya/Final/TheyaneshProjectFinal_LastADNI_1/src/main.py \
  --output_dir /N/slate/thjaya/Final1/ADNIResultsFinal292026 \
  --outer_folds 10 \
  --inner_folds 2 \
  --max_workers 4
```

### Argument Description

* `--data_path`: Path to the PET CSV dataset
* `--script_path`: Path to the main training script
* `--output_dir`: Directory where results will be saved
* `--outer_folds`: Number of outer cross-validation folds
* `--inner_folds`: Number of inner cross-validation folds
* `--max_workers`: Number of parallel worker processes

---

## Model Evaluation / Testing

After training completes, select the **best-performing fold** based on:

* AIC
* BIC
* Validation classification loss

Once the best fold is identified, use the corresponding trained model and GMM parameters to run inference on the dataset.

### Run Model in Test Mode

```bash
python main.py \
  --mode test \
  --model /N/slate/thjaya/Final1/ADNIResultsFinal292026/Results_train_outer_fold=0_type=0_k=3_lr=0.05_reassign=1.0/model.pth.tar \
  --gmmmodel /N/slate/thjaya/Final1/ADNIResultsFinal292026/Results_train_outer_fold=0_type=0_k=3_lr=0.05_reassign=1.0/gm_params.mat \
  --data /N/slate/thjaya/Final/pet.csv \
  --output_dir /N/slate/thjaya/Final1/ADNIResultsFinal292026 \
  --batch 32 \
  --workers 1
```

### Explanation

* The model and GMM parameters are loaded from the **best training fold**.
* The full dataset is passed again for evaluation/inference.
* This step produces final predictions and evaluation outputs using the selected optimal model.

---

## Citation

This codebase has been adapted and modified for Tau-PET data.

The original methodology has been extended and customized to support Tau-PET–based analysis, training, and evaluation. If you use or adapt this code, please cite the original work:

---

```bibtex
@article{KANG2024120737,
  title   = {Disentangling brain atrophy heterogeneity in Alzheimer's disease: A deep self-supervised approach with interpretable latent space},
  journal = {NeuroImage},
  volume  = {297},
  pages   = {120737},
  year    = {2024},
  issn    = {1053-8119},
  doi     = {https://doi.org/10.1016/j.neuroimage.2024.120737},
  url     = {https://www.sciencedirect.com/science/article/pii/S1053811924002301},
  author  = {Sohyun Kang and Sung-Woo Kim and Joon-Kyung Seong}
}
```

---

## Notes

* Ensure the virtual environment is activated before running any commands.
* Paths should be updated according to your local or cluster directory structure.
