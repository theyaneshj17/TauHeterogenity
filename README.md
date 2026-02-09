# TauHeterogenity


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

## Notes

* Ensure the virtual environment is activated before running any commands.
* Paths should be updated according to your local or cluster directory structure.
