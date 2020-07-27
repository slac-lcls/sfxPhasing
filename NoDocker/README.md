# Instruction for building a conda environment for SFX phasing on Cori

While the long-term goal is to use Shifter, we provide here a way to do without it, using a Conda environment only.

## Installation

### 1. install conda
```bash
cd sfxPhasing/Docker/conda.local
. sites/default.sh
./install_conda.sh
chmod +x env.sh
chmod +x env.local
sync
./env.sh
```

### 2. create the sfx conda environment
```bash
cd sfxPhasing
source Docker/conda.local/env.sh
conda create --name sfx --file requirements.txt
conda env list
```

### 3. activate the sfx conda environment
```bash
source activate sfx
conda env list
```

### 4. install crystallography libraries

Work in progres...

## Usage
```bash
source /.../sfxPhasing/Docker/conda.local/env.local
source activate sfx
source /.../ccp4-7.1/bin/ccp4.setup-sh
```
