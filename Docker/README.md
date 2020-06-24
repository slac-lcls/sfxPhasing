# Instructions for building a shifter image for sfx workflow on cori


## Getting Started
These instructions will let you build your own shifter image for sfx workflow and push it on spin registry so it can be deployed on one of the cori nodes

### Prerequisites

phenix version: 1.16 Freely available after requesting a download password: https://www.phenix-online.org/download/

ccp4 version: 7.0.077 (It includes SHELXC version:2016/1, SHELXD version:2013/2, and Crank2 version:2.0.229) Freely available: http://www.ccp4.ac.uk/download

datasets: download from here  https://stanford.box.com/s/8s55n5smgxb7gt00ihzoj241qgk1bocl

nersc-compatible conda installation: included in conda.local/

docker running on your machine.

In order for the docker build to work you need all downloads (phenix,ccp4,datasets, sfx code) in the same directory as the Dockerfile.base

```
```

### Building the docker image  and uploading it to nersc shifter registry

The Dockerfile will handle the installation of all software components. For ccp4 a custom path is required under 
```
CCP4_SCR="/tmp/`whoami | tr ' \\\\' _`" 
```
the Dockerfile will automatically create that path under 
```
img/tmp/your_nersc_username
```
 make necessary changes in the ccp4 config file. To build the image simply run:
```
docker build -f Dockerfile.base --build-arg nersc_user=your_nersc_username -t sfx .
```
Once building is complete:
Get your image id:
```
docker image list | grep sfx
```
Tag your image:
```
docker tag image_id registry.services.nersc.gov/your_nersc_username/sfx
```
Push image on spin registry (you need an authentication token for that). Instructions on how to obtain the token can be found here: https://docs.nersc.gov/development/shifter/how-to-use/

```
docker push registry.services.nersc.gov/your_nersc_username/sfx
``` 

## Deploying the image interactively on cori:
Login to cori and pull the image from spin registry:
```
shifterimg pull registry.services.nersc.gov/your_nersc_username/sfx
```
request an interactive node while passing the image as a parameter
```
salloc -N 1 -q interactive -A your_project_account -C hardware_node_type --image=docker:registry.services.nersc.gov/your_nersc_username/sfx -t 02:00:00
``` 
login to shifter image:
```
shifter /bin/bash
```
activate sfx conda env:
```
source /img/conda.local/env.local
source acticate sfx
```
load ccp4 and phenix params:
```
source /img/phenix/phenix-1.18.2-3874/phenix_env.sh
source /img/ccp4/ccp4-7.0/bin/ccp4.setup-sh
```
