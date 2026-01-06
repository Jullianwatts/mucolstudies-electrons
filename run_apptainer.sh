#apptainer build k4toroid.sif docker://madbaron/k4test-ubuntu:latest
apptainer run --no-home -B /ospool/uc-shared/project/futurecolliders/data:/data -B /scratch/$USER -B /home/$USER /scratch/trholmes/mucol/k4toroid.sif
