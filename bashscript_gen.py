import os

#######################
### User parameters ###
#######################

RUN = False  # Set to True to submit the job immediately

# Name of the PLAMS script to be run
script = "gwrose_script_chlorophyl.py"

# Define the molecules, methods, and dimer separations
#molecules = ["ethylene"]#,"pyrene", "pyridine"]
#molecules = ["ethylene"]
molecules = ["chlorophyl"]
methods = ["qsGW"]#,"qsGW","TDDFT"]
#distances = [3.0, 4.0, 5.0, 6.0, 8.0, 10.0]  # Corresponds to dimer separations
#distances = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 , 12.0 , 13.0 , 14.0]
distances = [3.0]

# Output directory
cluster_info = "/scistor/tc/huw587" #Input of cluster 

output_name = "chl_dim_qsgw_occvirt"
outdir = f"{cluster_info}/RoseADF_Notebooks/"
outdirmain = f"{cluster_info}/{output_name}"
outdir1 = f"{outdirmain}/workdirs"
outdir_coupling = f"{outdirmain}/coupling_mat"


# SLURM Job parameters
njobs = 2               # Max concurrent jobs
time = "10-00:00:00"      # Walltime
mem = ""             # Memory allocation
partition = "tc"         # Partition namex

# Compute the total number of jobs
num_jobs = len(molecules) * len(methods) * len(distances)

# Additional setup commands
additional_lines = [
    'module purge\n',
    'module load shared\n',
    'module load 2024\n',
    'module load slurm\n',
    'module load intel/mkl/64/2020/4.304\n',
    'module load Python/3.11.5-GCCcore-13.2.0\n', 
    'export SCM_TMPDIR="$TMPDIR"\n',
    'export AMSHOME="/scistor/tc/afr790/amshome"\n',
    'export AMSBIN="/scistor/tc/afr790/amshome/bin"\n',
    'export AMSRESOURCES="$AMSHOME/atomicdata"\n',
    'export SCMLICENSE="/scistor/tc/huw587/ADF_trunk/license.txt"\n',
    f'mkdir -pv {outdir}\n'
]

# Create directories
current_dir = os.getcwd()
job_dir = os.path.join(current_dir, 'jobs')
out_dir = os.path.join(current_dir, 'out')
err_dir = os.path.join(current_dir, 'err')

os.makedirs(job_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)
os.makedirs(err_dir, exist_ok=True)

# Generate SLURM job script
job_file = os.path.join(job_dir, f"{output_name}.job")

with open(job_file, 'wt') as fh:
    fh.writelines(
        [
            "#!/bin/bash\n",
            f"#SBATCH --job-name=multi_calc\n",
            f"#SBATCH --output={os.path.join(out_dir, 'multi_calc_%A_%a.out')}\n",
            f"#SBATCH --error={os.path.join(err_dir, 'multi_calc_%A_%a.err')}\n",
            f"#SBATCH --time={time}\n",
           # f"#SBATCH --mem={mem}\n",
            f"#SBATCH --partition={partition}\n",
            f"#SBATCH --ntasks-per-node=32\n",
            f"#SBATCH --ntasks-per-core=1\n",
            f"#SBATCH -N 1\n",
            f"#SBATCH --mail-type=NONE\n",
            f"#SBATCH --array=0-{num_jobs - 1}%{njobs}\n",
        ] + additional_lines + [
            "echo $TMPDIR\n",
            "unset SLURM_JOB_ID\n",
            "cd $TMPDIR\n",
            f"cp -rfv {outdir}/{script} {outdir}/geometry $TMPDIR\n",
            "\n"
            "# Compute molecule, method, and separation index\n",
            f"molecule_index=$((SLURM_ARRAY_TASK_ID / ({len(methods)} * {len(distances)})))\n",
            f"method_index=$(( (SLURM_ARRAY_TASK_ID / {len(distances)}) % {len(methods)} ))\n",
            f"dist_index=$((SLURM_ARRAY_TASK_ID % {len(distances)}))\n",
            f"molecules=({' '.join(molecules)})\n",
            f"methods=({' '.join(methods)})\n",
            f"distances=({' '.join(map(str, distances))})\n",
            "molecule=${molecules[$molecule_index]}\n",
            "method=${methods[$method_index]}\n",
            "distance=${distances[$dist_index]}\n",
            "\n",
            f"$AMSBIN/amspython -u {script} $molecule $distance $method\n",
            f"mv plams_workdir plams_workdir_${{molecule}}_${{method}}_${{distance}}\n",
            f"mkdir -pv {outdir1}/$molecule/$method\n",
            f"mkdir -pv {outdir_coupling}/$molecule/$method\n",
            f"cp -rfv plams_workdir_${{molecule}}_${{method}}_${{distance}} {outdir1}/$molecule/$method\n",
            f"cp -rfv output* {outdir_coupling}/$molecule/$method\n"
        ]
    )

if RUN:
    os.system("sbatch " + job_file)


