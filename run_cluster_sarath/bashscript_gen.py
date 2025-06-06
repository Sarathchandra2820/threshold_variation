import os


#######################
### User parameters ###
#######################

RUN = False  # Set to True to submit the job immediately

# Name of the PLAMS script to be run
script = "tv_benzene.py"
start_av,end_av,step_av = 2,6,12
start_ft,end_ft,step_ft = 8,15,12

def linspace_list(start, stop, num):
    if num == 1:
        return [start]
    step = (stop - start) / (num - 1)
    return [round(start + i * step, 2) for i in range(num)]

molecules = ["bezene"]
add_virt_threshold = linspace_list(start_av,end_av,step_av)
frag_threshold = linspace_list(start_ft,end_ft,step_ft)

print(add_virt_threshold)
print(frag_threshold)

# # # Output directory
# cluster_info = "/scistor/tc/ehc732" #Input of cluster  threshold_variation 

# output_name = "threshold_output"
# outdir = f"{cluster_info}/threshold_variation/"
# outdirmain = f"{cluster_info}/{output_name}"
# outdir1 = f"{outdirmain}/workdirs"
# outdir_coupling = f"{outdirmain}/coupling_mat"


# # SLURM Job parameters
# njobs = 2               # Max concurrent jobs
# time = "10-00:00:00"      # Walltime
# mem = ""             # Memory allocation
# partition = "tc"         # Partition namex

# # Compute the total number of jobs
# num_jobs = len(molecules) * len(add_virt_threshold) * len(frag_threshold)

# # Additional setup commands
# additional_lines = [
#     'module purge\n',
#     'module load shared\n',
#     'module load 2024\n',
#     'module load slurm\n',
#     'module load intel/mkl/64/2020/4.304\n',
#     'module load Python/3.11.5-GCCcore-13.2.0\n', 
#     'export SCM_TMPDIR="$TMPDIR"\n',
#     'export AMSHOME="/scistor/tc/afr790/amshome"\n',
#     'export AMSBIN="/scistor/tc/afr790/amshome/bin"\n',
#     'export AMSRESOURCES="$AMSHOME/atomicdata"\n',
#     'export SCMLICENSE="/scistor/tc/huw587/ADF_trunk/license.txt"\n',
#     f'mkdir -pv {outdir}\n'
# ]

# # Create directories
# current_dir = os.getcwd()
# job_dir = os.path.join(current_dir, 'jobs')
# out_dir = os.path.join(current_dir, 'out')
# err_dir = os.path.join(current_dir, 'err')

# os.makedirs(job_dir, exist_ok=True)
# os.makedirs(out_dir, exist_ok=True)
# os.makedirs(err_dir, exist_ok=True)

# # Generate SLURM job script
# job_file = os.path.join(job_dir, f"{output_name}.job")

# with open(job_file, 'wt') as fh:
#     fh.writelines(
#         [
#             "#!/bin/bash\n",
#             f"#SBATCH --job-name=multi_calc\n",
#             f"#SBATCH --output={os.path.join(out_dir, 'multi_calc_%A_%a.out')}\n",
#             f"#SBATCH --error={os.path.join(err_dir, 'multi_calc_%A_%a.err')}\n",
#             f"#SBATCH --time={time}\n",
#            # f"#SBATCH --mem={mem}\n",
#             f"#SBATCH --partition={partition}\n",
#             f"#SBATCH --ntasks-per-node=32\n",
#             f"#SBATCH --ntasks-per-core=1\n",
#             f"#SBATCH -N 1\n",
#             f"#SBATCH --mail-type=NONE\n",
#             f"#SBATCH --array=0-{num_jobs - 1}%{njobs}\n",
#         ] + additional_lines + [
#             "echo $TMPDIR\n",
#             "unset SLURM_JOB_ID\n",
#             "cd $TMPDIR\n",
#             f"cp -rfv {outdir}/run_cluster_luna/{script} {outdir}/geometry $TMPDIR\n",
#             "\n"
#             "# Compute molecule, method, and separation index\n",
#             f"molecule_index=$((SLURM_ARRAY_TASK_ID / ({len(add_virt_threshold)} * {len(frag_threshold)})))\n",
#             f"av_index=$(( (SLURM_ARRAY_TASK_ID / {len(frag_threshold)}) % {len(add_virt_threshold)} ))\n",
#             f"ft_index=$((SLURM_ARRAY_TASK_ID % {len(frag_threshold)}))\n",
#             f"molecules=({' '.join(molecules)})\n",
#             f"av=({' '.join(float(add_virt_threshold))})\n",
#             f"ft=({' '.join(map(str, float(frag_threshold)))})\n",
#             "molecule=${molecules[$molecule_index]}\n",
#             "add_virt_thr=${av[$av_index]}\n",
#             "frag_thr=${ft[$ft_index]}\n",
#             "\n",
#             f"$AMSBIN/amspython -u {script} $molecule $add_virt_thr $frag_thr\n",
#             f"mv plams_workdir plams_workdir_${{molecule}}_${{add_virt_thr}}_${{frag_thr}}\n",
#             f"mkdir -pv {outdir1}/$molecule/$add_virt_thr\n",
#             f"mkdir -pv {outdir_coupling}/$molecule/$add_virt_thr\n",
#             f"cp -rfv plams_workdir_${{molecule}}_${{add_virt_thr}}_${{frag_thr}} {outdir1}/$molecule/$add_virt_thr\n",
#             f"cp -rfv output* {outdir_coupling}/$molecule/$add_virt_thr\n"
#         ]
#     )

# if RUN:
#     os.system("sbatch " + job_file)


