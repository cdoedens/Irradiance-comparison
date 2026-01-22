import subprocess
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

script_to_run = '041_weather_features'

dataset = 'barra-r2'
wf = 'cycmask'


# Generate a unique file name based on iteration
joboutdir = '/home/548/cd3022/repos/Irradiance-comparisons/scripts/qsub_jobs/'

for year in range(2016, 2023):
    job_script_filename = joboutdir + f'{script_to_run}___{dataset}__{wf}__{year}.qsub'

    # Open the file for writing
    with open(job_script_filename, "w") as f3:
        f3.write('#!/bin/bash \n')
        f3.write('#PBS -l walltime=1:00:00 \n')
        f3.write('#PBS -l mem=96GB \n')
        f3.write('#PBS -l ncpus=48 \n')
        f3.write('#PBS -l jobfs=10GB \n')
        f3.write('#PBS -l storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74+gdata/su28+scratch/er8 \n')
        f3.write('#PBS -l other=hyperthread \n')
        f3.write('#PBS -q normal \n')
        f3.write('#PBS -P er8 \n')
        f3.write(f'#PBS -o /home/548/cd3022/repos/Irradiance-comparisons/logs/{script_to_run}_{dataset}_{wf}_{year}.oe \n')
        f3.write('#PBS -j oe \n')
        f3.write('cd /home/548/cd3022/repos/Irradiance-comparisons \n')
        f3.write('source ./env.sh \n')
        f3.write(f'python3 scripts/{script_to_run}.py {dataset} {wf} {year}\n')


    # Submit the generated script to the job scheduler (PBS) using qsub
    try:
        # Run the qsub command and submit the script
        subprocess.run(['qsub', job_script_filename], check=True)
        print(f"Job script {job_script_filename} submitted successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job script {job_script_filename}: {e}")