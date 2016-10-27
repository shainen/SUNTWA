from string import Template
from os import getcwd

project=getcwd().split('/')[-1]

runname=getcwd().split('/')[-2]

topic=getcwd().split('/')[-3]

qsubfile = Template("""
#!/bin/sh
#PBS -j y
#PBS -l mem=400mb
#PBS -V
#PBS -m ea
#PBS -M shainen@gmail.com

RUN_NAME=${rname}
SCRATCH_DIR=/projectnb/twambl/$$RUN_NAME/r$$SGE_TASK_ID
LOCAL_DIR=/project/twambl/${tp}

mkdir -p $$SCRATCH_DIR

# Copy them from this directory to data director
cd $$LOCAL_DIR/$$RUN_NAME
cp -r ${prj} $$SCRATCH_DIR/

# Move to the data directory
cd $$SCRATCH_DIR/

# Run the script
cd ${prj}
echo "SeedRandom[$$SGE_TASK_ID]" > randomSeed.wl
cd ..
time math -script ${prj}/runTWA.wl

# Remove the now-useless files
rm -r ${prj} 
""")

with open("../job"+runname+".qsub", "w") as f:
#f=open("../job"+runname+".qsub", "w")
    f.write(qsubfile.substitute(rname=runname,prj=project,tp=topic))

