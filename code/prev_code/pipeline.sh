pipeline.sh

# 2023 data
# wget -r -np -R "index.html*" --http-user=user --http-password=nkAHPsUFgX https://cobb.sr.unh.edu/managed/231117_A01346_0123_BHH2GFDRX3_16Mer111723-BB-Tt_MassDiatoms_2023/reads
## re-run with qc samples
code/qiime2_denoise.sh  \
    BB-Tt_MassDiatoms_2023 \
    data/algae/runs \
    raw-data/cobb.sr.unh.edu/managed/231117_A01346_0123_BHH2GFDRX3_16Mer111723-BB-Tt_MassDiatoms_2023/reads \
    4 \
    rbcl \
    paired &> data/runlogs/runlog.BB-Tt_MassDiatoms_2023

code/qiime2_hybrid-learn.sh \
    BB-Tt_MassDiatoms_2023 \
    data/algae/runs \
    raw-data/cobb.sr.unh.edu/managed/231117_A01346_0123_BHH2GFDRX3_16Mer111723-BB-Tt_MassDiatoms_2023/reads \
    12 \
    rbcl &> data/runlogs/runlog.BB-Tt_MassDiatoms_2023.out

# wget -r -np -R "index.html*" --http-user=user --http-password=CmUdvuRPWN https://cobb-data.sr.unh.edu/projects/251028_A01346_0193_AHCMCLDRX7_AMP-102825-TetraTech-MassDiatoms-2025/reads
## re-run with qc samples
code/qiime2_denoise.sh  \
    MassDiatoms-2025 \
    algae/runs \
    cobb-data.sr.unh.edu/projects/251028_A01346_0193_AHCMCLDRX7_AMP-102825-TetraTech-MassDiatoms-2025/reads \
    4 \
    rbcl \
    paired &> algae/runlogs/runlog.MassDiatoms-2025

# copy old data to new location
# 2023
cp -r BB-Tt_MassDiatoms_2023 /home/users/jtm1171/algae/runs/
# 2024
cp -r AlgaeME-rbcLNX031523 /home/users/jtm1171/algae/runs/










