#!/bin/zsh
echo "0.03 0.04 0.01\n0.05 0.02 0.01" > sample_pts
~/local/bin/interactive_emulator interactive_mode multi_snapshot_file < sample_pts 
## now repeat with cov-mat-estimation
#~/local/bin/interactive_emulator interactive_mode multi_snapshot_file --covmat < sample_pts
