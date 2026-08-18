# test:
# N.B.  reinstall ATL1415 to make sure run_arctic* is up to date
run_arctic_prelim.sh default_args/latest_release.txt default_args/discover.txt  default_args/quarterly.txt IS
run_arctic_matched.sh default_args/latest_release.txt default_args/discover.txt  default_args/quarterly.txt IS
run_arctic_mosaic.sh default_args/latest_release.txt default_args/discover.txt  default_args/quarterly.txt IS
run_arctic_to_nc.sh default_args/latest_release.txt default_args/discover.txt  default_args/quarterly.txt IS




# For real:

# prelim
bash scripts/run_arctic_prelim.sh default_args/latest_release.txt default_args/discover.txt  default_args/quarterly.txt

# check for errors or timeouts
for j in RA IS CN CS SV; do echo $j; slurm_run_status.py $j"_prelim"; done

#matched:
bash scripts/run_arctic_matched.sh default_args/latest_release.txt default_args/discover.txt  default_args/quarterly.txt

for j in RA IS CN CS SV; do echo $j; slurm_run_status.py $j"_matched"; done


bash scripts/run_arctic_mosaic.sh default_args/latest_release.txt default_args/discover.txt  default_args/quarterly.txt
bash scripts/run_arctic_to_nc.sh default_args/latest_release.txt default_args/discover.txt  default_args/quarterly.txt


# monthly: 
run_arctic_prelim.sh default_args/latest_release.txt default_args/discover.txt  default_args/monthly.txt
run_arctic_matched.sh default_args/latest_release.txt default_args/discover.txt  default_args/monthly.txt 
run_arctic_mosaic.sh default_args/latest_release.txt default_args/discover.txt  default_args/monthly.txt
run_arctic_to_nc.sh default_args/latest_release.txt default_args/discover.txt  default_args/monthly.txt
