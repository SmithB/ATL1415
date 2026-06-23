# edit the run_arctic* to run only SV for a test

# test:

run_arctic_prelim.sh default_args/rel_005_0329.txt default_args/discover.txt  default_args/quarterly.txt IS
run_arctic_matched.sh default_args/rel_005_0329.txt default_args/discover.txt  default_args/quarterly.txt IS
run_arctic_mosaic.sh default_args/rel_005_0329.txt default_args/discover.txt  default_args/quarterly.txt IS


# For real:

# prelim
bash scripts/run_arctic_prelim.sh default_args/rel_005_0329.txt default_args/discover.txt  default_args/quarterly.txt

# check for errors or timeouts
for j in RA IS CN CS SV; do echo $j; pushd $j"_prelim"; grep -i canc *; grep -i kill *; log_parser.py .; popd; done

#matched:
bash scripts/run_arctic_matched.sh default_args/rel_005_0329.txt default_args/discover.txt  default_args/quarterly.txt

for j in RA IS CN CS SV; do echo $j; pushd $j"_matched"; grep -i canc *; grep -i kill *; log_parser.py .; popd; done


bash scripts/run_arctic_mosaic.sh default_args/rel_005_0329.txt default_args/discover.txt  default_args/quarterly.txt

bash scripts/run_arctic_to_nc.sh default_args/rel_005_0329.txt default_args/discover.txt  default_args/quarterly.txt
