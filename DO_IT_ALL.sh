set cmonkey_jobid = 322191 ## job id of running cmonkey2 job array

python egrin2-tools/postproc/fimo.py --user dreiss --csh --organism_name mtu \
    --seqs_file cache/Mycobacterium_tuberculosis_H37Rv_NC_000962.2

python egrin2-tools/postproc/coding_fracs.py --user dreiss --csh --organism_name mtu \
    --num_runs 300 --features_file cache/Mycobacterium_tuberculosis_H37Rv_features

python egrin2-tools/postproc/tomtom.py --user dreiss --csh --prefix mtu-out- \
    --dir . --targetdir tomtom_out

rpl 'qsub ' "qsub -q baliga -hold_jid $cf_jobid " qsub_tomtom.sh

set fimo_jobid = `qsub -q baliga -hold_jid $cmonkey_jobid qsub_fimo.sh | awk '{print $3}'`

set cf_jobid = `qsub -q baliga -hold_jid $fimo_jobid qsub_coding_fracs.sh | awk '{print $3}'`

source qsub_tomtom.sh

set tt_jobid = `qstat -f -u dreiss -j '*tomtom*' | grep 'job_number' | sort -n -k2 | tail -1 | awk '{print $2}'`

## need to ping the queue and wait until the last tomtom job is done...
set ntomtomjobs = `qstat -f -u dreiss -j '*tomtom*' | grep job_number | wc -l`
while ( $ntomtomjobs != "0" )
  echo "waiting... $ntomtomjobs ..."
  sleep 10
end

python egrin2-tools/postproc/motif_clustering.py --input_dir tomtom_out
