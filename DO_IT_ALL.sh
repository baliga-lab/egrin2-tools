set cmonkey_jobid = 322191 ## job id of running cmonkey2 job array

set ORGANISM_LONG = "Mycobacterium_tuberculosis_H37Rv"
set ORGANISM_SHORT = "mtu"
set ORGANISM_GENOME = "${ORGANISM_LONG}_NC_000962.2"
set COMMONARGS = "--user dreiss --csh"

python egrin2-tools/postproc/fimo.py ${COMMONARGS} --organism "${ORGANISM_SHORT}" \
    --genome cache/${ORGANISM_GENOME}

set fimo_jobid = `qsub -q baliga -hold_jid $cmonkey_jobid qsub_fimo.sh | awk '{print $3}' | cut -d'.' -f1`
echo $fimo_jobid

python egrin2-tools/postproc/coding_fracs.py ${COMMONARGS} --organism "${ORGANISM_SHORT}" \
    --features cache/${ORGANISM_LONG}_features

set cf_jobid = `qsub -q baliga -hold_jid $fimo_jobid qsub_coding_fracs.sh | awk '{print $3}' | cut -d'.' -f1`
echo $cf_jobid

python egrin2-tools/postproc/tomtom.py ${COMMONARGS} --prefix "${ORGANISM_SHORT}-out-"

rpl 'qsub ' "qsub -q baliga -hold_jid $cf_jobid " qsub_tomtom.sh

source qsub_tomtom.sh

set tt_jobid = `qstat -f -u dreiss -j '*tomtom*' | grep 'job_number' | sort -n -k2 | tail -1 | awk '{print $2}'`
echo $tt_jobid

## need to ping the queue and wait until the last tomtom job is done...
set ntomtomjobs = `qstat -f -u dreiss -j '*tomtom*' | grep job_number | wc -l`
while ( $ntomtomjobs != "0" )
  echo "waiting... $ntomtomjobs ..."
  sleep 10
end

## This is NOT run on the SGE cluster -- it is run locally so run it on osiris, not aegir!
python egrin2-tools/postproc/motif_clustering.py --input_dir tomtom_out
