if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters, exiting."
    exit 1
fi

jid1=$(sbatch separate_channel_id.sh $1)
jid1=${jid1##* }
echo $jid1
jid2=$(sbatch --dependency=afterok:$jid1 detrend_id.sh $1)
jid2=${jid2##* }
echo $jid2
jid3=$(sbatch --dependency=afterok:$jid2 hemo_correct_id.sh $1)
jid3=${jid3##* }
echo $jid3
jid4=$(sbatch --dependency=afterok:$jid3 post_hemo_id.sh $1)
jid4=${jid4##* }
echo $jid4
