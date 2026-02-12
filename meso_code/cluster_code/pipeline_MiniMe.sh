if [ "$#" -ne 2 ]; then

    echo "Illegal number of parameters, exiting."

    exit 2

fi

jid1=$(sbatch separate_channel.sh $1 $2)
jid1=${jid1##* }
echo $jid1
jid2=$(sbatch --dependency=afterok:$jid1 hemo_correct.sh $1 $2)
jid2=${jid2##* }
echo $jid2


