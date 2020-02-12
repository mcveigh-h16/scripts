cat $1 | awk '{printf "select NX_id,locus from NX_hist where NX_id = \"%s\"\ngo -m bcp\n",$1;}' | sqsh -S REFTRACK -D RefTrack -U anyone -P 'allowed'
