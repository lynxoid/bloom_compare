# run experiments
BIN=bin

make all

READS=/Users/lynxoid/dev/seq_data/hs_alt_CHM1_1.1_chr20_100K.fa
# repeat every run 10 times
for method in vallentine partow
do
	echo $method
	for i in {1..2}
	do
		# clean the cache every time
		time $BIN/$method $READS 20
	done
done
