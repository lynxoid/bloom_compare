BIN=bin
LIB=-L /usr/local/lib/ # -L $(HOME)/tbb2017_20160916oss/lib/

all: create partow vallentine marcais

create:
	mkdir -p $(BIN)
# compile all 3 variants: Vallentine, Partow, Marcais
partow:
	clang++ -O3 -o $(BIN)/$@ -I include/partow/ src/bloom_filter_example01.cpp
vallentine:
	clang++ -O3 -o $(BIN)/$@ -I include/vallentine/ src/vallentine.cpp -lbf $(LIB)
marcais:
	clang++ -O3 -o $(BIN)/$@ -I include/marcais/ src/marcais.cpp
run:
	READS=xyz.fasta
	# repeat every run 10 times
	for method in vallentine partow marcais
	do
		echo $method
		for i in {1..10}
		do
			# clean the cache every time
			$(BIN)/$(METHOD) $READS
		done
	done
