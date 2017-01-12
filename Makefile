# compile all 3 variants: Vallentine, Partow, Marcais
clang++ -O3 -o bin/partow -I include/partow/ src/bloom_filter_example01.cpp

clang++ -O3 -o bin/vallentine src/vallentine.cpp -l libbf

clang++ -O3 -o bin/marcais -I include/marcais/ src/marcais.cpp

# repeat every run 10 times
for i in {1..10}
do
	# clean the cache every time
	bin/vallentine $READS
done
