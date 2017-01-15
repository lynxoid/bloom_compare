BIN=bin
LIB=-L /usr/local/lib/ # -L $(HOME)/tbb2017_20160916oss/lib/
FLAGS=-O3 -std=c++11

all: create partow vallentine

create:
	mkdir -p $(BIN)
# compile all 3 variants: Vallentine, Partow, Marcais
partow:
	clang++ $(FLAGS) -o $(BIN)/$@ -I include/partow/ src/bloom_filter_example01.cpp
vallentine:
	clang++ $(FLAGS) -o $(BIN)/$@ -I include/vallentine/ src/vallentine.cpp -lbf $(LIB)
marcais:
	clang++ $(FLAGS) -o $(BIN)/$@ -I include/marcais/ src/marcais.cpp
