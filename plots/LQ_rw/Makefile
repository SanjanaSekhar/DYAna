CXXFLAGS=$(shell root-config --cflags)
GLIBS=$(shell root-config --glibs)


all: LHA_test add_pdf_weight

#LQ_rw: LQ_rw_denom.C
	#g++ $(CXXFLAGS) $^ -o $@ -I$(HOME)/local/include/ $(GLIBS) `lhapdf-config --cflags --ldflags`

add_pdf_weight: add_pdf_weight.C
	g++ $(CXXFLAGS) $^ -o $@ -O 0 -g -I$(HOME)/local/include/ $(GLIBS) `lhapdf-config --cflags --ldflags`

LHA_test: LHA_test.C
	g++ $(CXXFLAGS) $^ -o $@ -I$(HOME)/local/include/ $(GLIBS) `lhapdf-config --cflags --ldflags`
