CXXFLAGS=$(shell root-config --cflags)
GLIBS=$(shell root-config --glibs)


all: LHA_test add_pdf_weight

add_pdf_weight: add_pdf_weight.C
	g++ $(CXXFLAGS) $^ -o $@ -I$(HOME)/local/include/ $(GLIBS) `lhapdf-config --cflags --ldflags`

LHA_test: LHA_test.C
	g++ $(CXXFLAGS) $^ -o $@ -I$(HOME)/local/include/ $(GLIBS) `lhapdf-config --cflags --ldflags`
