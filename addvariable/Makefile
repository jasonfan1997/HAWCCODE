CPP=g++
#CFLAGS=-O2 -Wall -g -std=c++0x
CFLAGS=-O2 -std=c++0x
targets: add

add: add.cpp xgboost_classifier.o xgboost_classifier.h
	$(CPP) add.cpp xgboost_classifier.cpp -I${XCDFROOT}/include -lxcdf -L${XCDFROOT}/lib -I${HAWC_INSTALL}/include/ -lhawcnest -L${HAWC_INSTALL}/lib -I${BOOSTROOT}/include/ -L${BOOSTROOT}/lib -o addproba.exe $(CFLAGS)
	
xgboost_classifier.o: xgboost_classifier.h xgboost_classifier.cpp
	$(CPP) $(CFLAGS) -c xgboost_classifier.cpp
	
	
.PHONY: clean

clean:
	rm -vf addproba.exe