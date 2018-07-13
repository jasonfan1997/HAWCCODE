CPP=g++
CFLAGS= -O2 -Wall -g -std=c++0x

targets: add

add: add.cpp xgboost_classifier.o xgboost_classifier.h
	$(CPP) add.cpp xgboost_classifier.cpp -I${XCDFROOT}/include -lxcdf -L${XCDFROOT}/lib -o addproba.exe $(CFLAGS)
	
xgboost_classifier.o: xgboost_classifier.h xgboost_classifier.cpp
	$(CPP) $(CFLAGS) -c xgboost_classifier.cpp
	
	
.PHONY: clean

clean:
	rm -vf addproba.exe