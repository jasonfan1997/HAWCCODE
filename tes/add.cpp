#include <iostream>

//#include "TMath.h"

#include <fstream>

#include <string>
#include <vector>
#include <xcdf/XCDF.h>
#include <xcdf/XCDFField.h>
#include <xcdf/XCDFFile.h>
#include <xcdf/utility/XCDFUtility.h>
#include "xgboost_classifier.h"
#include <chrono>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <iterator>
using namespace std;

double proba(float logproba)
{
	return exp(logproba)/(1+exp(logproba));
	
}
int main(int argc, char** argv) {
  //char* name=getFileName(argv[0]);
  XCDFFile rec;
  rec.Open(argv[1], "r");
  XCDFFile recFile(argv[2],"w");
  XCDFUnsignedIntegerField  a = rec.GetUnsignedIntegerField("rec.nHit");
  XCDFFloatingPointField b = rec.GetFloatingPointField("rec.CxPE40");
  XCDFFloatingPointField c = rec.GetFloatingPointField("rec.PINC");
  XCDFFloatingPointField d = rec.GetFloatingPointField("rec.logNNEnergyV2");
  XCDFFloatingPointField e = rec.GetFloatingPointField("rec.disMax");
  XCDFFloatingPointField f = rec.GetFloatingPointField("rec.LDFAmp");
  XCDFFloatingPointField g = rec.GetFloatingPointField("rec.LDFChi2");
  XCDFUnsignedIntegerField  h = rec.GetUnsignedIntegerField("rec.nChAvail");
  XCDFUnsignedIntegerField i = rec.GetUnsignedIntegerField("rec.nHitSP20");
  double aa=0;
  double bb=0;
  double cc=0;
  double dd=0;
  double ee=0;
  double ff=0;
  double gg=0;
  double hh=0;
  double ii=0;
  std::vector<std::vector<float>> samples;


  FieldCopyBuffer buffer(recFile);
  XCDFFloatingPointField newfield = recFile.AllocateFloatingPointField("rec.proba", 0.000000001);
  std::set<std::string> fields;

  GetFieldNamesVisitor getFieldNamesVisitor(fields);
  rec.ApplyFieldVisitor(getFieldNamesVisitor);

  SelectFieldVisitor selectFieldVisitor(rec, fields, buffer);
  rec.ApplyFieldVisitor(selectFieldVisitor);

  //ifstream file ( argv[3]);
  string value;
   while ( rec.Read() )
  {
	std::vector<float> sample (7, 0);
	aa=*a;
	bb=*b;
	cc=*c;
	dd=*d;
	ee=*e;
	ff=*f;
	gg=*g;
	hh=*h;
	ii=*i;
	if(aa==0||bb==0){
		sample[0]=0;
	}
	else{
	sample[0]=log10(bb/aa);
	}
	sample[1]=cc;
	sample[2]=dd;
	sample[3]=ee;
	sample[4]=ff;
	sample[5]=gg;
	if(aa==0||ii==0){
		sample[6]=0;
	}
	else{
	sample[6]=ii/hh;
	}
	
	
	std::vector<float> res = xgb_classify(sample);
	double probability= proba(res[0]);
	newfield << probability;
	
    buffer.CopyData();

     recFile.Write();
   }
   rec.Close();
    recFile.Close();

  return 0;
}
