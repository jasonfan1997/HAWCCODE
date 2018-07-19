#include <iostream>

//#include "TMath.h"

#include <fstream>

#include <string>
#include <vector>
#include <xcdf/XCDF.h>
#include <xcdf/XCDFField.h>
#include <xcdf/XCDFFile.h>
#include <xcdf/utility/XCDFUtility.h>
#include <hawcnest/CommandLineConfigurator.h>
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
	  CommandLineConfigurator cl(
					"Append probability to XCDF "
					"files."
					);
	cl.AddPositionalOption<std::vector<std::string> >(
					"i", "Input XCDF files.  Default is stdin."
					);
	cl.AddOption<std::string>(
					"out-file,o", "/dev/stdout", "Output XCDF file."
					);
	if (!cl.ParseCommandLine(argc, argv)) return EXIT_FAILURE;

	std::vector<std::string> inFileNames;
	if (cl.HasFlag("i"))
			inFileNames = cl.GetArgument<std::vector<std::string> >("i");
	else inFileNames.push_back("/dev/stdin");

	XCDFFile recFile(cl.GetArgument<std::string>("out-file").c_str(), "w");
	FieldCopyBuffer buffer(recFile);
	XCDFFloatingPointField newfield = recFile.AllocateFloatingPointField("rec.proba", 0.000000001);
	  double nnHit=0;
	  double CCxPE40=0;
	  double PPINC=0;
	  double llogNNEnergyV2=0;
	  double ddisMax=0;
	  double LLDFAmp=0;
	  double LLDFChi2=0;
	  double nnChAvail=0;
	  double nnHitSP20=0;
	for(int i=0,n=inFileNames.size();i<n;i++)
	{
		
		  XCDFFile rec;
		  rec.Open(inFileNames[i], "r");
		  
		  XCDFUnsignedIntegerField  nHit = rec.GetUnsignedIntegerField("rec.nHit");
		  XCDFFloatingPointField CxPE40 = rec.GetFloatingPointField("rec.CxPE40");
		  XCDFFloatingPointField PINC = rec.GetFloatingPointField("rec.PINC");
		  XCDFFloatingPointField logNNEnergyV2 = rec.GetFloatingPointField("rec.logNNEnergyV2");
		  XCDFFloatingPointField disMax = rec.GetFloatingPointField("rec.disMax");
		  XCDFFloatingPointField LDFAmp = rec.GetFloatingPointField("rec.LDFAmp");
		  XCDFFloatingPointField LDFChi2 = rec.GetFloatingPointField("rec.LDFChi2");
		  XCDFUnsignedIntegerField  nChAvail = rec.GetUnsignedIntegerField("rec.nChAvail");
		  XCDFUnsignedIntegerField nHitSP20 = rec.GetUnsignedIntegerField("rec.nHitSP20");

		  //FieldCopyBuffer buffer(recFile);
		  
		  std::set<std::string> fields;

		  GetFieldNamesVisitor getFieldNamesVisitor(fields);
		  rec.ApplyFieldVisitor(getFieldNamesVisitor);
		  fields.erase("rec.proba");
		  SelectFieldVisitor selectFieldVisitor(rec, fields, buffer);
		  rec.ApplyFieldVisitor(selectFieldVisitor);

		  //ifstream file ( argv[3]);
		  string value;
		   while ( rec.Read() )
		  {
			std::vector<float> sample (7, 0);
			nnHit=*nHit;
			CCxPE40=*CxPE40;
			PPINC=*PINC;
			llogNNEnergyV2=*logNNEnergyV2;
			ddisMax=*disMax;
			LLDFAmp=*LDFAmp;
			LLDFChi2=*LDFChi2;
			nnChAvail=*nChAvail;
			nnHitSP20=*nHitSP20;
			if(nnHit==0||CCxPE40==0){
				sample[0]=0;
			}
			else{
			sample[0]=log10(CCxPE40/nnHit);
			}
			sample[1]=PPINC;
			sample[2]=llogNNEnergyV2;
			sample[3]=ddisMax;
			sample[4]=LLDFAmp;
			sample[5]=LLDFChi2;
			if(nnChAvail==0||nnHitSP20==0){
				sample[6]=0;
			}
			else{
			sample[6]=nnHitSP20/nnChAvail;
			}			
					
			std::vector<float> res = xgb_classify(sample);
			double probability= proba(res[0]);
			newfield << probability;
			
			buffer.CopyData();

			 recFile.Write();
		   }
		   rec.Close();
	}
    recFile.Close();

  return 0;
}
