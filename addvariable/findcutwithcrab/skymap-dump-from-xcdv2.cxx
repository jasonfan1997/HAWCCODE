#include <xcdf/utility/XCDFUtility.h>
#include <xcdf/utility/Expression.h>

#include <data-structures/time/TimeStamp.h>
#include <hawcnest/HAWCNest.h>
#include <track-fitter/ZenithAlignment.h>
#include <astro-service/StdAstroService.h>
#include <detector-service/ConfigDirDetectorService.h>
#include <hawcnest/CommandLineConfigurator.h>

#include <TH2F.h>
#include <TFile.h>

#include <fstream>
#include <set>
#include <iostream>

#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <fitshandle.h>
#include <healpix-cxx/cxxsupport/vec3.h>

typedef Healpix_Map<double> HMap;

using namespace std;
bool fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}
double Delta(double ra1, double dec1, double ra2, double dec2) {
  double x1 = cos(ra1)*sin(3.1415926536/2.-dec1);
  double y1 = sin(ra1)*sin(3.1415926536/2.-dec1);
  double z1 = cos(3.1415926536/2.-dec1);

  double x2 = cos(ra2)*sin(3.1415926536/2.-dec2);
  double y2 = sin(ra2)*sin(3.1415926536/2.-dec2);
  double z2 = cos(3.1415926536/2.-dec2);
  
  double a = acos(x1*x2+y1*y2+z1*z2);
  return a;
}

double lima(double r_on, double r_off, double alpha)
{
  double signi,diff,term1,term2,term3;

  signi = 0.;
  term1 = ((1.+alpha)/alpha)*(r_on /(r_on + r_off));
  term2 =  (1.+alpha)       *(r_off/(r_on + r_off));

  if (term1 > 0. && term2 > 0.) {
    term3 = r_on  * log(term1) + r_off * log(term2);
    if (term3 < 0.) term3 = 0.;
    signi = sqrt(2.) * sqrt(term3);
    diff = r_on - alpha * r_off;
    if (diff < 0.)  signi=(-1)*signi;
  }

  return signi;
}


int main(int argc, char** argv){
  CommandLineConfigurator cl;

  cl.AddOption<string>("input", "HAWC rec XCDF file");

  cl.AddOption<string>("zenith-alignment-file,z", "", "zenith aligment XML");
  cl.AddOption<string>("output,o", "", "output ROOT file");
  cl.AddOption<string>("report,r", "", "report file");

  cl.AddFlag("healpix","healpix all-sky map instead of ROOT output");

  cl.AddOption<double>("azimuth-correction-value,a",0., "azimuth correction value in degree");

  cl.AddOption<double>("pinc", "pinc cut");
  cl.AddFlag("dopinc", "do pinc cut");

  cl.AddOption<double>("dash", "dash cut");
  cl.AddFlag("dodash", "do dash cut");

  cl.AddOption<double>("age", "age cut");
  cl.AddFlag("doage", "do age cut");

  cl.AddOption<double>("licompact", "licompact cut");
  cl.AddFlag("dolicompact", "do licompact cut");
  
  cl.AddPositionalOption<vector<double> >("proba", "proba cut");
  cl.AddFlag("doproba", "do proba cut");
  
  cl.AddOption<double>("ldfchi2", "LDF Chi2 cut");
  cl.AddFlag("doldfchi2", "do LDF Chi2 cut");
  
  cl.AddOption<string>("customname", "custom name");
  cl.AddOption<double>("customcut", "custom cut value");
  cl.AddFlag("docustom", "do custom cut");
  
  cl.AddOption<double>("fidu", "Fiducial Area cut");
  cl.AddFlag("dofidu", "do Fiducial Area cut");

  cl.AddFlag("select", "Save GPS times to a file called selected.times for use with aerie-apps-xcdf-time-selection");

  cl.AddOption<double>("costhmin", "Minimum 1./cos(th) Cut");
  cl.AddOption<double>("costhmax", "Maximum 1./cos(th) Cut");
  cl.AddFlag("docosth", "do 1./cos(th) cut");

  cl.AddFlag("poisson", "Possion statistics when computing significance (default is Gaussian)");

  if (!cl.ParseCommandLine(argc, argv)) {
    return 1;
  }

  string recfile = cl.GetArgument<string>("input");
  string report = "report";
  string output = "output";
  if (cl.HasFlag("report")) {
    report = cl.GetArgument<string>("report");
    printf("report file set to: %s\n",report.c_str());
  }
  if (cl.HasFlag("output")) {
    output = cl.GetArgument<string>("output");
    printf("output file set to: %s\n",output.c_str());
  }
  bool docustom =  false;
  string customcutname;
  double customcut=0;
  if (cl.HasFlag("docustom")) {
	 docustom=true;
	customcutname = cl.GetArgument<string>("customname");
	customcut=cl.GetArgument<double>("customcut");
	printf("customcutname: %s\n",output.c_str());
  }
  
  
  bool doproba =  false;
  vector<double> proba_cut_vec;
  double proba_cut_value=0;
  if(cl.HasFlag("doproba")){
	doproba = true;
	proba_cut_vec = cl.GetArgument<vector<double> >("proba");
  }
  else
	  proba_cut_vec.push_back(0.0);
	for(int i=0,n=proba_cut_vec.size();i<n;i++)
	{
	  proba_cut_value=proba_cut_vec[i];
  // go get stuff out of XCDF file
	  XCDFFile f(recfile.c_str(),"r");

	  XCDFFloatingPointField ra_x = f.GetFloatingPointField("rec.ra");
	  
	  XCDFFloatingPointField dec_x = f.GetFloatingPointField("rec.dec");
	  XCDFFloatingPointField theta_x = f.GetFloatingPointField("rec.zenithAngle");
	  
	  XCDFFloatingPointField phi_x = f.GetFloatingPointField("rec.azimuthAngle");

	  XCDFFloatingPointField pinc_x = f.GetFloatingPointField("rec.PINC");
	  XCDFFloatingPointField proba_x = f.GetFloatingPointField("rec.proba");
	  XCDFFloatingPointField ldfchi2_x = f.GetFloatingPointField("rec.LDFChi2");
	  //XCDFFloatingPointField dash_x = f.GetFloatingPointField("rec.DASH");
	  XCDFFloatingPointField age_x = f.GetFloatingPointField("rec.LDFAge");
	  XCDFUnsignedIntegerField nhitsp20_x = f.GetUnsignedIntegerField("rec.nHitSP20");
	  XCDFFloatingPointField cxpe40_x = f.GetFloatingPointField("rec.CxPE40");

	  XCDFUnsignedIntegerField gpsSec_x = f.GetUnsignedIntegerField("rec.gpsSec");
	  XCDFUnsignedIntegerField gpsNanosec_x = f.GetUnsignedIntegerField("rec.gpsNanosec");
	  XCDFUnsignedIntegerField runid_x = f.GetUnsignedIntegerField("rec.runID");

	  XCDFUnsignedIntegerField coreFiduScale_x = f.GetUnsignedIntegerField("rec.coreFiduScale");
	  XCDFUnsignedIntegerField angleFitStatus_x = f.GetUnsignedIntegerField("rec.angleFitStatus");
	  XCDFUnsignedIntegerField nChAvail_x = f.GetUnsignedIntegerField("rec.nChAvail");
	  XCDFUnsignedIntegerField nChTot_x = f.GetUnsignedIntegerField("rec.nChTot");
	  XCDFFloatingPointField custom_x = f.GetFloatingPointField(customcutname.c_str());

	  // make a Nest instance to get access to HAWC services
	  HAWCNest nest;
	  nest.Service<ConfigDirDetectorService>("det");
	  nest.Service<ZenithAlignment>("zenithAlignment")
		("dataLocation", "zenith-SPhigh-on-liff-geom-single.xml");
	  nest.Service<StdAstroService>("astroService");

	  nest.Configure();
	  ZenithAlignment& alignment = 
		GetService<ZenithAlignment>("zenithAlignment");
	  AstroService& astro = 
		GetService<AstroService>("astroService");
	  DetectorService& detSv    = GetService<DetectorService>("det");

	  //  TH2F* sky = new TH2F("skymap","skymap",100,80,87,100,19,25);
	  TH2D* sky = new TH2D("skymap","sky",400,79.63,87.63,400,18.,26.);
	  TH2D* allsky = new TH2D("allskymap","sky",1200,0.,360.,300,-25.,65.);
	  TH2D* allsky2 = new TH2D("allskymap2","sky",360,0.,360.,90,-25.,65.);
	  TH1D* evexcess = new TH1D("evexcess","Excess by radial bin",200,0.,5.);
	  TH1D* evbins   = new TH1D("evbins","Number of map bins in radial bin",200,0.,5.);
	  TH1D* evdensity= new TH1D("evdensity","Event Density",200,0.,5.);
	  TH1D* gadensity= new TH1D("gadensity","Gamma-Ray Density (Background Subtracted)",200,0.,5.);
	  TH1D* evint    = new TH1D("evint","Integral Event Count",200,0.,5.);
	  TH1D* bckint   = new TH1D("bckint","Integral Background Count",200,0.,5.);
	  TH1D* exint    = new TH1D("xvint","Integral Event Excess",200,0.,5.);
	  TH1D* evsig    = new TH1D("evsig","Significance of Excess",200,0.,5.);
	  TH1D* evrsq    = new TH1D("evrsq","Rsq",200,0.,1.);

	  HMap skymap(1024,RING,SET_NSIDE);

	  bool dopinc =  false;
	  double pinc_cut_value=0;
	  if(cl.HasFlag("dopinc")){
		dopinc = true;
		pinc_cut_value = cl.GetArgument<double>("pinc");
	  }


	  bool doldfchi2 =  false;
	  double ldfchi2_cut_value=0;
	  if(cl.HasFlag("doldfchi2")){
		doldfchi2 = true;
		ldfchi2_cut_value = cl.GetArgument<double>("ldfchi2");
	  }

	  bool dofidu = false;
	  unsigned int  fidu_cut_value=80;
	  if(cl.HasFlag("dofidu")){
		dofidu = true;
		fidu_cut_value = (unsigned int) abs(cl.GetArgument<double>("fidu"));
	  }

	  bool doselect = false;
	  FILE *fpselect = NULL;
	  if(cl.HasFlag("select")){
		doselect = true;
		fpselect = fopen("selected.times","w");
	  }

	  bool docosth = false;  // note: I all the varaibel costh, but it is really 1/cos(th)
	  double costhmincut=1.;    // 0 deg
	  double costhmaxcut=1.41;  // 45 deg
	  if(cl.HasFlag("docosth")){
		docosth = true;
		costhmincut = (double) abs(cl.GetArgument<double>("costhmin"));
		costhmaxcut = (double) abs(cl.GetArgument<double>("costhmax"));
	  }

	  bool poisson  = false;
	  if(cl.HasFlag("poisson")){
		poisson = true;
	  }

	  bool dodash =  false;
	  double dash_cut_value=0;
	  if(cl.HasFlag("dodash")){
		dodash = true;
		dash_cut_value = cl.GetArgument<double>("dash");
	  }

	  bool dolicompact =  false;
	  double licompact_cut_value=0;
	  if(cl.HasFlag("dolicompact")){
		dolicompact = true;
		licompact_cut_value = cl.GetArgument<double>("licompact");
	  }

	  bool doage =  false;
	  double age_cut_value=0;
	  if(cl.HasFlag("doage")){
		doage = true;
		age_cut_value = cl.GetArgument<double>("age");
	  }
	  int count = 0;

	  // loop over events in XCD file
	  while(f.Read()){
		Vector direction;
		direction.SetRThetaPhi(1.0,theta_x.At(0),phi_x.At(0));

		TimeStamp t((unsigned)gpsSec_x.At(0));
		ModifiedJulianDate mjd(t);
		alignment.Align(t,direction);

		const det::Detector& det  = detSv.GetDetector(t);
		const LatLonAlt& loc = det.GetLatitudeLongitudeHeight();

		
		// compute J2000 epoch RA,Dec
		EquPoint equ;
		astro.Loc2Equ(t,loc,direction,equ);
		EquPoint j2kPosn(equ);
		astro.PrecessFromEpochToJ2000(mjd,j2kPosn);

		// apply quality cuts
		bool quality = true;

		if(coreFiduScale_x.At(0) > fidu_cut_value)
		  quality = false;

		double costhval=1./cos(theta_x.At(0));
		if (costhval<costhmincut || costhval>costhmaxcut) 
		  quality = false;

		//printf("%d %d\n", fidu_cut_value, (int) coreFiduScale_x.At(0));

		if(angleFitStatus_x.At(0) != 0)
		  quality = false;

		if(nChAvail_x.At(0) < 700)
		  quality = false;

		if((double)nChAvail_x.At(0)/(double)nChTot_x.At(0) < 0.9)
		  quality = false;

		// apply gamma/hadron cuts
		bool isgamma = true;

		if(dopinc){
		  if(pinc_x.At(0) > pinc_cut_value)
		isgamma = false;
		}
		if(doproba){
		  if(proba_x.At(0) < proba_cut_value)
		isgamma = false;
		}

		if(doldfchi2){
		  if(ldfchi2_x.At(0) > ldfchi2_cut_value)
		isgamma = false;
		}
		if(docustom){
		  if(custom_x.At(0) < customcut)
		isgamma = false;
		}

		//if(dodash){
		//  if(dash_x.At(0) < dash_cut_value)
		//	isgamma = false;
		//}

		if(doage){
		  if(age_x.At(0) > age_cut_value)
		isgamma = false;
		}

		if(dolicompact){
		  double licompact = log10(cxpe40_x.At(0)/(double)nhitsp20_x.At(0));
		  if(licompact > licompact_cut_value )
		isgamma = false;
		}

		double racrab  = 83.6332083*(3.1415926536/180.);
		double deccrab = 22.0144444*(3.1415926536/180.);
		if(quality && isgamma){
		  count = count + 1;
		  sky->Fill(j2kPosn.GetRA()*180./3.1415926536,j2kPosn.GetDec()*180./3.1415926536);
		  allsky->Fill(j2kPosn.GetRA()*180./3.1415926536,j2kPosn.GetDec()*180./3.1415926536);
		  allsky2->Fill(j2kPosn.GetRA()*180./3.1415926536,j2kPosn.GetDec()*180./3.1415926536);
		  pointing pt(3.1415926536/2. - j2kPosn.GetDec(),j2kPosn.GetRA());
		  int pixel = skymap.ang2pix(pt);
		  skymap[pixel] = skymap[pixel] + 1;
		  double raev = j2kPosn.GetRA();
		  double deev = j2kPosn.GetDec();
		  double r = Delta(racrab,deccrab,raev,deev)*(180./3.1415926536);
		  if (doselect && r<0.25) { 
			printf("writing event: %d %llu %llu\n",(unsigned int)runid_x.At(0),
				   (unsigned long long int)gpsSec_x.At(0),(unsigned long long int)gpsNanosec_x.At(0));
			fprintf(fpselect,"%llu,%llu\n",(unsigned long long int)gpsSec_x.At(0),(unsigned long long int)gpsNanosec_x.At(0));
		  }
		}

	  } // end loop pber events

	  // make histogram of radial distribution from Crab Position
	  double racrab  = 83.6332083*(3.1415926536/180.);
	  double deccrab = 22.0144444*(3.1415926536/180.);
	  double bck     = 0.;
	  double bckarea = 0.;
	  for (int i=1;i<=sky->GetXaxis()->GetNbins();i++) {
		for (int j=1;j<=sky->GetYaxis()->GetNbins();j++) {
		  double rabin  = sky->GetXaxis()->GetBinCenter(i)*(3.1415926536/180.);
		  double decbin = sky->GetYaxis()->GetBinCenter(j)*(3.1415926536/180.);
		  double r = Delta(racrab,deccrab,rabin,decbin)*(180./3.1415926536);
		  double entries = sky->GetBinContent(i,j);
		  double binarea = sky->GetXaxis()->GetBinWidth(i)*sky->GetYaxis()->GetBinWidth(j)*cos(decbin);
		  evexcess->Fill(r,entries);
		  evrsq->Fill(r*r,entries);
		  evbins->Fill(r,1.*binarea);
		  //printf("%5.2f %5.2f %5.2f %5.2f %5.2f\n",racrab,deccrab,rabin,decbin,r);
		  // take background from region from 2-4 deg from crab
		  if (r>=2.0 && r<=4.0) {
			bck+=entries;
			bckarea+=binarea;
		  }
		}
	  }
	  double bckden = bck/bckarea;
	  
	  // open report file
	  printf("opening file: %s\n",report.c_str());
	  FILE *fp = fopen(report.c_str(),"a");
	  fprintf(fp,"-----\n");
	  fprintf(fp,"File                       : %s\n",recfile.c_str());
	  if (dopinc) 
		fprintf(fp,"PINC Cut                   : %6.2f\n",pinc_cut_value);
	  if(doproba)
		fprintf(fp,"proba Cut          : %6.6f\n",proba_cut_value);  
	  if(docustom)
		fprintf(fp,"custom name          : %s\n",customcutname.c_str()); 
		fprintf(fp,"custom Cut          : %6.6f\n",customcut); 
	  if (dolicompact) 
		fprintf(fp,"LICompactness              : %6.2f\n",licompact_cut_value);
	  if (doldfchi2) 
		fprintf(fp,"LDF Chi2 Cut               : %6.2f\n",ldfchi2_cut_value);
	  if (dofidu) 
		fprintf(fp,"Fiducial Core Cut          : %6d\n",fidu_cut_value);
	  if (docosth) 
		fprintf(fp,"1/cos(th) cut              : %6.3f - %6.3f\n",costhmincut,costhmaxcut);
	  fprintf(fp,"Bck (ev)                   : %8.2f\n",bck);
	  fprintf(fp,"Bck Area (deg sq)          : %8.2f\n",bckarea);
	  fprintf(fp,"Bck Density (bck ev/deg^2) : %8.2f\n",bckden);

	  // now plot gamma-ray density vs radius
	  double intevt = 0.;
	  double intbck = 0.;
	  for (int i=1;i<=evexcess->GetXaxis()->GetNbins();i++) {
		double ev = evexcess->GetBinContent(i);    // number of evts in the radial slice
		double bn = evbins->GetBinContent(i);      // number of bins in the radial slice
		if (bn>0.) {
		  evdensity->SetBinContent(i,ev/bn);
		  evdensity->SetBinError(i,(ev/bn)/sqrt(ev));
		  gadensity->SetBinContent(i,(ev/bn)-bckden);
		} else {
		  evdensity->SetBinContent(i,0);
		  gadensity->SetBinContent(i,0);
		}
		intevt += ev;
		intbck += bn*bckden;
		double z;
		if (poisson) {
		  z = sqrt(2.*( intevt*log(intevt/intbck) - (intevt-intbck) )); // poisson z
		} else {
		  z = (intevt-intbck)/sqrt(intbck); // excess/sqrt(bkg) 
		}
		//printf("lima: %f %f %f %f %f\n",intevt,bck,bn/bckarea,bck*bn/bckarea,z);
		evint->SetBinContent(i,intevt);
		bckint->SetBinContent(i,intbck);
		exint->SetBinContent(i,intevt-intbck);
		evsig->SetBinContent(i,z);
		//printf("%f %f %f %f %f\n",ev,bn,intevt,intbck,z);
	  }
	  // find optimal bin size
	  double maxsig = evsig->GetMaximum();
	  int    ioptcut = evsig->GetMaximumBin();
	  double optcut = evsig->GetXaxis()->GetBinUpEdge(ioptcut);
	  double optexcess = exint->GetBinContent(ioptcut);
	  double optbck    = bckint->GetBinContent(ioptcut);
	  double totexcess = exint->GetBinContent(ioptcut*3);     // excess at 3x opt cut
	  int    iminoptcut = 0;
	  double minoptcut = 100.; 
	  double minoptcutexcess = 0.;
	  int    imaxoptcut = 0;
	  double maxoptcut = 0.;
	  double maxoptcutexcess = 0.;
	  double frac68 = 0.;              // fraction in bin closest to 68%
	  double signif68 = 0.;            // signif when using bin that contains 68% of excess
	  double cut68 = 0.;               // bin size which give 68% containment
	  for (int i=1;i<=evsig->GetXaxis()->GetNbins();i++) {
		double binedge = evsig->GetXaxis()->GetBinUpEdge(i);
		if (evsig->GetBinContent(i)>0.9*maxsig) {
		  if (binedge<minoptcut) { minoptcut=binedge; iminoptcut=i; minoptcutexcess=exint->GetBinContent(i); }
		  if (binedge>maxoptcut) { maxoptcut=binedge; imaxoptcut=i; maxoptcutexcess=exint->GetBinContent(i); }
		}
		double frac = exint->GetBinContent(i)/totexcess;
		if (fabs(frac68-0.68)>fabs(frac-0.68) && binedge<1.8) {
		  frac68   = frac;
		  signif68 = evsig->GetBinContent(i);
		  cut68    = binedge;
		}
	  }
	  fprintf(fp,"Peak Significance (sigmas) : %6.2f",maxsig);
	  if (poisson) {
		fprintf(fp," (Pois)\n");
	  } else {
		fprintf(fp," (Gaus)\n");
	  }
	  fprintf(fp,"Optimal Bin Size (deg)     : %6.3f\n",optcut);
	  fprintf(fp,"Optimal Bin Excess (ev)    : %8.2f\n",optexcess);
	  fprintf(fp,"Optimal Bin Bkg (ev)       : %8.2f\n",optbck);
	  fprintf(fp,"Optimal Bin S/B Ratio      : %8.2f\n",optexcess/optbck);
	  fprintf(fp,"Optimal Bin Fraction       : %6.2f\n",100.*optexcess/totexcess);
	  fprintf(fp,"Optimal Bin Range (90 pct) : %6.2f - %6.2f\n",minoptcut,maxoptcut);
	  fprintf(fp,"Optimal Bin Range Fraction : %6.2f - %6.2f\n",100*minoptcutexcess/totexcess,100*maxoptcutexcess/totexcess);
	  fprintf(fp,"68 Pct Bin                 : %6.3f\n",cut68);
	  fprintf(fp,"Fraction in 68 Pct Bin     : %6.2f\n",frac68*100);
	  fprintf(fp,"Significance of 68 Pct Bin : %6.2f",signif68);
	  if (poisson) {
		fprintf(fp," (Pois)\n");
	  } else {
		fprintf(fp," (Gaus)\n");
	  }
	  fprintf(fp,"-----\n");
	  fclose(fp);
	  if (doselect) fclose(fpselect);

	  if(cl.HasFlag("healpix")){
		if(!fexists(output.c_str()))
		{
		arr<std::string> colname(1);
		colname[0] = "sensitivity";
		fitshandle out;
		out.create(output.c_str());
		prepare_Healpix_fitsmap(out,skymap, PLANCK_FLOAT64, colname);
		out.write_column(1,skymap.Map());
		out.close();
		}
	  } else {
		if(!fexists(output.c_str())){
		evrsq->GetXaxis()->SetTitle("R sq (deg sq)");
		
		TFile* out = new TFile(output.c_str(),"RECREATE");
		sky->Write();
		allsky->Write();
		allsky2->Write();
		evexcess->Write();
		evbins->Write();
		evdensity->Write();
		gadensity->Write();
		evint->Write();
		evsig->Write();
		evrsq->Write();
		out->Write();
		}
	  }
	}
}
