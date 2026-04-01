// C++
#include <iostream>
#include <fstream>
#include <string>

// ROOT
#include "TApplication.h"
#include <TROOT.h>
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"

// My own files
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/include/OSmethod.h>
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/src/setTDRStyle.cpp>


using namespace std;

int main( int argc, char *argv[] )
{
   setTDRStyle();

   TString YEAR = "2024";
   
   TString Data;
   if ( YEAR == "2022" ) { Data = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraCD_preEE_SKIMMED.root"; }
   if ( YEAR == "2022EE" ) { Data = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2022_Data/Data_eraEFG_postEE_SKIMMED.root"; }
   if ( YEAR == "2023preBPix" )  { Data = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraC_preBPix_SKIMMED.root"; }
   if ( YEAR == "2023postBPix" )  { Data = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2023_Data/Data_eraD_postBPix_SKIMMED.root"; }
   if ( YEAR == "2024" )  { Data = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/2024_Data/ZZ4lAnalysis_SKIMMED.root"; }

   TString path = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/" + YEAR + "_MC";
   
   TString file_name = "/ZZ4lAnalysis_SKIMMED.root";
   
   TString WZ      = path + "/WZto3LNu"       + file_name;
   TString ZZ      = path + "/ZZTo4l"        + file_name;
   TString ttbar   = path + "/TTto2L2Nu"      + file_name;
   TString DY      = path + "/DYJetsToLL"     + file_name;

   TString DataMCOSFile      = "DataMC_OS_" + YEAR + ".root";
   TString HistosOSFile      = "Histos_OS_" + YEAR + ".root";
   TString FakeRatesOSFile   = "FakeRates_OS_" + YEAR + ".root";
   TString ZXHistosOSFile    = "ZXHistos_OS_" + YEAR + ".root";

   bool SubtractWZ = true;
   bool Remove_NegBins_FR = true;
   bool Remove_NegBins_ZX = true;

   float pT_bins[] = {5, 7, 10, 20, 30, 40, 50, 80};

   OSmethod *os = new OSmethod();
   //os->SetLumi(35.92); //2016
   //os->SetLumi(41.53);   //2017
   //os->SetLumi(59.74); //2018
   if ( YEAR == "2022" ) { os->SetLumi(7.9804);}
   if ( YEAR == "2022EE" ) { os->SetLumi(26.6728);}
   if ( YEAR == "2023preBPix" ) { os->SetLumi(18.06265);}
   if ( YEAR == "2023postBPix" ) { os->SetLumi(9.693);}
   if ( YEAR == "2024" ) { os->SetLumi(108.822);}

   ///////////////////////////////////
   // Fill control histos           //
   ///////////////////////////////////
   cout << "** Fill control histos" << endl;
   os->FillDataMCPlots(Data);
   os->FillDataMCPlots(WZ);
   os->FillDataMCPlots(ZZ);
   os->FillDataMCPlots(ttbar);
   os->FillDataMCPlots(DY);
   os->SaveDataMCHistos(DataMCOSFile);
   
   ///////////////////////////////////
   // Fill passing/failing histos  //
   ///////////////////////////////////
   cout << "** Fill passing/failing histos" << endl;
   os->FillFRHistos(Data);
   os->FillFRHistos(WZ);
   os->SaveFRHistos(HistosOSFile, SubtractWZ, Remove_NegBins_FR);

   ///////////////////////////////////
   // Calculate fake rates          //
   ///////////////////////////////////
   cout << "** Calculate fake rates" << endl;
   os->GetFRHistos(HistosOSFile);
   os->Set_pT_binning(8, pT_bins);
   os->ProduceFakeRates(FakeRatesOSFile);

   ///////////////////////////////////
   // Fill ZX contributions histos  //
   ///////////////////////////////////
   cout << "** Fill ZX contributions histos" << endl;
   os->MakeHistogramsZX(Data, FakeRatesOSFile);
   os->MakeZXMCContribution(ZZ, FakeRatesOSFile);
   os->SaveZXHistos(ZXHistosOSFile, Remove_NegBins_ZX);

   ///////////////////////////////////
   // Plot control plots            //
   ///////////////////////////////////
   cout << "** Plot control plots" << endl;
   os->GetZXHistos(ZXHistosOSFile);
   os->PrintZXYields();
   os->GetDataMCHistos(DataMCOSFile);
   os->PlotDataMC("M4l", "Plots");
   os->PlotDataMC_2P2F( "M4l", "Plots" );
   os->PlotDataMC_3P1F( "M4l", "Plots" );

   ///////////////////////////////////
   // Plot Z+X plots                //
   ///////////////////////////////////
   cout << "** Plot Z+X plots" << endl;
   os->GetZXHistos(ZXHistosOSFile);
   os->PlotZXContributions("Plots");
   os->FitZX("Plots");
	
   cout << "** DONE" << endl;

   delete os;
}
