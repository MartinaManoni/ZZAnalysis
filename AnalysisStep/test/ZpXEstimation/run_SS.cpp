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
#include <ZZAnalysis/AnalysisStep/test/ZpXEstimation/include/SSmethod.h>
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

   TString path = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/HIG-25-015/RunIII_byZ1Z2/Moriond26_JES/" + YEAR +"_MC";
   
   TString file_name  = "/ZZ4lAnalysis_SKIMMED.root";
   
   TString WZ      = path + "/WZto3LNu"       + file_name;
   TString ZZ      = path + "/ZZTo4l"        + file_name;
   TString ttbar   = path + "/TTto2L2Nu"      + file_name;
   TString DY      = path + "/DYJetsToLL"     + file_name;

   TString DataMCSSFile      = "DataMC_SS_" + YEAR + ".root";
   TString HistosSSFile      = "Histos_SS_" + YEAR + ".root";
   TString FakeRatesSSFile   = "FakeRates_SS_" + YEAR + ".root";
   TString ZXHistosSSFile    = "ZXHistos_SS_" + YEAR + ".root";

   bool SubtractWZ = true;
   bool Remove_NegBins_FR = true;
   bool SubtractMCContribution = true;
	
   float pT_bins[] = {5, 7, 10, 20, 30, 40, 50, 80};
   SSmethod *ss = new SSmethod();
   //ss->SetLumi(35.92); // 2016 lumi
   //ss->SetLumi(41.53); // 2017 lumi
   //ss->SetLumi(59.74); // 2018 lumi
   if ( YEAR == "2022" ) { ss->SetLumi(7.9804);}
   if ( YEAR == "2022EE" ) { ss->SetLumi(26.6728);}
   if ( YEAR == "2023preBPix" ) { ss->SetLumi(18.06265);}
   if ( YEAR == "2023postBPix" ) { ss->SetLumi(9.693);}
   if ( YEAR == "2024" ) { ss->SetLumi(108.822);}

   ///////////////////////////////////
   // Fill control histos           //
   ///////////////////////////////////
   cout << "** Fill control histos" << endl;
   ss->FillDataMCPlots(Data);
   ss->FillDataMCPlots(WZ);
   ss->FillDataMCPlots(ZZ);
   ss->FillDataMCPlots(ttbar);
   ss->FillDataMCPlots(DY);
   ss->SaveDataMCHistos(DataMCSSFile);
   
   ///////////////////////////////////
   // Fill passing/failling histos  //
   ///////////////////////////////////
   ss->FillFRHistos(Data);
   ss->FillFRHistos(WZ);
   ss->SaveFRHistos(HistosSSFile, SubtractWZ, Remove_NegBins_FR);

   ///////////////////////////////////
   // Calculate fake rates          //
   ///////////////////////////////////
   ss->GetFRHistos(HistosSSFile);
   ss->Set_pT_binning(8, pT_bins);
   ss->ProduceFakeRates(FakeRatesSSFile, Data);

   ///////////////////////////////////
   // Calculate OS/SS ratios        //
   ///////////////////////////////////
   ss->Calculate_SSOS_Ratio( Data, ZZ, SubtractMCContribution);

   ///////////////////////////////////
   // Fill ZX contributions histos  //
   ///////////////////////////////////
   ss->MakeHistogramsZX(Data, FakeRatesSSFile);
   ss->SaveZXHistos(ZXHistosSSFile);

   ///////////////////////////////////
   // Plot control plots            //
   ///////////////////////////////////
   ss->GetDataMCHistos(DataMCSSFile);
   ss->PlotDataMC( "M4l", "Plots" );

   ///////////////////////////////////
   // Plot and fit Z+X              //
   ///////////////////////////////////
   ss->GetZXHistos(ZXHistosSSFile);
   ss->PlotZX("M4l", "Plots");
   ss->FitZX("M4l", "Plots");

   delete ss;
}
