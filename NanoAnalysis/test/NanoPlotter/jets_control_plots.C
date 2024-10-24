#include "CMS_lumi.C"
#include <tuple>
#include <vector>

void DrawRatioPlot(string name, TCanvas *c, TH1D *data, TH1D *MC, TH1D*MCUp, TH1D*MCDn, double _lumi){

        c->Divide(0,2,0,0);
        // TPad * pad1 = new TPad("pad1","pad1", 0, 0.3, 1, 1.0);
        // pad1->Draw();
        c->cd(1);

        gPad->SetBottomMargin(0.02);
        gPad->SetTopMargin(0.18);
        gPad->SetLeftMargin(0.10);

        // MC->Scale( data->Integral() / MC->Integral() );

        double max = 0.;
        for (int bin = 0; bin < MC->GetSize() - 2; bin++){
          if ( MC->GetBinContent(bin + 1) > max) max = MC->GetBinContent(bin + 1);
        }
        MC->SetMaximum(1.4*max);

        MC->SetFillColor(kOrange + 1);
        MC->GetYaxis()->SetTitle("Events/0.2");
        MC->Draw("HIST");
        MC->GetXaxis()->SetLabelSize(0);
        MC->GetYaxis()->SetTitleSize(0.07);
        MC->GetYaxis()->SetLabelSize(0.07);
        MC->GetYaxis()->SetTitleOffset(0.7);

        data->SetLineColor(kBlack);
        data->SetMarkerStyle(20);
        data->SetMarkerSize(0.6);
        data->Draw("p E1 X0 SAME");

        // c->Update();
        //
        // c->cd();

        // TPad * pad2 = new TPad("pad2","pad2", 0, 0.05, 1, 0.3);
        // pad2->SetGridy();
        // pad2->Draw();
        c->cd(2);

        gPad->SetBottomMargin(0.20);
        gPad->SetTopMargin(0.02);
        gPad->SetLeftMargin(0.10);

        TH1F * Ratio = new TH1F("Ratio","Ratio", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());
        TH1F * sigma_up = new TH1F("sigma_up","sigma_up", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());
        TH1F * sigma_dn = new TH1F("sigma_dn","sigma_dn", MC->GetSize() - 2, MC->GetXaxis()->GetXmin(), MC->GetXaxis()->GetXmax());

        for(int bin = 0; bin < MC->GetSize() - 2; bin++){
          if (MC->GetBinContent(bin + 1) == 0){
                  Ratio   ->SetBinContent(bin + 1, 0);
                  Ratio   ->SetBinError(bin + 1, 0);
                  sigma_up->SetBinContent(bin + 1, 0);
                  sigma_dn->SetBinContent(bin + 1, 0);
          }else{
                  Ratio   ->SetBinContent(bin + 1, float(data->GetBinContent(bin + 1))/MC->GetBinContent(bin + 1));
                  Ratio   ->SetBinError  (bin + 1, 1./pow((data->GetBinContent(bin + 1)),2));
                  sigma_up->SetBinContent(bin + 1, float(MCUp->GetBinContent(bin + 1))/MC->GetBinContent(bin + 1));
                  sigma_dn->SetBinContent(bin + 1, float(MCDn->GetBinContent(bin + 1))/MC->GetBinContent(bin + 1));
          }
        }

        double temp;
        for(int bin = 0; bin < Ratio->GetSize() - 2; bin++){
          temp = Ratio->GetBinContent(bin + 1);
          Ratio->SetBinContent( bin + 1, temp - 1);
        }

        for(int bin = 0; bin < sigma_up->GetSize() - 2; bin++){
          temp = sigma_up->GetBinContent(bin + 1);
          sigma_up->SetBinContent( bin + 1, temp - 1);
        }

        for(int bin = 0; bin < sigma_dn->GetSize() - 2; bin++){
          temp = sigma_dn->GetBinContent(bin + 1);
          sigma_dn->SetBinContent( bin + 1, temp - 1);
        }

        sigma_up->SetFillColor(kGray);
        sigma_up->SetLineColor(kGray);
        sigma_up->SetMaximum(2.0);
        sigma_up->SetMinimum(-2.0);
        sigma_up->GetXaxis()->SetTitle("Leading jet #eta");
        sigma_up->GetYaxis()->SetTitle("(Data/MC)-1");
        sigma_up->GetYaxis()->SetTitleOffset(0.6);
        sigma_up->Draw("HIST");
        sigma_up->GetXaxis()->SetLabelSize(0.07);
        sigma_up->GetXaxis()->SetTitleSize(0.07);
        sigma_up->GetYaxis()->SetLabelSize(0.07);
        sigma_up->GetYaxis()->SetTitleSize(0.07);

        sigma_dn->SetFillColor(kGray);
        sigma_dn->SetLineColor(kGray);
        sigma_dn->Draw("HIST SAME");

        gPad->RedrawAxis();

        Ratio->SetMarkerStyle(20);
        Ratio->SetMarkerSize(0.6);
        Ratio->Draw("p E1 X0 SAME");

        c->cd(1);
        TLegend * leg = new TLegend(0.74,0.45,0.95,0.75);
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.07);
        leg->AddEntry(MC, "DY + t#bar{t}  MC", "f" );
        leg->AddEntry(data, "Data", "p");
        leg->AddEntry(sigma_up, "JEC Uncertainty", "f");
        leg->Draw();

        c->Update();

        c->SaveAs((TString)name+".pdf");
        c->SaveAs((TString)name+".png");
        c->SaveAs((TString)name+".C");

        c->Clear();
}


double computeLeadJetUp(ROOT::RVec<float> &_JetEta, ROOT::RVec<float> &_JetPt_JER, ROOT::RVec<float> &_JetSigma){
  double highest = 0;
  int i_up = 0;
  double variation = 0;
  for(int i=0; i<_JetEta.size(); i++){
    variation = _JetPt_JER.at(i)*(1+_JetSigma.at(i));
    if (variation > highest && variation > 30){
            highest = variation;
            i_up = i;
    }
  }
  return _JetEta.at(i_up);

}

double computeLeadJetDn(ROOT::RVec<float> &_JetEta, ROOT::RVec<float> &_JetPt_JER, ROOT::RVec<float> &_JetSigma){
  double highest = 0;
  int i_up = 0;
  double variation = 0;
  for(int i=0; i<_JetEta.size(); i++){
    variation = _JetPt_JER.at(i)*(1-_JetSigma.at(i));
    if (variation > highest && variation > 30){
            highest = variation;
            i_up = i;
    }
  }
  return _JetEta.at(i_up);
}




tuple<ROOT::RDF::RResultPtr<TH1D>,ROOT::RDF::RResultPtr<TH1D>,ROOT::RDF::RResultPtr<TH1D>> Histo_MC(string fpath, double _lumi){

  ROOT::RDataFrame rdf ("ZTree/candTree", fpath);
  TFile * f             = new TFile((TString)fpath,"READ");
  TH1F  * t1            = (TH1F*)f->Get("ZTree/Counters");
  double gen_sumWeights = t1->GetBinContent(1);
  auto skim_rdf         = rdf.Filter("Zsel>0 && JetEta.size()>0")
                             .Define("weight", Form("1000 * xsec * overallEventWeight * %f / %f", _lumi, gen_sumWeights))
                             .Define("leadEta", "JetEta.at(0)")
                             .Define("leadEtaUp",   computeLeadJetUp, {"JetEta", "JetPt_JERUp",   "JetSigma"})
                             .Define("leadEtaDown", computeLeadJetDn, {"JetEta", "JetPt_JERDown", "JetSigma"});

  auto hist_nominal = skim_rdf.Histo1D({"leadEta_MC",     "leadEta_MC",     47, -4.7, 4.7}, "leadEta",     "weight");
  auto hist_up      = skim_rdf.Histo1D({"leadEtaUp_MC",   "leadEtaUp_MC",   47, -4.7, 4.7}, "leadEtaUp",   "weight");
  auto hist_dn      = skim_rdf.Histo1D({"leadEtaDown_MC", "leadEtaDown_MC", 47, -4.7, 4.7}, "leadEtaDown", "weight");

  return std::make_tuple(hist_nominal, hist_up, hist_dn);

}




void makePlot(){

  string fDY   = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2017/DYJetsToLL_M50/ZZ4lAnalysis.root";
  string fTT   = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2017/TTTo2L2Nu/ZZ4lAnalysis.root";
  string fdata = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200430_LegacyRun2/Data_2017/AllData/ZZ4lAnalysis.root";
  double lumi  = 41.5;

  ROOT::RDF::RResultPtr<TH1D> hist_DY;
  ROOT::RDF::RResultPtr<TH1D> hist_DY_up;
  ROOT::RDF::RResultPtr<TH1D> hist_DY_dn;
  tie(hist_DY,hist_DY_up,hist_DY_dn) = Histo_MC(fDY, lumi);

  ROOT::RDF::RResultPtr<TH1D> hist_ttbar;
  ROOT::RDF::RResultPtr<TH1D> hist_ttbar_up;
  ROOT::RDF::RResultPtr<TH1D> hist_ttbar_dn;
  tie(hist_ttbar,hist_ttbar_up,hist_ttbar_dn) = Histo_MC(fTT, lumi);


  ROOT::RDataFrame rdf_data("ZTree/candTree", fdata);
  auto hist_data  = rdf_data.Filter("Zsel>0 && JetEta.size()>0")
                            .Define("leadEta", "JetEta.at(0)")
                            .Histo1D({"leadEta_data", "leadEta_data", 47, -4.7, 4.7}, "leadEta");

  hist_DY    -> Add(hist_ttbar.GetPtr());
  hist_DY_up -> Add(hist_ttbar_up.GetPtr());
  hist_DY_dn -> Add(hist_ttbar_dn.GetPtr());

  std::cout << "Saving distributions into root file ..." << std::endl;
  TFile * outfile = new TFile("test_histos.root", "RECREATE");
  outfile    -> cd();
  hist_data  -> Write();
  hist_DY    -> Write();
  hist_DY_up -> Write();
  hist_DY_dn -> Write();
  outfile    -> Close();

  TCanvas * canvas = new TCanvas();
  DrawRatioPlot("test", canvas,
                hist_data.GetPtr(),
                hist_DY.GetPtr(), hist_DY_up.GetPtr(), hist_DY_dn.GetPtr(),
                lumi);

}