#include <iostream>
#include <TMath.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <cmath>
#include <TStyle.h>
#include <TLegend.h>


void selection() {
    TFile ifile("newroot.root", "read");
    TTree* datatree = (TTree*)ifile.Get("h10");
    int nph; float eph[45]; float thetaph[45]; float phiph[45];

    datatree->SetBranchAddress("nph", &nph);
    datatree->SetBranchAddress("eph", &eph);
    datatree->SetBranchAddress("thetaph", thetaph);
    datatree->SetBranchAddress("phiph", phiph);

    std::vector<int> candidates;
    TH1D* hist_inv_mass = new TH1D("hist_inv_mass","inv mass, GeV", 50, 0.09, 0.21);
    TH1D* hist_angle = new TH1D("hist_angle", "angle", 50, 0, TMath::Pi());

    for (int i=0; i<datatree->GetEntries(); i++) {
        datatree->GetEntry(i);
        int candidate_num = 0;
        std::vector<double> inv_mass_vec;
        for (int j = 0; j < nph - 1; j++) {
            for (int k = j+1; k < nph; k++) {
                TVector3 v1(sin(thetaph[j])*cos(phiph[j]), sin(thetaph[j])*sin(phiph[j]), cos(thetaph[j]));
                TVector3 v2(sin(thetaph[k])*cos(phiph[k]), sin(thetaph[k])*sin(phiph[k]), cos(thetaph[k]));
                auto ang = v1.Angle(v2);
                double m_inv = sqrt(2. * eph[j] * eph[k] * (1. - cos(ang)));
                hist_angle->Fill(ang);
                if (0.1 <= m_inv && m_inv <= 0.2){
                    candidate_num++;
                    inv_mass_vec.push_back(m_inv);
                }
            }
        }
        candidates.push_back(candidate_num);
        if (candidate_num == 2){
            for (double val : inv_mass_vec){
                hist_inv_mass->Fill(val);
            }
        }
    }
    
    TCanvas *canvas = new TCanvas("canvas", "title.png", 1200, 900);

    canvas->Divide(1, 2);

    canvas->cd(1);
    hist_inv_mass->SetTitle("hist_inv_mass");
    hist_inv_mass->SetName("hist_inv_mass"); 
    hist_inv_mass->SetLineWidth(2);
    hist_inv_mass->SetXTitle("inv mass, MeV");
    hist_inv_mass->SetYTitle("N");
    hist_inv_mass->Draw();

    canvas->cd(2);
    hist_angle->SetTitle("hist_angle");
    hist_angle->SetName("hist_angle"); 
    hist_angle->SetLineWidth(2);
    hist_angle->SetXTitle("angle");
    hist_angle->SetYTitle("N");
    hist_angle->Draw();

    TFile* ofile = new TFile("selected.root", "recreate");
    hist_inv_mass->Write();
    hist_angle->Write();
}