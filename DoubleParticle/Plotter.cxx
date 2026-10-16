// ----------------------------------------------------------------------------
// Script for plotting results from DoubleParticleBenchmark.cxx
// (e + gamma, Expected / Remnant clusters).
//
// For each region ONE multi-page PDF is created (no single-plot PDFs).
//
// Run:     root -l -b -q PlottingDoubleParticle.cxx
// Input:   Plots/<label>.root
// Output:  Plots/<label>/All_<label>.pdf
//
// All settings that need editing are in the "SETTINGS" section below.
// ----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TEfficiency.h>

using namespace std;

// ============================== SETTINGS ====================================
struct Sample {
  string label;       // file: Plots/<label>.root, directory: Plots/<label>/
  string title;       // title displayed on plots
  double eta[2];      // eta range
};

const vector<Sample> samples = {
  // label                  title                              eta
  {"ElectronPhotonBarrel",  "e^{-}+#gamma Barrel ECal",        {-1.8,  1.5}},
  {"ElectronPhotonNEndCap", "e^{-}+#gamma N Endcap ECal",      {-3.6, -1.7}},
  {"ElectronPhotonPEndCap", "e^{-}+#gamma P Endcap ECal",      { 1.3,  3.6}},
};

const int    kColElec   = kAzure + 1;   // electron (Expected clusters)
const int    kColPhot   = kRed + 1;     // photon (Remnant clusters)
const int    kColTot    = kGreen + 2;   // sum
const double kTextSize  = 0.042;
const bool   kHeader    = true;         // "ePIC Benchmark | ..." header (false for a paper)
const int    kCanvasW   = 800;
const int    kCanvasH   = 700;
// ============================================================================


// ---------- helper functions -------------------------------------------------
vector<TCanvas*> gCanvases;   // all pages of the current region

// Print all pages into ONE multi-page PDF, then clean up.
void FlushCombined(const string& file) {
  if (gCanvases.empty()) return;
  gCanvases.front()->Print((file + "[").c_str());
  for (auto c : gCanvases) c->Print(file.c_str());
  gCanvases.back()->Print((file + "]").c_str());
  for (auto c : gCanvases) delete c;
  gCanvases.clear();
}

void SetStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillStyle(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.3, "Y");
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(64);
}

string Fmt(const char* fmt, double a, double b = 0) { return string(Form(fmt, a, b)); }

template <typename T>
T* Get(TFile* f, const char* name) {
  TObject* o = f->Get(name);
  if (!o) cerr << "BRAK obiektu: " << name << " w " << f->GetName() << endl;
  return dynamic_cast<T*>(o);
}

TCanvas* NewCanvas(const Sample& s, const string& name, bool is2D = false) {
  auto c = new TCanvas((name + "_" + s.label).c_str(), "", kCanvasW, kCanvasH);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.14);
  c->SetTopMargin(kHeader ? 0.07 : 0.03);
  c->SetRightMargin(is2D ? 0.15 : 0.05);
  return c;
}

// Only collects the page - the PDF is written once in FlushCombined().
void Save(TCanvas* c) { gCanvases.push_back(c); }

void Header(const string& title) {
  if (!kHeader) return;
  TLatex t;
  t.SetNDC();
  t.SetTextFont(62);
  t.SetTextSize(0.045);
  t.DrawLatex(gPad->GetLeftMargin(), 0.94, ("ePIC Benchmark | " + title).c_str());
}

void Text(const vector<string>& lines, double x = 0.58, double y = 0.88) {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(kTextSize);
  for (auto& l : lines) { t.DrawLatex(x, y, l.c_str()); y -= 1.3 * kTextSize; }
}

void Style1D(TH1* h, int col) {
  h->SetTitle("");
  h->SetLineColor(col);
  h->SetLineWidth(2);
  h->SetFillColorAlpha(col, 0.3);
}

void VLine(double x, TH1* h) {
  auto l = new TLine(x, 0, x, h->GetMaximum());
  l->SetLineStyle(2);
  l->SetLineColor(kGray + 2);
  l->Draw();
}


// ---------- plot types -------------------------------------------------------

// Standard 1D histogram. stats -> Mean/Std/Entries; vref >= 0 -> dashed vertical line.
void DrawHist(const Sample& s, TH1* h, int col, const char* xt, const char* yt,
              const string& name, bool stats = true, double vref = -1) {
  if (!h) return;
  auto c = NewCanvas(s, name);
  Style1D(h, col);
  h->GetXaxis()->SetTitle(xt);
  h->GetYaxis()->SetTitle(yt);
  h->SetMaximum(1.35 * h->GetMaximum());
  h->Draw("hist");
  if (vref >= 0) VLine(vref, h);
  Header(s.title);
  if (stats)
    Text({Fmt("Mean = %.3f", h->GetMean()),
          Fmt("Std = %.3f",  h->GetStdDev()),
          Fmt("Entries = %.0f", h->GetEntries())});
  Save(c);
}

// Two 1D histograms on one plot (electron vs photon).
void DrawOverlay(const Sample& s, TH1* a, const char* la, TH1* b, const char* lb,
                 const char* xt, const char* yt, const string& name, double vref = -1) {
  if (!a || !b) return;
  auto c = NewCanvas(s, name);
  Style1D(a, kColElec);
  Style1D(b, kColPhot);
  a->GetXaxis()->SetTitle(xt);
  a->GetYaxis()->SetTitle(yt);
  a->SetMaximum(1.45 * std::max(a->GetMaximum(), b->GetMaximum()));
  a->Draw("hist");
  b->Draw("hist same");
  if (vref >= 0) VLine(vref, a);
  Header(s.title);
  auto leg = new TLegend(0.58, 0.76, 0.94, 0.90);
  leg->SetTextSize(kTextSize);
  leg->AddEntry(a, Form("%s (#mu=%.3f, #sigma=%.3f)", la, a->GetMean(), a->GetStdDev()), "f");
  leg->AddEntry(b, Form("%s (#mu=%.3f, #sigma=%.3f)", lb, b->GetMean(), b->GetStdDev()), "f");
  leg->Draw();
  Save(c);
}

// 2D map, optionally with a profile and/or a diagonal y = x.
void Draw2D(const Sample& s, TH2* h, TProfile* prof, const char* xt, const char* yt,
            const string& name, bool diagonal = false) {
  if (!h) return;
  auto c = NewCanvas(s, name, true);
  h->SetTitle("");
  h->GetXaxis()->SetTitle(xt);
  h->GetYaxis()->SetTitle(yt);
  h->Draw("COLZ");
  if (diagonal) {
    double lo = std::max(h->GetXaxis()->GetXmin(), h->GetYaxis()->GetXmin());
    double hi = std::min(h->GetXaxis()->GetXmax(), h->GetYaxis()->GetXmax());
    auto l = new TLine(lo, lo, hi, hi);
    l->SetLineStyle(2); l->SetLineColor(kBlack); l->Draw();
  }
  if (prof) {
    prof->SetLineColor(kBlack);  prof->SetMarkerColor(kBlack);
    prof->SetMarkerStyle(20);    prof->SetMarkerSize(0.9); prof->SetLineWidth(2);
    prof->Draw("same");
  }
  Header(s.title);
  if (prof) Text({"Black: mean profile"}, 0.50, 0.88);
  Save(c);
}

// Profiles (up to 3) on one plot.
void DrawProfiles(const Sample& s, vector<TProfile*> p, vector<int> cols, vector<string> labels,
                  const char* xt, const char* yt, const string& name,
                  double ymin, double ymax, double yref = -1) {
  auto c = NewCanvas(s, name);
  bool first = true;
  auto leg = new TLegend(0.50, 0.76, 0.94, 0.90);
  leg->SetTextSize(kTextSize);
  for (size_t i = 0; i < p.size(); ++i) {
    if (!p[i]) continue;
    p[i]->SetTitle("");
    p[i]->SetLineColor(cols[i]); p[i]->SetMarkerColor(cols[i]);
    p[i]->SetMarkerStyle(20 + i); p[i]->SetMarkerSize(1.0); p[i]->SetLineWidth(2);
    p[i]->GetXaxis()->SetTitle(xt);
    p[i]->GetYaxis()->SetTitle(yt);
    p[i]->SetMinimum(ymin);
    p[i]->SetMaximum(ymax);
    p[i]->Draw(first ? "E1" : "E1 same");
    leg->AddEntry(p[i], labels[i].c_str(), "lp");
    first = false;
  }
  if (first) { delete c; return; }
  if (yref >= 0) {
    auto l = new TLine(p[0] ? p[0]->GetXaxis()->GetXmin() : 0, yref,
                       p[0] ? p[0]->GetXaxis()->GetXmax() : 1, yref);
    l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw();
  }
  leg->Draw();
  Header(s.title);
  Save(c);
}

// Efficiencies (up to 2) on one plot.
void DrawEff(const Sample& s, vector<TEfficiency*> e, vector<int> cols, vector<string> labels,
             const char* xt, const string& name, double xlo, double xhi) {
  auto c = NewCanvas(s, name);
  auto fr = c->DrawFrame(xlo, 0.0, xhi, 1.25);
  fr->GetXaxis()->SetTitle(xt);
  fr->GetYaxis()->SetTitle("Efficiency");
  auto leg = new TLegend(0.45, 0.76, 0.94, 0.90);
  leg->SetTextSize(kTextSize);
  bool any = false;
  for (size_t i = 0; i < e.size(); ++i) {
    if (!e[i]) continue;
    e[i]->SetLineColor(cols[i]); e[i]->SetMarkerColor(cols[i]);
    e[i]->SetMarkerStyle(20 + i); e[i]->SetMarkerSize(1.0); e[i]->SetLineWidth(2);
    e[i]->Draw("P same");
    leg->AddEntry(e[i], labels[i].c_str(), "lp");
    any = true;
  }
  if (!any) { delete c; return; }
  auto l = new TLine(xlo, 1.0, xhi, 1.0);
  l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw();
  leg->Draw();
  Header(s.title);
  Save(c);
}


// ---------- all plots for one sample ----------------------------------------
void DrawSample(const Sample& s, TFile* f) {
  gSystem->mkdir(("Plots/" + s.label).c_str(), true);
  gCanvases.clear();
  const string pdf = "Plots/All_" + s.label + ".pdf";

  // --- load ---
  auto hEMCe = Get<TH1D>(f, "h_EMC_Elec");   auto hEMCg = Get<TH1D>(f, "h_EMC_Phot");
  auto hDRMC = Get<TH1D>(f, "h_DeltaR_MC_e_g");
  auto hNExp = Get<TH1D>(f, "h_NExpected");  auto hNRem = Get<TH1D>(f, "h_NRemnant");
  auto hEExp = Get<TH1D>(f, "h_EExp");       auto hERem = Get<TH1D>(f, "h_ERem");
  auto hERemL = Get<TH1D>(f, "h_ERemLead");
  auto rE  = Get<TH1D>(f, "h_Ratio_Elec");   auto rG  = Get<TH1D>(f, "h_Ratio_Phot");
  auto rGL = Get<TH1D>(f, "h_Ratio_PhotLead"); auto rT = Get<TH1D>(f, "h_Ratio_Tot");
  auto eOp = Get<TH1D>(f, "h_ExpOverPTrack");
  auto pOm = Get<TH1D>(f, "h_PTrackOverMC");
  auto rOp = Get<TH1D>(f, "h_RemOverPTrack");
  auto m2e = Get<TH2D>(f, "h2_EExp_vs_EMC"); auto m2g = Get<TH2D>(f, "h2_ERem_vs_EMC");
  auto hRRE = Get<TH1D>(f, "h_RRemExp");     auto hDEP = Get<TH2D>(f, "h_DeltaEtaPhi");
  auto pE  = Get<TProfile>(f, "p_RatioElec_vs_dR");
  auto pG  = Get<TProfile>(f, "p_RatioPhot_vs_dR");
  auto pT  = Get<TProfile>(f, "p_RatioTot_vs_dR");
  auto pEP = Get<TProfile>(f, "p_ExpOverP_vs_dR");
  auto h2E = Get<TH2D>(f, "h2_RatioElec_vs_dR");
  auto h2G = Get<TH2D>(f, "h2_RatioPhot_vs_dR");
  auto efE  = Get<TEfficiency>(f, "eff_Elec_vs_dR");
  auto efG  = Get<TEfficiency>(f, "eff_Phot_vs_dR");
  auto efEE = Get<TEfficiency>(f, "eff_Elec_vs_E");
  auto efGE = Get<TEfficiency>(f, "eff_Phot_vs_E");
  auto efQE = Get<TEfficiency>(f, "effE_Elec_vs_dR");
  auto efQG = Get<TEfficiency>(f, "effE_Phot_vs_dR");

  // --- MC truth ---
  DrawOverlay(s, hEMCe, "Electron", hEMCg, "Photon", "E_{MC} [GeV]", "Events", "c_mc_energy");
  DrawHist(s, hDRMC, kGray + 2, "#DeltaR_{MC}(e,#gamma)", "Events", "c_mc_deltaR", false);

  // --- Multiplicities ---
  DrawHist(s, hNExp, kColElec, "N_{expected clusters}", "Events", "c_n_expected", false);
  DrawHist(s, hNRem, kColPhot, "N_{remnant clusters}",  "Events", "c_n_remnant",  false);

  // --- Reconstructed energies ---
  DrawHist(s, hEExp,  kColElec, "E_{expected} [GeV]",       "Events", "c_E_expected");
  DrawHist(s, hERem,  kColPhot, "#SigmaE_{remnant} [GeV]",  "Events", "c_E_remnant");
  DrawHist(s, hERemL, kColPhot, "E_{remnant, lead} [GeV]",  "Events", "c_E_remnant_lead");

  // --- Energy response (reco / MC) ---
  DrawOverlay(s, rE, "E_{exp}/E_{MC,e}", rG, "#SigmaE_{rem}/E_{MC,#gamma}",
              "E_{reco} / E_{MC}", "Events", "c_ratio_elec_phot", 1.0);
  DrawHist(s, rE,  kColElec, "E_{exp} / E_{MC,e}",                 "Events", "c_ratio_elec", true, 1.0);
  DrawHist(s, rG,  kColPhot, "#SigmaE_{rem} / E_{MC,#gamma}",      "Events", "c_ratio_phot", true, 1.0);
  DrawHist(s, rGL, kColPhot, "E_{rem,lead} / E_{MC,#gamma}",       "Events", "c_ratio_phot_lead", true, 1.0);
  DrawHist(s, rT,  kColTot,  "(E_{exp}+#SigmaE_{rem}) / (E_{MC,e}+E_{MC,#gamma})", "Events", "c_ratio_total", true, 1.0);

  // --- Comparison with the track ---
  DrawHist(s, eOp, kColElec, "E_{exp} / p_{Tr}",      "Events", "c_Eexp_over_pTr", true, 1.0);
  DrawHist(s, pOm, kColElec, "p_{Tr} / E_{MC,e}",     "Events", "c_pTr_over_EMC", true, 1.0);
  DrawHist(s, rOp, kColPhot, "#SigmaE_{rem} / p_{Tr}", "Events", "c_Erem_over_pTr", true);

  // --- Reco vs MC energy ---
  Draw2D(s, m2e, nullptr, "E_{MC,e} [GeV]",      "E_{exp} [GeV]",            "c_2D_Eexp_vs_EMC", true);
  Draw2D(s, m2g, nullptr, "E_{MC,#gamma} [GeV]", "#SigmaE_{rem} [GeV]",      "c_2D_Erem_vs_EMC", true);

  // --- Remnant <-> Expected in eta-phi (reco level) ---
  Draw2D(s, hDEP, nullptr, "#Delta#eta (rem - exp)", "#Delta#phi [rad]", "c_delta_eta_phi");
  DrawHist(s, hRRE, kGray + 2, "R_{#eta#phi}(remnant, expected)", "Events", "c_R_rem_exp", false);

  // --- Dependence on the e-gamma distance ---
  Draw2D(s, h2E, pE, "#DeltaR_{MC}(e,#gamma)", "E_{exp} / E_{MC,e}",             "c_2D_ratio_elec_vs_dR");
  Draw2D(s, h2G, pG, "#DeltaR_{MC}(e,#gamma)", "#SigmaE_{rem} / E_{MC,#gamma}",  "c_2D_ratio_phot_vs_dR");
  DrawProfiles(s, {pE, pG}, {kColElec, kColPhot},
               {"Electron: E_{exp}/E_{MC,e}", "Photon: #SigmaE_{rem}/E_{MC,#gamma}"},
               "#DeltaR_{MC}(e,#gamma)", "< E_{reco} / E_{MC} >", "c_profile_ratio_vs_dR", 0.0, 1.6, 1.0);
  DrawProfiles(s, {pT}, {kColTot}, {"(E_{exp}+#SigmaE_{rem}) / (E_{MC,e}+E_{MC,#gamma})"},
               "#DeltaR_{MC}(e,#gamma)", "< E_{reco,tot} / E_{MC,tot} >", "c_profile_total_vs_dR", 0.0, 1.6, 1.0);
  DrawProfiles(s, {pEP}, {kColElec}, {"E_{exp} / p_{Tr}"},
               "#DeltaR_{MC}(e,#gamma)", "< E_{exp} / p_{Tr} >", "c_profile_Eexp_over_pTr_vs_dR", 0.0, 1.6, 1.0);

  // --- Efficiencies ---
  DrawEff(s, {efE, efG}, {kColElec, kColPhot}, {"Electron (expected)", "Photon (remnant)"},
          "#DeltaR_{MC}(e,#gamma)", "c_eff_vs_dR", 0.0, 0.5);
  DrawEff(s, {efEE, efGE}, {kColElec, kColPhot}, {"Electron (expected)", "Photon (remnant)"},
          "E_{MC} [GeV]", "c_eff_vs_E", 0.0, 15.0);
  DrawEff(s, {efQE, efQG}, {kColElec, kColPhot},
          {"Electron: |E/E_{MC}-1| < 0.2", "Photon: |E/E_{MC}-1| < 0.2"},
          "#DeltaR_{MC}(e,#gamma)", "c_eff_energy_vs_dR", 0.0, 0.5);

  // ONE multi-page PDF for this region.
  FlushCombined(pdf);
  cout << "  -> " << pdf << endl;
}


// ---------- main ------------------------------------------------------------
void Plotter() {
  SetStyle();
  gROOT->SetBatch(kTRUE);
  for (const auto& s : samples) {
    TFile* f = TFile::Open(("Plots/" + s.label + ".root").c_str(), "READ");
    if (!f || f->IsZombie()) { cerr << "Nie moge otworzyc Plots/" << s.label << ".root\n"; continue; }
    cout << "Rysuje: " << s.label << endl;
    DrawSample(s, f);
    f->Close();
  }
}