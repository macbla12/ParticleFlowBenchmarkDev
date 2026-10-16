// ----------------------------------------------------------------------------
// Macro for analyzing e + gamma events using the Expected / Remnant clusters
//
//   *ExpectedClusters -> electron (cluster built from the track projection)
//   *RemnantClusters  -> photon   (what remains after the electron is removed)
//
// No own cluster-track matching is done here. The track is used only as a
// reference for the energy (p_track taken from the track projection that
// the Expected cluster sits on).
//
// Usage:
//   root -b -q DoubleParticleBenchmark.cxx++
// ----------------------------------------------------------------------------

#define BENCHMARKELECTRONPHOTON_CXX

#include <cmath>
#include <string>
#include <limits>
#include <iostream>
#include <array>
#include <vector>
#include <algorithm>

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TMath.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>

#include <edm4eic/Cluster.h>
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/TrackPoint.h>
#include <edm4eic/TrackSegment.h>
#include <edm4eic/TrackSegmentCollection.h>

#include <edm4hep/Vector3f.h>
#include <edm4hep/utils/vector_utils.h>
#include <edm4hep/MCParticleCollection.h>

using std::string;
using std::cout;
using std::endl;
using std::numeric_limits;

struct Options {
  string in_file;
  string out_file;
  string expectedClusters;          // electron   (e.g. EcalEndcapNExpectedClusters)
  string remnantClusters;           // photon     (e.g. EcalEndcapNRemnantClusters)
  string projections;               // track projections (for p_track)
  uint64_t system;                  // calorimeter system ID (for projections)
  bool debug;
  std::array<double, 2> etaRange;
  std::array<double, 2> ptRange;
};

const Options ElectronPhotonNEndCapOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/ElectronGamma2_reco.root",
  "Plots/ElectronPhotonNEndCap.root",
  "EcalEndcapNExpectedClusters",
  "EcalEndcapNRemnantClusters",
  "CalorimeterTrackProjections",
  103,
  false,
  {-3.6, -1.7},
  {0.0, 10.0}
};

const Options ElectronPhotonPEndCapOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/ElectronGamma2_reco.root",
  "Plots/ElectronPhotonPEndCap.root",
  "EcalEndcapPExpectedClusters",
  "EcalEndcapPRemnantClusters",
  "CalorimeterTrackProjections",
  102,
  false,
  {1.3, 3.6},
  {0.0, 10.0}
};

const Options ElectronPhotonBarrelOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/ElectronGamma2_reco.root",
  "Plots/ElectronPhotonBarrel.root",
  "EcalBarrelExpectedClusters",
  "EcalBarrelRemnantClusters",
  "CalorimeterTrackProjections",
  101,
  false,
  {-1.8, 1.5},
  {0.0, 20.0}
};

static double DeltaPhi(double phi1, double phi2) {
  double dphi = phi2 - phi1;
  if (dphi > TMath::Pi()) dphi -= 2.0 * TMath::Pi();
  else if (dphi <= -TMath::Pi()) dphi += 2.0 * TMath::Pi();
  return dphi;
}

// p_track from the projection closest (in eta-phi) to a given position.
// Returns -1 if there is no projection for this system / surface.
static double TrackMomentumAt(const edm4eic::TrackSegmentCollection& segments,
                              const edm4hep::Vector3f& pos, uint64_t system, uint64_t surface) {
  const double eta = edm4hep::utils::eta(pos);
  const double phi = edm4hep::utils::angleAzimuthal(pos);
  double best = numeric_limits<double>::max();
  double p    = -1.0;
  for (edm4eic::TrackSegment segment : segments) {
    for (edm4eic::TrackPoint pr : segment.getPoints()) {
      if (pr.system != system) continue;
      if (system != 101 && pr.surface != surface) continue;
      const double d = std::hypot(eta - edm4hep::utils::eta(pr.position),
                                  DeltaPhi(phi, edm4hep::utils::angleAzimuthal(pr.position)));
      if (d < best) {
        best = d;
        p = std::sqrt(pr.momentum.x * pr.momentum.x + pr.momentum.y * pr.momentum.y +
                      pr.momentum.z * pr.momentum.z);
      }
    }
  }
  return p;
}

void DoubleParticleBenchmark(const Options& opt = ElectronPhotonBarrelOptions, uint64_t s_use = 1) {

  gSystem->mkdir(gSystem->DirName(opt.out_file.c_str()), true);
  TFile* Output = new TFile(opt.out_file.c_str(), "RECREATE");

  const int    nDR   = 20;
  const double drMax = 0.5;

  // --- MC truth ---
  TH1D* h_EMC_Elec      = new TH1D("h_EMC_Elec", "E_{MC} electron;E_{MC,e} [GeV];Counts", 100, 0, 15);
  TH1D* h_EMC_Phot      = new TH1D("h_EMC_Phot", "E_{MC} photon;E_{MC,#gamma} [GeV];Counts", 100, 0, 15);
  TH1D* h_DeltaR_MC_e_g = new TH1D("h_DeltaR_MC_e_g", "#DeltaR_{MC}(e,#gamma);#DeltaR;Counts", 100, 0, drMax);

  // --- Multiplicities ---
  TH1D* h_NExpected = new TH1D("h_NExpected", "Expected clusters per event;N;Events", 10, -0.5, 9.5);
  TH1D* h_NRemnant  = new TH1D("h_NRemnant",  "Remnant clusters per event;N;Events", 10, -0.5, 9.5);

  // --- Energies ---
  TH1D* h_EExp   = new TH1D("h_EExp", "Expected cluster (electron);E_{exp} [GeV];Events", 100, 0, 15);
  TH1D* h_ERem   = new TH1D("h_ERem", "Remnant clusters, sum (photon);E_{rem} [GeV];Events", 100, 0, 15);
  TH1D* h_ERemLead = new TH1D("h_ERemLead", "Leading remnant cluster (photon);E_{rem,lead} [GeV];Events", 100, 0, 15);

  // --- Response: reco / MC ---
  TH1D* h_Ratio_Sum_Elec = new TH1D("h_Ratio_Sum_Elec", "E_{exp,sum} / E_{MC,e};Ratio;Events", 100, 0, 2);
  TH1D* h_Ratio_Elec = new TH1D("h_Ratio_Elec", "E_{exp} / E_{MC,e};Ratio;Events", 100, 0, 2);
  TH1D* h_Ratio_Phot = new TH1D("h_Ratio_Phot", "E_{rem} / E_{MC,#gamma};Ratio;Events", 100, 0, 2);
  TH1D* h_Ratio_PhotLead = new TH1D("h_Ratio_PhotLead", "E_{rem,lead} / E_{MC,#gamma};Ratio;Events", 100, 0, 2);
  TH1D* h_Ratio_Tot  = new TH1D("h_Ratio_Tot", "(E_{exp}+E_{rem}) / (E_{MC,e}+E_{MC,#gamma});Ratio;Events", 100, 0, 2);

  // --- Comparison to the track ---
  TH1D* h_ExpOverPTrack  = new TH1D("h_ExpOverPTrack", "E_{exp} / p_{track};Ratio;Events", 100, 0, 2);
  TH1D* h_PTrackOverMC   = new TH1D("h_PTrackOverMC", "p_{track} / E_{MC,e};Ratio;Events", 100, 0, 2);
  TH1D* h_RemOverPTrack  = new TH1D("h_RemOverPTrack", "E_{rem} / p_{track};Ratio;Events", 100, 0, 2);

  // --- 2D: reco vs MC ---
  TH2D* h2_EExp_vs_EMC = new TH2D("h2_EExp_vs_EMC", "Electron: E_{exp} vs E_{MC};E_{MC} [GeV];E_{exp} [GeV]", 100, 0, 15, 100, 0, 15);
  TH2D* h2_ERem_vs_EMC = new TH2D("h2_ERem_vs_EMC", "Photon: E_{rem} vs E_{MC};E_{MC} [GeV];E_{rem} [GeV]", 100, 0, 15, 100, 0, 15);

  // --- Eta-phi distance remnant <-> expected (reco level) ---
  TH1D* h_RRemExp     = new TH1D("h_RRemExp", "R_{#eta#phi}(remnant, expected);R;Events", 200, 0, drMax);
  TH2D* h_DeltaEtaPhi = new TH2D("h_DeltaEtaPhi", "Remnant - expected;#Delta#eta;#Delta#phi", 200, -0.3, 0.3, 200, -0.3, 0.3);

  // --- Dependence on the e-gamma distance (MC level) ---
  TProfile* p_RatioElec_vs_dR = new TProfile("p_RatioElec_vs_dR", "<E_{exp}/E_{MC,e}> vs #DeltaR_{MC}(e,#gamma);#DeltaR_{MC};<E_{exp}/E_{MC,e}>", nDR, 0, drMax);
  TProfile* p_RatioPhot_vs_dR = new TProfile("p_RatioPhot_vs_dR", "<E_{rem}/E_{MC,#gamma}> vs #DeltaR_{MC}(e,#gamma);#DeltaR_{MC};<E_{rem}/E_{MC,#gamma}>", nDR, 0, drMax);
  TProfile* p_RatioTot_vs_dR  = new TProfile("p_RatioTot_vs_dR", "<E_{tot}^{reco}/E_{tot}^{MC}> vs #DeltaR_{MC}(e,#gamma);#DeltaR_{MC};ratio", nDR, 0, drMax);
  TProfile* p_ExpOverP_vs_dR  = new TProfile("p_ExpOverP_vs_dR", "<E_{exp}/p_{track}> vs #DeltaR_{MC}(e,#gamma);#DeltaR_{MC};<E_{exp}/p_{track}>", nDR, 0, drMax);
  TH2D* h2_RatioPhot_vs_dR    = new TH2D("h2_RatioPhot_vs_dR", "E_{rem}/E_{MC,#gamma} vs #DeltaR_{MC};#DeltaR_{MC};E_{rem}/E_{MC,#gamma}", nDR, 0, drMax, 100, 0, 2);
  TH2D* h2_RatioElec_vs_dR    = new TH2D("h2_RatioElec_vs_dR", "E_{exp}/E_{MC,e} vs #DeltaR_{MC};#DeltaR_{MC};E_{exp}/E_{MC,e}", nDR, 0, drMax, 100, 0, 2);

  // --- Reconstruction efficiency vs distance / energy ---
  // "found" = at least one Expected (electron) / Remnant (photon) cluster
  TEfficiency* eff_Elec_vs_dR = new TEfficiency("eff_Elec_vs_dR", "Electron (expected cluster) found;#DeltaR_{MC}(e,#gamma);Efficiency", nDR, 0, drMax);
  TEfficiency* eff_Phot_vs_dR = new TEfficiency("eff_Phot_vs_dR", "Photon (remnant cluster) found;#DeltaR_{MC}(e,#gamma);Efficiency", nDR, 0, drMax);
  TEfficiency* eff_Elec_vs_E  = new TEfficiency("eff_Elec_vs_E", "Electron (expected cluster) found;E_{MC,e} [GeV];Efficiency", 30, 0, 15);
  TEfficiency* eff_Phot_vs_E  = new TEfficiency("eff_Phot_vs_E", "Photon (remnant cluster) found;E_{MC,#gamma} [GeV];Efficiency", 30, 0, 15);
  // "good" = found and energy within 20% of MC
  TEfficiency* effE_Phot_vs_dR = new TEfficiency("effE_Phot_vs_dR", "Photon: found and |E_{rem}/E_{MC}-1|<0.2;#DeltaR_{MC}(e,#gamma);Efficiency", nDR, 0, drMax);
  TEfficiency* effE_Elec_vs_dR = new TEfficiency("effE_Elec_vs_dR", "Electron: found and |E_{exp}/E_{MC}-1|<0.2;#DeltaR_{MC}(e,#gamma);Efficiency", nDR, 0, drMax);

  cout << "\n  Start macro: Expected/Remnant cluster analysis (e + gamma)" << endl;

  podio::ROOTReader reader = podio::ROOTReader();
  reader.openFile(opt.in_file);

  uint64_t nFrames = reader.getEntries(podio::Category::Event);
  if (opt.debug) nFrames = std::min<uint64_t>(nFrames, 20);
  cout << "  Processing " << nFrames << " frames..." << endl;

  uint64_t nUsed = 0, nElecFound = 0, nPhotFound = 0;

  for (uint64_t iFrame = 0; iFrame < nFrames; ++iFrame) {
    auto frame = podio::Frame(reader.readNextEntry(podio::Category::Event));

    auto& expected     = frame.get<edm4eic::ClusterCollection>(opt.expectedClusters);
    auto& remnant      = frame.get<edm4eic::ClusterCollection>(opt.remnantClusters);
    auto& segments     = frame.get<edm4eic::TrackSegmentCollection>(opt.projections);
    auto& mc_particles = frame.get<edm4hep::MCParticleCollection>("MCParticles");

    h_NExpected->Fill(expected.size());
    h_NRemnant->Fill(remnant.size());

    // 1. MC truth e and gamma
    TLorentzVector mcElec, mcPhot;
    bool foundElec = false, foundPhot = false;

    for (edm4hep::MCParticle mc_part : mc_particles) {
      if (mc_part.getGeneratorStatus() != 1) continue;
      const int pdg = mc_part.getPDG();
      auto mom = mc_part.getMomentum();
      const double eta_mc = edm4hep::utils::eta(mom);
      if (eta_mc < opt.etaRange[0] || eta_mc > opt.etaRange[1]) continue;

      if (std::abs(pdg) == 11 && !foundElec) {
        mcElec.SetPxPyPzE(mom.x, mom.y, mom.z, mc_part.getEnergy());
        foundElec = true;
      } else if (pdg == 22 && !foundPhot) {
        mcPhot.SetPxPyPzE(mom.x, mom.y, mom.z, mc_part.getEnergy());
        foundPhot = true;
      }
    }
    if (!foundElec || !foundPhot) {
      if (opt.debug) cout << "Frame " << iFrame << ": no MC e-gamma pair. Skipping." << endl;
      continue;
    }
    ++nUsed;

    const double eMCe = mcElec.E();
    const double eMCg = mcPhot.E();
    const double dRMC = std::hypot(mcElec.Eta() - mcPhot.Eta(), DeltaPhi(mcElec.Phi(), mcPhot.Phi()));

    h_EMC_Elec->Fill(eMCe);
    h_EMC_Phot->Fill(eMCg);
    h_DeltaR_MC_e_g->Fill(dRMC);

    // 2. Electron = Expected cluster (the most energetic if there are several)
    bool   hasExp = false;
    double eExp = 0.0, eExpSum = 0.0;
    edm4hep::Vector3f posExp{0, 0, 0};
    for (edm4eic::Cluster c : expected) {
      eExpSum += c.getEnergy();
      if (!hasExp || c.getEnergy() > eExp) {
        eExp = c.getEnergy();
        posExp = c.getPosition();
        hasExp = true;
      }
    }

    // 3. Photon = Remnant clusters (sum + leading)
    bool   hasRem = false;
    double eRemSum = 0.0, eRemLead = 0.0;
    edm4hep::Vector3f posRemLead{0, 0, 0};
    for (edm4eic::Cluster c : remnant) {
      eRemSum += c.getEnergy();
      if (!hasRem || c.getEnergy() > eRemLead) {
        eRemLead = c.getEnergy();
        posRemLead = c.getPosition();
        hasRem = true;
      }
    }
    hasExp = hasExp && eExp > 0;
    hasRem = hasRem && eRemSum > 0;

    // 4. Efficiencies (denominator = all MC e-gamma events)
    eff_Elec_vs_dR->Fill(hasExp, dRMC);
    eff_Phot_vs_dR->Fill(hasRem, dRMC);
    eff_Elec_vs_E->Fill(hasExp, eMCe);
    eff_Phot_vs_E->Fill(hasRem, eMCg);
    if (hasExp) ++nElecFound;
    if (hasRem) ++nPhotFound;

    effE_Elec_vs_dR->Fill(hasExp && std::abs(eExp / eMCe - 1.0) < 0.2, dRMC);
    effE_Phot_vs_dR->Fill(hasRem && std::abs(eRemSum / eMCg - 1.0) < 0.2, dRMC);

    // 5. Electron energy response + comparison to the track
    double pTrack = -1.0;
    if (hasExp) {
      h_EExp->Fill(eExpSum);
      h_Ratio_Sum_Elec->Fill(eExpSum/ eMCe);
      h_Ratio_Elec->Fill(eExp / eMCe);
      h2_EExp_vs_EMC->Fill(eMCe, eExp);
      p_RatioElec_vs_dR->Fill(dRMC, eExp / eMCe);
      h2_RatioElec_vs_dR->Fill(dRMC, eExp / eMCe);

      pTrack = TrackMomentumAt(segments, posExp, opt.system, s_use);
      if (pTrack > 0) {
        h_ExpOverPTrack->Fill(eExp / pTrack);
        h_PTrackOverMC->Fill(pTrack / eMCe);
        p_ExpOverP_vs_dR->Fill(dRMC, eExp / pTrack);
      }
    }

    // 6. Photon energy response
    if (hasRem) {
      h_ERem->Fill(eRemSum);
      h_ERemLead->Fill(eRemLead);
      h_Ratio_Phot->Fill(eRemSum / eMCg);
      h_Ratio_PhotLead->Fill(eRemLead / eMCg);
      h2_ERem_vs_EMC->Fill(eMCg, eRemSum);
      p_RatioPhot_vs_dR->Fill(dRMC, eRemSum / eMCg);
      h2_RatioPhot_vs_dR->Fill(dRMC, eRemSum / eMCg);
      if (pTrack > 0) h_RemOverPTrack->Fill(eRemSum / pTrack);
    }

    // 7. Remnant <-> Expected distance in eta-phi (reco level)
    if (hasExp && hasRem) {
      const double dEta = edm4hep::utils::eta(posRemLead) - edm4hep::utils::eta(posExp);
      const double dPhi = DeltaPhi(edm4hep::utils::angleAzimuthal(posExp), edm4hep::utils::angleAzimuthal(posRemLead));
      h_DeltaEtaPhi->Fill(dEta, dPhi);
      h_RRemExp->Fill(std::hypot(dEta, dPhi));

      const double ratioTot = (eExp + eRemSum) / (eMCe + eMCg);
      h_Ratio_Tot->Fill(ratioTot);
      p_RatioTot_vs_dR->Fill(dRMC, ratioTot);
    }

    if (opt.debug) {
      cout << "Frame " << iFrame << ": dR_MC=" << dRMC
           << " | E_MC(e,g)=" << eMCe << "," << eMCg
           << " | E_exp=" << eExp << " | E_rem(sum)=" << eRemSum
           << " | p_track=" << pTrack << endl;
    }
  }  // end event loop

  cout << "  Events with MC e-gamma pair: " << nUsed << endl;
  cout << "  Expected (electron) found:   " << nElecFound << "/" << nUsed << endl;
  cout << "  Remnant (photon) found:      " << nPhotFound << "/" << nUsed << endl;

  Output->Write();
  Output->Close();

  cout << "  Analysis complete! Histograms saved to: " << opt.out_file << endl;
}