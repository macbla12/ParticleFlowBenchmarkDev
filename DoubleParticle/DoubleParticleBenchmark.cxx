// ----------------------------------------------------------------------------
// Macro for analyzing cluster-track matching in e + gamma events
//
// Electron-cluster matching follows the same approach as in
// BenchmarkTrackCluster.cxx:
//   - for EACH cluster, find the nearest track projection (across all
//     segments and points, filtered by system / surface),
//   - if dR(cluster, projection) < MaxREtaPhi, the cluster is "electron-like",
//   - all remaining clusters are treated as "photon-like".
//
// Usage:
//   root -b -q DoubleParticleBenchmark.cxx++
// ----------------------------------------------------------------------------

#define BENCHMARKELECTRONPHOTON_CXX

#include <cmath>
#include <string>
#include <limits>
#include <optional>
#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <iomanip>

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TMath.h>

#include <podio/Frame.h>
#include <podio/ROOTReader.h>
#include <podio/CollectionBase.h>

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
using std::vector;
using std::optional;
using std::numeric_limits;

struct Options {
  string in_file;
  string out_file;
    string clusters;                  // Cluster collection name (e.g. EcalEndcapNClusters)
    string projections;               // Track projection collection name
    uint64_t system;                  // Calorimeter system ID
  bool debug;
    double MaxREtaPhi;                // dR cut for track-cluster matching (electron cluster)
  std::array<double, 2> etaRange;
  std::array<double, 2> ptRange;
};

const Options ElectronPhotonNEndCapOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/ElectronGamma_reco.root",
  "Plots/ElectronPhotonNEndCap.root",
  "EcalEndcapNClusters",
  "CalorimeterTrackProjections",
  103,
  false,
  0.05,
  {-3.6, -1.7},
  {0.0, 10.0}
};

const Options ElectronPhotonPEndCapOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/ElectronGamma_reco.root",
  "Plots/ElectronPhotonPEndCap.root",
  "EcalEndcapPClusters",
  "CalorimeterTrackProjections",
  102, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.04,
  {1.3, 3.6},
  {0.0, 10.0}
};

const Options ElectronPhotonBarrelOptions = {
  "/run/media/epic/Data/Benchmarking/Pairs/ElectronGamma_reco.root",
  "Plots/ElectronPhotonBarrel.root",
  "EcalBarrelClusters",
  "CalorimeterTrackProjections",
  101, //< see https://github.com/eic/epic/blob/main/compact/definitions.xml
  false,
  0.04,
  {-1.8, 1.5},
  {0.0, 20.0}
};

double DeltaPhi(double phi1, double phi2) {
    double dphi = phi2 - phi1;
    if (dphi > TMath::Pi()) dphi -= 2.0 * TMath::Pi();
    else if (dphi <= -TMath::Pi()) dphi += 2.0 * TMath::Pi();
    return dphi;
}

// The function name matches the file name so ROOT can run it with `root -b -q DoubleParticleBenchmark.cxx++`.
void DoubleParticleBenchmark(const Options& opt = ElectronPhotonNEndCapOptions, uint64_t s_use = 1) {

    TFile *Output = new TFile(opt.out_file.c_str(), "RECREATE");

    // --- MC truth histograms ---
    TH1D* h_EMC_Elec      = new TH1D("h_EMC_Elec", "E_{MC} Electron;E_{MC, e} [GeV];Counts", 100, 0, 15);
    TH1D* h_EMC_Phot      = new TH1D("h_EMC_Phot", "E_{MC} Photon;E_{MC, #gamma} [GeV];Counts", 100, 0, 15);
    TH1D* h_DeltaR_MC_e_g = new TH1D("h_DeltaR_MC_e_g", "#DeltaR_{MC}(e, #gamma);#DeltaR;Counts", 100, 0, 0.5);

    // --- Track-cluster matching histograms (as in BenchmarkTrackCluster) ---
    TH1D* h_REtaPhi       = new TH1D("h_REtaPhi", "REtaPhiHist;r;Events", 1000, 0, 0.5);
    TH2D* h_DeltaEtaPhi   = new TH2D("h_DeltaEtaPhi", "DeltaEtaPhiHist;#Delta#eta;#Delta#phi", 1000, -0.3, 0.3, 1000, -0.3, 0.3);

    // --- Cluster energy histograms ---
    TH1D* h_EClus_Elec    = new TH1D("h_EClus_Elec", "Sum E_{Clus} Matched to Track (Electron);E_{Clus, e} [GeV];Events", 100, 0, 15);
    TH1D* h_EClus_Phot    = new TH1D("h_EClus_Phot", "Sum E_{Clus} Unmatched (Photon Candidate);E_{Clus, #gamma} [GeV];Events", 100, 0, 15);

    // --- Comparison histograms (detector response) ---
    TH1D* h_Ratio_Elec    = new TH1D("h_Ratio_Elec", "E_{Clus, e} / E_{MC, e};Ratio;Events", 100, 0, 2);
    TH1D* h_Ratio_Phot    = new TH1D("h_Ratio_Phot", "E_{Clus, #gamma} / E_{MC, #gamma};Ratio;Events", 100, 0, 2);
    TH1D* h_Ratio_Tot     = new TH1D("h_Ratio_Tot", "(E_{Clus, e} + E_{Clus, #gamma}) / (E_{MC, e} + E_{MC, #gamma});Ratio;Events", 100, 0, 2);

    // --- E_cluster / p_track (as in EClusteroverETrack in BenchmarkTrackCluster) ---
    TH1D* h_EClusOverETrack = new TH1D("h_EClusOverETrack", "E_{Cl}/E_{Tr};Ratio E_{Cl}/E_{Tr};Events", 100, 0, 2);

    TH2D* h2_EClus_vs_EMC_Elec = new TH2D("h2_EClus_vs_EMC_Elec", "Electron: E_{Clus} vs E_{MC};E_{MC} [GeV];E_{Clus} [GeV]", 100, 0, 15, 100, 0, 15);
    TH2D* h2_EClus_vs_EMC_Phot = new TH2D("h2_EClus_vs_EMC_Phot", "Photon: E_{Clus} vs E_{MC};E_{MC} [GeV];E_{Clus} [GeV]", 100, 0, 15, 100, 0, 15);

    // Cluster information per event
    TH1D* h_NClustersTot       = new TH1D("h_NClustersTot", "Total Clusters per Event;N_{clusters};Events", 10, -0.5, 9.5);
    TH1D* h_NClustersMatched   = new TH1D("h_NClustersMatched", "Clusters Matched to Track (Electron);N_{clusters};Events", 10, -0.5, 9.5);
    TH1D* h_NClustersUnmatched = new TH1D("h_NClustersUnmatched", "Unmatched Clusters (Photon Candidates);N_{clusters};Events", 10, -0.5, 9.5);

    cout << "\n  Start macro: Electron-Photon Cluster Separation Analysis!" << endl;

    podio::ROOTReader reader = podio::ROOTReader();
    reader.openFile(opt.in_file);

    uint64_t nFrames = reader.getEntries(podio::Category::Event);
    if (opt.debug) nFrames = 20;
    cout << "  Processing " << nFrames << " frames..." << endl;

    uint64_t nClustTotal   = 0;
    uint64_t nClustMatched = 0;

    for (uint64_t iFrame = 0; iFrame < nFrames; ++iFrame) {
        auto frame = podio::Frame(reader.readNextEntry(podio::Category::Event));

        auto& clusters     = frame.get<edm4eic::ClusterCollection>(opt.clusters);
        auto& segments     = frame.get<edm4eic::TrackSegmentCollection>(opt.projections);
        auto& mc_particles = frame.get<edm4hep::MCParticleCollection>("MCParticles");

        h_NClustersTot->Fill(clusters.size());

        // 1. Find the MC truth electron and photon
        TLorentzVector mcElec, mcPhot;
        bool foundElec = false;
        bool foundPhot = false;

        for (edm4hep::MCParticle mc_part : mc_particles) {
            if (mc_part.getGeneratorStatus() != 1 ) continue;

            int pdg = mc_part.getPDG();
            auto mom = mc_part.getMomentum();
            double eta_mc = edm4hep::utils::eta(mom);

            if (eta_mc < opt.etaRange[0] || eta_mc > opt.etaRange[1]) continue;

            if (std::abs(pdg) == 11 && !foundElec) { // Electron
                mcElec.SetPxPyPzE(mom.x, mom.y, mom.z, mc_part.getEnergy());
                foundElec = true;
            } else if (pdg == 22 && !foundPhot) { // Photon
                mcPhot.SetPxPyPzE(mom.x, mom.y, mom.z, mc_part.getEnergy());
                foundPhot = true;
            }
        }

        if (!foundElec || !foundPhot) {
            if (opt.debug) cout << "Frame " << iFrame << ": No MC-level e-gamma pair found. Skipping." << endl;
            continue;
        }

        // Fill the MC histograms
        h_EMC_Elec->Fill(mcElec.E());
        h_EMC_Phot->Fill(mcPhot.E());

        double deltaR_MC = hypot(mcElec.Eta() - mcPhot.Eta(), DeltaPhi(mcElec.Phi(), mcPhot.Phi()));
        h_DeltaR_MC_e_g->Fill(deltaR_MC);

        // 2. Match clusters to tracks, exactly as in BenchmarkTrackCluster:
        //    find the nearest projection for each cluster (across all segments)
        double eClusterElec     = 0.0;   // Total energy of clusters matched to a track
        double eClusterPhotSum  = 0.0;   // Total energy of the remaining clusters
        int matchedClusterCount   = 0;
        int unmatchedClusterCount = 0;

        for (edm4eic::Cluster cluster : clusters) {

            // Cluster eta/phi
            const double etaClust = edm4hep::utils::eta(cluster.getPosition());
            const double phiClust = edm4hep::utils::angleAzimuthal(cluster.getPosition());

            double distMatch = numeric_limits<double>::max();
            double dEta      = numeric_limits<double>::max();
            double dPhi      = numeric_limits<double>::max();

            optional<edm4eic::TrackPoint> match;

            for (edm4eic::TrackSegment segment : segments) {
                for (edm4eic::TrackPoint projection : segment.getPoints()) {
                    const bool isInSystem = (projection.system == opt.system);
                    const bool isAtFace   = (projection.surface == s_use);

                    if (opt.system == 101) {
                        if (!isInSystem) continue;
                    } else {
                        if (!isInSystem || !isAtFace) continue;
                    }

                    const double etaProject = edm4hep::utils::eta(projection.position);
                    const double phiProject = edm4hep::utils::angleAzimuthal(projection.position);

                    const double distProject = hypot(
                        etaClust - etaProject,
                        DeltaPhi(phiClust, phiProject)
                    );

                    if (distProject < distMatch) {
                        distMatch = distProject;
                        dEta      = etaClust - etaProject;
                        dPhi      = DeltaPhi(phiClust, phiProject);
                        match     = projection;
                    }
                }  // end point loop
            }  // end segment loop

            h_DeltaEtaPhi->Fill(dEta, dPhi);
            h_REtaPhi->Fill(distMatch);
            ++nClustTotal;

            if (match.has_value() && distMatch < opt.MaxREtaPhi) {
                // Cluster matched to a track -> electron
                eClusterElec += cluster.getEnergy();
                ++matchedClusterCount;
                ++nClustMatched;

                auto trackMom = match->momentum;
                const double pTrack = sqrt(trackMom.x * trackMom.x +
                                           trackMom.y * trackMom.y +
                                           trackMom.z * trackMom.z);
                // High energy limit pTrack = eTrack
                if (pTrack > 0) h_EClusOverETrack->Fill(cluster.getEnergy() / pTrack);

                if (opt.debug) {
                    cout << "      [Electron cluster] E = " << cluster.getEnergy()
                         << " | dR = " << distMatch << endl;
                }
            } else {
                // Remaining clusters -> photon
                eClusterPhotSum += cluster.getEnergy();
                ++unmatchedClusterCount;

                if (opt.debug) {
                    cout << "      [Photon cluster]   E = " << cluster.getEnergy()
                         << " | dR(min) = " << distMatch << endl;
                }
            }
        }  // end cluster loop

        // 3. Fill statistics and compare energies
        h_NClustersMatched->Fill(matchedClusterCount);
        h_NClustersUnmatched->Fill(unmatchedClusterCount);

        if (eClusterElec > 0) {
            h_EClus_Elec->Fill(eClusterElec);
            h_Ratio_Elec->Fill(eClusterElec / mcElec.E());
            h2_EClus_vs_EMC_Elec->Fill(mcElec.E(), eClusterElec);
        }

        if (eClusterPhotSum > 0) {
            h_EClus_Phot->Fill(eClusterPhotSum);
            h_Ratio_Phot->Fill(eClusterPhotSum / mcPhot.E());
            h2_EClus_vs_EMC_Phot->Fill(mcPhot.E(), eClusterPhotSum);
        }

        double totalMC = mcElec.E() + mcPhot.E();
        double totalClust = eClusterElec + eClusterPhotSum;
        if (totalMC > 0) {
            h_Ratio_Tot->Fill(totalClust / totalMC);
        }

    } // End of event loop

    cout << "  Finished frame loop: " << nClustMatched << "/" << nClustTotal << " clusters matched to tracks." << endl;

    // Write to the output file
    Output->Write();
    Output->Close();

    cout << "  Analysis complete! Histograms saved to: " << opt.out_file << endl;
}