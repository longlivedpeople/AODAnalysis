//=======================================================================================================================================================================================================================//ryAOD" 
//                                                                                                                                                                                                                       // 
//$$$$$$\ $$$$$$$$\  $$$$$$\   $$$$$$\                      $$\                                    $$\       $$\                            $$\  $$$$$$\                      $$\                     $$                 //
//\_$$  _|$$  _____|$$  __$$\ $$  __$$\                     $$ |                                   $$ |      \__|                           $$ |$$  __$$\                     $$ |                    \__|               //
//  $$ |  $$ |      $$ /  \__|$$ /  $$ |                    $$ |      $$$$$$\  $$$$$$$\   $$$$$$\  $$ |      $$\ $$\    $$\  $$$$$$\   $$$$$$$ |$$ /  $$ |$$$$$$$\   $$$$$$\  $$ |$$\   $$\  $$$$$$$\ $$\  $$$$$$$       //
//  $$ |  $$$$$\    $$ |      $$$$$$$$ |      $$$$$$\       $$ |     $$  __$$\ $$  __$$\ $$  __$$\ $$ |      $$ |\$$\  $$  |$$  __$$\ $$  __$$ |$$$$$$$$ |$$  __$$\  \____$$\ $$ |$$ |  $$ |$$  _____|$$ |$$  _____|     //
//  $$ |  $$  __|   $$ |      $$  __$$ |      \______|      $$ |     $$ /  $$ |$$ |  $$ |$$ /  $$ |$$ |      $$ | \$$\$$  / $$$$$$$$ |$$ /  $$ |$$  __$$ |$$ |  $$ | $$$$$$$ |$$ |$$ |  $$ |\$$$$$$\  $$ |\$$$$$$        // 
//  $$ |  $$ |      $$ |  $$\ $$ |  $$ |                    $$ |     $$ |  $$ |$$ |  $$ |$$ |  $$ |$$ |      $$ |  \$$$  /  $$   ____|$$ |  $$ |$$ |  $$ |$$ |  $$ |$$  __$$ |$$ |$$ |  $$ | \____$$\ $$ | \____$$       //
//$$$$$$\ $$ |      \$$$$$$  |$$ |  $$ |                    $$$$$$$$\\$$$$$$  |$$ |  $$ |\$$$$$$$ |$$$$$$$$\ $$ |   \$  /   \$$$$$$$\ \$$$$$$$ |$$ |  $$ |$$ |  $$ |\$$$$$$$ |$$ |\$$$$$$$ |$$$$$$$  |$$ |$$$$$$$  |     //
//\______|\__|       \______/ \__|  \__|                    \________|\______/ \__|  \__| \____$$ |\________|\__|    \_/     \_______| \_______|\__|  \__|\__|  \__| \_______|\__| \____$$ |\_______/ \__|\_______/      //
//                                                                                       $$\   $$ |                                                                               $$\   $$ |                             // 
//                                                                                       \$$$$$$  |                                                                               \$$$$$$  |                             //
//=======================================================================================================================================================================================================================//
//                                                                                                                                                                                                                       //
// Authors of the code: Celia Fernandez Madrazo                                                                                                                                                                          //
//                      Pablo Martinez Ruiz Del Arbol                                                                                                                                                                    //
//                                                                                                                                                                                                                       //
//=======================================================================================================================================================================================================================//
//                                                                                                                                                                                                                       //
// Description: Main analyzer                                                                                                                                                                                            //
//                                                                                                                                                                                                                       //
//=======================================================================================================================================================================================================================//



#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"


#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"



#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/TrimmedVertexFit/interface/TrimmedVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"


#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Candidate/interface/Candidate.h"


#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"



//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"



//#include "DataFormats/MuonReco/interface/MuonFwd.h" 

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"


//=======================================================================================================================================================================================================================//


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// FUNCTIONS ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

bool generalTrackPreselection(const reco::Track & track)
{

    // Return trie if the general track fulfils very soft requirements
    
    if (fabs(track.eta()) > 2.0) { return false; }
    if (track.pt() < 20) { return false; }

    return true;

}


bool photonPreselection(const reco::Photon & photon)
{

    // Return true if the photon fulfills with the analysis requirements and false instead

    if (fabs(photon.eta()) > 1.4442 && fabs(photon.eta()) < 1.566) { return false; } // narrow EB region to be defined
    if (photon.hadronicOverEm() > 0.05) { return false; }
    if (photon.isEE() && photon.full5x5_sigmaIetaIeta() > 0.034) { return false; }
    if (photon.isEB() && photon.full5x5_sigmaIetaIeta() > 0.012) { return false; }

    return true;

}


bool isLongLivedLepton(const reco::GenParticle &p)
{

    // Return true if the genparticle is one of the four displaced leptons of the event

    if (!( abs(p.pdgId()) == 11 || abs(p.pdgId()) == 13)){ return false; }
    if (abs(p.mother()->pdgId()) != 1000022){ return false; }

    return true;

}


/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// DATA DEFINITION //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////// BRANCHES /////////////////////////////////////

//-> EVENT INFO
Int_t Event_event;
Int_t Event_luminosityBlock;
Int_t Event_run;


//-> PRIMARY VERTEX SELECTION
Int_t nPV;
Float_t PV_vx;
Float_t PV_vy;
Float_t PV_vz;
Int_t PV_nTracks;
Int_t PV_fittingtracks;
Int_t PV_tracksSize;

// -> REFITTED PRIMARY VERTEX
Float_t RefittedPV_vx;
Float_t RefittedPV_vy;
Float_t RefittedPV_vz;


//-> BEAM SPOT
Float_t BeamSpot_x0;
Float_t BeamSpot_y0;
Float_t BeamSpot_z0;
Float_t BeamSpot_BeamWidthX;
Float_t BeamSpot_BeamWidthY;


//-> GENERAL TRACKS + DISPLACED TRACKS
const Int_t nTrackMax = 100000;
Int_t nTrack;
Float_t TrackSel_pt[nTrackMax];
Float_t TrackSel_eta[nTrackMax];
Float_t TrackSel_phi[nTrackMax];
Float_t TrackSel_dxy[nTrackMax];
Float_t TrackSel_dxyError[nTrackMax];
Float_t TrackSel_d0[nTrackMax];
Float_t TrackSel_d0Error[nTrackMax];
Float_t TrackSel_dz[nTrackMax];
Float_t TrackSel_dzError[nTrackMax];
Float_t TrackSel_vx[nTrackMax];
Float_t TrackSel_vy[nTrackMax];
Float_t TrackSel_vz[nTrackMax];
Int_t TrackSel_isDisplaced[nTrackMax];
Int_t TrackSel_trackIndex[nTrackMax];
Int_t TrackSel_numberOfValidTrackerHits[nTrackMax];
Int_t TrackSel_numberOfValidPixelHits[nTrackMax];
Int_t TrackSel_numberOfValidPixelBarrelHits[nTrackMax];
Int_t TrackSel_numberOfValidPixelEndcapHits[nTrackMax];
Int_t TrackSel_numberOfValidStripHits[nTrackMax];
Int_t TrackSel_numberOfValidStripTIBHits[nTrackMax];
Int_t TrackSel_numberOfValidStripTIDHits[nTrackMax];
Int_t TrackSel_numberOfValidStripTOBHits[nTrackMax];
Int_t TrackSel_numberOfValidStripTECHits[nTrackMax];


//-> PHOTONS
const Int_t nPhotonMax = 100;
Int_t nPhoton;
Int_t nPhotonOriginal;
Float_t PhotonSel_et[nPhotonMax];
Float_t PhotonSel_eta[nPhotonMax];
Float_t PhotonSel_phi[nPhotonMax];
Float_t PhotonSel_hadronicOverEm[nPhotonMax];
Float_t PhotonSel_full5x5_sigmaIetaIeta[nPhotonMax];
Int_t PhotonSel_isEB[nPhotonMax];
Int_t PhotonSel_isEE[nPhotonMax];

//-> MUON TRIGGER OBJECTS
const Int_t nMuonTriggerObjectMax = 500;
Int_t nMuonTriggerObject;
Float_t MuonTriggerObjectSel_pt[nMuonTriggerObjectMax];
Float_t MuonTriggerObjectSel_eta[nMuonTriggerObjectMax];
Float_t MuonTriggerObjectSel_phi[nMuonTriggerObjectMax];


//-> GENPARTICLES (DISPLACED LEPTONS)
const Int_t nGenParticleMax = 500;
Int_t nGenParticle;
Float_t GenParticleSel_pt[nGenParticleMax];
Float_t GenParticleSel_eta[nGenParticleMax];
Float_t GenParticleSel_phi[nGenParticleMax];
Int_t GenParticleSel_pdgId[nGenParticleMax];


//-> TRIGGER INFO
const std::string muonPathName = "HLT_DoubleMu43NoFiltersNoVtx_v3"; // default
Bool_t HLT_DoubleMu43NoFiltersNoVtx_v3;

// -> ELECTRON CANDIDATES
const Int_t nElectronCandidateMax = 100;
Int_t nElectronCandidate;
Float_t ElectronCandidate_pt[nElectronCandidateMax];
Float_t ElectronCandidate_et[nElectronCandidateMax];
Float_t ElectronCandidate_eta[nElectronCandidateMax];
Float_t ElectronCandidate_phi[nElectronCandidateMax];
Int_t ElectronCandidate_photonIdx[nElectronCandidateMax];
Int_t ElectronCandidate_trackIdx[nElectronCandidateMax];


//-> MUON CANDIDATES
const Int_t nMuonCandidateMax = 100;
Int_t nMuonCandidate;
Float_t MuonCandidate_pt[nMuonCandidateMax];
Float_t MuonCandidate_eta[nMuonCandidateMax];
Float_t MuonCandidate_phi[nMuonCandidateMax];
Int_t MuonCandidate_muonTriggerObjectIdx[nMuonCandidateMax];
Int_t MuonCandidate_trackIdx[nMuonCandidateMax];


/////////////////////////////////////// OUTPUT //////////////////////////////////////

TFile *file_out;
TTree *tree_out;


//=======================================================================================================================================================================================================================//
class AODAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit AODAnalysis(const edm::ParameterSet&);
      ~AODAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      std::string output_filename;
      edm::ParameterSet parameters;


      edm::EDGetTokenT<edm::View<reco::Vertex> > thePrimaryVertexCollection;
      edm::EDGetTokenT<edm::View<reco::Track> > theGeneralTrackCollection;
      edm::EDGetTokenT<edm::View<reco::Track> > theDisplacedTrackCollection;
      edm::EDGetTokenT<reco::BeamSpot> theBeamSpot;
      edm::EDGetTokenT<edm::View<reco::Photon> > thePhotonCollection;
      edm::EDGetTokenT<edm::View<reco::GenParticle> >  theGenParticleCollection;


      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<trigger::TriggerEvent> triggerSummary_;


      //HLTConfigProvider hltConfig_;


};
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
AODAnalysis::AODAnalysis(const edm::ParameterSet& iConfig)
{
   usesResource("TFileService");
   
   parameters = iConfig;


   thePrimaryVertexCollection = consumes<edm::View<reco::Vertex> >  (parameters.getParameter<edm::InputTag>("PrimaryVertexCollection"));
   theGeneralTrackCollection = consumes<edm::View<reco::Track> > (parameters.getParameter<edm::InputTag>("GeneralTrackCollection"));
   theDisplacedTrackCollection = consumes<edm::View<reco::Track> > (parameters.getParameter<edm::InputTag>("DisplacedTrackCollection"));
   theBeamSpot = consumes<reco::BeamSpot>  (parameters.getParameter<edm::InputTag>("BeamSpot"));
   thePhotonCollection = consumes<edm::View<reco::Photon> > (parameters.getParameter<edm::InputTag>("PhotonCollection"));
   theGenParticleCollection = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("GenParticleCollection"));


   triggerBits_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("bits"));
   triggerSummary_ = consumes<trigger::TriggerEvent> (parameters.getParameter<edm::InputTag>("triggerSummary"));

}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
AODAnalysis::~AODAnalysis()
{

}
//=======================================================================================================================================================================================================================//



//=======================================================================================================================================================================================================================//
void AODAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////// MAIN CODE /////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////



   //////////////////////////////// GET THE COLLECTIONS ////////////////////////////////
   
   //-> Handle declaration
   edm::Handle<edm::View<reco::Vertex> > primaryVertices;
   edm::Handle<edm::View<reco::Track> > generalTracks;
   edm::Handle<edm::View<reco::Track> > displacedTracks;
   edm::Handle<reco::BeamSpot> beamSpot;
   edm::Handle<edm::View<reco::Photon> > photons;
   edm::Handle<edm::View<reco::GenParticle> > genParticles;

   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<trigger::TriggerEvent> trigSummary;

   //-> Get by Token from the file
   iEvent.getByToken(thePrimaryVertexCollection, primaryVertices);
   iEvent.getByToken(theGeneralTrackCollection, generalTracks);
   iEvent.getByToken(theDisplacedTrackCollection, displacedTracks);
   iEvent.getByToken(theBeamSpot, beamSpot);
   iEvent.getByToken(thePhotonCollection, photons);
   iEvent.getByToken(theGenParticleCollection, genParticles);

   iEvent.getByToken(triggerBits_, triggerBits);
   iEvent.getByToken(triggerSummary_, trigSummary);


   /////////////////////////////////// EVENT INFO //////////////////////////////////////


   Event_event = iEvent.id().event();
   Event_run = iEvent.id().run();
   Event_luminosityBlock = iEvent.id().luminosityBlock();


   //////////////////////////////////// BEAM SPOT //////////////////////////////////////

   reco::BeamSpot beamSpotObject = *beamSpot;
   BeamSpot_x0 = beamSpotObject.x0();
   BeamSpot_y0 = beamSpotObject.y0();
   BeamSpot_z0 = beamSpotObject.z0();
   BeamSpot_BeamWidthX = beamSpotObject.BeamWidthX();
   BeamSpot_BeamWidthY = beamSpotObject.BeamWidthY();



   ///////////////////////////////// TRACK COLLECTION //////////////////////////////////

   std::vector<reco::Track> trackCollection;

   std::vector<float> displacedPt;
   std::vector<float> displacedEta;
   std::vector<float> displacedPhi;

   std::vector<float> allTrackPt; // std vector to sort the tracks by pt
   std::vector<int> isDisplaced; // 0 for generalTrack collection and 1 for displacedTrack collection
   std::vector<int> trackIndex; // index of the track to obtain the values

   int nTrackCollection_general = 0;
   int nTrackCollection_displaced = 0;
   int nGeneralTrack = generalTracks->size();
   int nDisplacedTrack = displacedTracks->size();


   // Include first the displaced Collection because is a smaller collection (easier to look for overlaps)

   // Loop over displaced tracks collection (include all):

   for (size_t i = 0; i < displacedTracks->size(); i++){

       const reco::Track track = (*displacedTracks)[i];

       trackCollection.push_back(track);

       displacedPt.push_back(track.pt()); displacedEta.push_back(track.eta()); displacedPhi.push_back(track.phi());
       allTrackPt.push_back(track.pt()); isDisplaced.push_back(1); trackIndex.push_back(i);

       nTrackCollection_displaced++;       
   }

   // Loop over general tracks collection (skip those that are also in the displaced collection):
   
   
   for (size_t i = 0; i < generalTracks->size(); i++){

       const reco::Track track = (*generalTracks)[i];       
       
       if(std::find(displacedPt.begin(), displacedPt.end(), track.pt()) != displacedPt.end() && std::find(displacedEta.begin(), displacedEta.end(), track.eta()) != displacedEta.end() && std::find(displacedPhi.begin(), displacedPhi.end(), track.phi()) != displacedPhi.end()) {
           continue;
       }
       

       trackCollection.push_back(track);
       allTrackPt.push_back(track.pt()); isDisplaced.push_back(0); trackIndex.push_back(i);
       nTrackCollection_general++;
   }


   std::vector<int> iTrack; 
   for (size_t i = 0; i < trackCollection.size(); i++) {iTrack.push_back(i);}

   // sort all tracks by pt
   std::sort( std::begin(iTrack), std::end(iTrack), [&](int i1, int i2){ return allTrackPt.at(i1) < allTrackPt.at(i2); });
  // std::sort( std::begin(trackIndex), std::end(trackIndex), [&](int i1, int i2){ return allTrackPt.at(i1) < allTrackPt.at(i2); });
   //std::sort( std::begin(isDisplaced), std::end(isDisplaced), [&](int i1, int i2){ return allTrackPt.at(i1) < allTrackPt.at(i2); });



   ////////////////////////////////// GENERAL TRACKS //////////////////////////////////

   std::vector<int> iT; // track selection indexes

   for (size_t i = 0; i < trackCollection.size(); i++){

       reco::Track track = trackCollection.at(i);

       if (generalTrackPreselection(track)) { iT.push_back(i); }

   }


   // Loop over the preselected tracks
   for (size_t i = 0; i < iT.size(); i++){

       reco::Track track = trackCollection.at(iT.at(i)); // get the track
       
       // Get the variables
       TrackSel_pt[i] = track.pt();
       TrackSel_eta[i] = track.eta();
       TrackSel_phi[i] = track.phi();
       TrackSel_dxy[i] = track.dxy();
       TrackSel_d0[i] = track.d0();
       TrackSel_dxyError[i] = track.dxyError();
       TrackSel_d0Error[i] = track.d0Error();
       TrackSel_dz[i] = track.dz();
       TrackSel_dzError[i] = track.dzError();
       TrackSel_vx[i] = track.vx();
       TrackSel_vy[i] = track.vy();
       TrackSel_vz[i] = track.vz(); 

       // Track nature information
       TrackSel_isDisplaced[i] = isDisplaced.at(iT.at(i));
       TrackSel_trackIndex[i] = trackIndex.at(iT.at(i));


       // Hit info
       const reco::HitPattern &hits = track.hitPattern();

       TrackSel_numberOfValidTrackerHits[i] = hits.numberOfValidTrackerHits();
       TrackSel_numberOfValidPixelHits[i] = hits.numberOfValidPixelHits();
       TrackSel_numberOfValidPixelBarrelHits[i] = hits.numberOfValidPixelBarrelHits();
       TrackSel_numberOfValidPixelEndcapHits[i] = hits.numberOfValidPixelEndcapHits();
       TrackSel_numberOfValidStripHits[i] = hits.numberOfValidStripHits();
       TrackSel_numberOfValidStripTIBHits[i] = hits.numberOfValidStripTIBHits();
       TrackSel_numberOfValidStripTIDHits[i] = hits.numberOfValidStripTIDHits();
       TrackSel_numberOfValidStripTOBHits[i] = hits.numberOfValidStripTOBHits();
       TrackSel_numberOfValidStripTECHits[i] = hits.numberOfValidStripTECHits();

   }


   nTrack = iT.size(); // number of selected tracks




   ////////////////////////////// PRIMARY VERTEX FEATURES //////////////////////////////


   nPV = primaryVertices->size();
   const reco::Vertex *vertex = &((*primaryVertices)[0]);

   PV_vx = vertex->position().x();
   PV_vy = vertex->position().y();
   PV_vz = vertex->position().z();


   PV_nTracks = vertex->nTracks(); // number of tracks with weights > 0.5
   PV_tracksSize= vertex->tracksSize(); // number of tracks associated to the vertex


   ////////////////////////////////// PHOTON FEATURES //////////////////////////////////
   
   std::vector<int> iP; // photon indexes

   // Loop over the photons
   for (size_t i = 0; i < photons->size(); i++){

       const reco::Photon & photon = (*photons)[i];

       if (photonPreselection(photon)) { iP.push_back(i); }

   }


   // Sort preselected photon indexes by et
   std::sort( std::begin(iP), std::end(iP), [&](int i1, int i2){ return photons->at(i1).et() < photons->at(i2).et(); });


   // Loop over the preselected photons
   for (size_t i = 0; i < iP.size(); ++i){

       const reco::Photon & photon = (*photons)[iP.at(i)];

       PhotonSel_et[i] = photon.et();
       PhotonSel_eta[i] = photon.eta();
       PhotonSel_phi[i] = photon.phi();
       PhotonSel_hadronicOverEm[i] = photon.hadronicOverEm();
       PhotonSel_full5x5_sigmaIetaIeta[i] = photon.full5x5_sigmaIetaIeta();
       PhotonSel_isEB[i] = photon.isEB();
       PhotonSel_isEE[i] = photon.isEE();


   }

   nPhoton = iP.size(); // number of preselected photons
   nPhotonOriginal = photons->size(); // original number of photons in AOD


   ////////////////////////////// MUON TRIGGER OBJECTS /////////////////////////////////

   // -> TRIGGER INFO

   // Access the trigger names of the event
   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerBits); 

   std::string muonPathName = "HLT_DoubleMu43NoFiltersNoVtx_v3"; // default

   HLT_DoubleMu43NoFiltersNoVtx_v3 = triggerBits->accept(trigNames.triggerIndex(muonPathName));  



   // -> MUON OBJECTS
   std::vector<int> iMT; // muon trigger object indexes

   
   // Access the muon objects
   trigger::size_type filterIndex = trigSummary->filterIndex(edm::InputTag("hltL3fL1sMu5EG20orMu20EG15L1f5L2NVf16L3NoFiltersNoVtxFiltered48","", "HLT")); // the filters has to be defined 
   trigger::TriggerObjectCollection triggerObjects = trigSummary->getObjects();

   // filterIndex only less if the filter is present
   if( !(filterIndex >= trigSummary->sizeFilters()) ) {

       const trigger::Keys& keys = trigSummary->filterKeys( filterIndex );

       // Loop over the keys of the objects associated to the filter
       for(unsigned short key : keys) {
       
           iMT.push_back(key);
       
       }

   }


   // Sort muon object triggers by pt:
   std::sort(std::begin(iMT), std::end(iMT), [&](int i1, int i2){ return triggerObjects[i1].pt() < triggerObjects[i2].pt(); });


   // Store the muon trigger objects
   for (size_t i = 0; i < iMT.size(); i++){

       trigger::TriggerObject foundObject = triggerObjects[iMT.at(i)];

       MuonTriggerObjectSel_pt[i] = foundObject.pt();
       MuonTriggerObjectSel_eta[i] = foundObject.eta();
       MuonTriggerObjectSel_phi[i] = foundObject.phi();

   }

   nMuonTriggerObject = iMT.size(); // number of found muon trigger objects


   //////////////////////////////// GENPARTICLE FEATURES ///////////////////////////////
   
   std::vector<int> iGP;

   
   for(size_t i = 0; i < genParticles->size(); i++){

        const reco::GenParticle &genparticle = (*genParticles)[i];

        if (isLongLivedLepton(genparticle)){ iGP.push_back(i); }


   }

   // Loop over the selected gen particles
   for(size_t i = 0; i < iGP.size(); i++){

       const reco::GenParticle &genparticle = (*genParticles)[iGP.at(i)];

       GenParticleSel_pdgId[i] = genparticle.pdgId();

       // Get the last genparticle (to avoid radiative effects):
       if (genparticle.numberOfDaughters() > 0){

           const reco::Candidate *d = genparticle.daughter(0);
           while(d->numberOfDaughters()> 0 && d->daughter(0)->pdgId() == d->pdgId()){ d = d->daughter(0); }

           GenParticleSel_pt[i] = d->pt();
           GenParticleSel_eta[i] = d->eta();
           GenParticleSel_phi[i] = d->phi();

       } else {

           GenParticleSel_pt[i] = genparticle.pt();
           GenParticleSel_eta[i] = genparticle.eta();
           GenParticleSel_phi[i] = genparticle.phi();

       }


   }

   nGenParticle = iGP.size();



   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////// LEPTON CANDIDATES ///////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////

   int e = 0; // index of the electron candidate
   int m = 0; // index of the muon candidate
   float deltaR_min; // variable ot found the minimum deltaR

   std::vector<int> matched_clusters; // clusted that are already matched
   std::vector<int> matched_triggerObjects; // trigger that are already matched

   int m_cluster; // cluster saved iterator
   int m_triggerObject; // trigger saved iterator


   // Loop over the track collection to do a lepton matching
   for (size_t i = 0; i < iT.size(); i++){

       reco::Track track = trackCollection.at(iT.at(i));

       // Matching variables initialization
       deltaR_min = 10.;
       m_cluster = -99;
       m_triggerObject = -99;


       // Electron matching (loop over the photons) REVISIT
       for (size_t jp = 0; jp < iP.size(); jp++){

           const reco::Photon &photon = (*photons)[iP.at(jp)]; // get the photon candidate

           float deltaPhi = fabs(photon.phi() - track.phi());
           float deltaEta = fabs(photon.eta() - track.eta());
           float deltaR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

           if (deltaR < 0.1 && deltaR < deltaR_min){

               m_cluster = jp; // matched cluster
               deltaR_min = deltaR; // deltaR minimization

           }
       }


       // Muon matching (loop over the muon objects) REVISIT
       for (size_t jm = 0; jm < iMT.size(); jm++){

           trigger::TriggerObject muonTriggerObject = triggerObjects[iMT.at(jm)];

           float deltaPhi = fabs(muonTriggerObject.phi() - track.phi());
           float deltaEta = fabs(muonTriggerObject.eta() - track.eta());
           float deltaR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

           if (deltaR < 0.1 && deltaR < deltaR_min){

               m_triggerObject = jm; m_cluster = -99; // matched muon object, unmatch photon object
               deltaR_min = deltaR;

           }
       }

       
       // FINAL LEPTON MATCHING

       // Case 0: No clusters or trigger objects are matched
       if (m_cluster == -99 && m_triggerObject == -99){ continue;}

       // Case 1: Muon candidate found
       if (m_cluster == -99 && m_triggerObject != -99){

           // Check if the muon is already matched:
           if(std::find(matched_triggerObjects.begin(), matched_triggerObjects.end(), m_triggerObject) != matched_triggerObjects.end()){ continue; }

           // If not, define the final muon candidate:
           //trigger::TriggerObject muonTriggerObject = triggerObjects[iMT.at(m_triggerObject)];

           MuonCandidate_pt[m] = track.pt();
           MuonCandidate_phi[m] = track.phi();
           MuonCandidate_eta[m] = track.eta();
           MuonCandidate_muonTriggerObjectIdx[m] = m_triggerObject;
           MuonCandidate_trackIdx[m] = i;

           matched_triggerObjects.push_back(m_triggerObject); // add it to the list of matched muons
           m++; // Next muon candidate

           continue;
       }

       // Case 2: Electron candidate found
       if (m_cluster != -99 && m_triggerObject == -99){

           // Cherk if the cluster is already matched
           if(std::find(matched_clusters.begin(), matched_clusters.end(), m_cluster) != matched_clusters.end()){ continue; }

           // If not, define the final electron candidate
           const reco::Photon photon = (*photons)[iP.at(m_cluster)];

           ElectronCandidate_pt[e] = track.pt();
           ElectronCandidate_et[e] = photon.et();
           ElectronCandidate_phi[e] = track.phi();
           ElectronCandidate_eta[e] = track.eta();
           ElectronCandidate_photonIdx[e] = m_cluster;
           ElectronCandidate_trackIdx[e] = i;

           matched_clusters.push_back(m_cluster); // add it to the list of matched electrons
           e++; // Next electron candidate

       }

   }


   nMuonCandidate = m;
   nElectronCandidate = e;




   /////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////// VERTEX REFITTING /////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////



   // Define a transient track builder: 
   edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);

   // std vector with the tracks participating in the refitting
   std::vector<reco::TransientTrack> refit_tracks; // tracks for refitting

 
   // Loop over the PV tracks:
   for (auto tv=vertex->tracks_begin(); tv!=vertex->tracks_end(); tv++){

       // Transient track definition:
       const reco::TrackRef trackRef = tv->castTo<reco::TrackRef>();
       reco::TransientTrack  transientTrack = theTransientTrackBuilder->build(trackRef); 
       transientTrack.setBeamSpot(beamSpotObject);
  
       // HERE: should be a selection of tracks (ignore those matched to leptons)

       // Append the transient track:
       refit_tracks.push_back(transientTrack);

   }



   // Reffit the vertex with the selected tracks stored in refit_tracks
   if (refit_tracks.size() > 1){

       // AdaptiveVertexFitter definition with a track significance cutoff of 2.5:
       AdaptiveVertexFitter  theFitter(GeometricAnnealing(2.5));

       // Vertex refitting:
       TransientVertex myVertex = theFitter.vertex(refit_tracks);


       // Get the refitted vertex values:
       if (myVertex.isValid()){

           RefittedPV_vx = myVertex.position().x();
           RefittedPV_vy = myVertex.position().y();
           RefittedPV_vz = myVertex.position().z();

       }

   }


   PV_fittingtracks = refit_tracks.size();
   
   
   /////////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////// FILL THE TREE ///////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////
   tree_out->Fill();


}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void AODAnalysis::beginJob()
{

    // Output file definition
    output_filename = parameters.getParameter<std::string>("nameOfOutput");
    file_out = new TFile(output_filename.c_str(), "RECREATE");
    file_out->cd();

    std::cout << "the file is created" << std::endl;
    
    // Output Tree definition
    tree_out = new TTree("Events", "Events");

    std::cout << "The tree is created" << std::endl;

    ///////////////////////////////// EVENT INFO BRANCHES ///////////////////////////////

    tree_out->Branch("Event_event", &Event_event, "Event_event/I");
    tree_out->Branch("Event_run", &Event_run, "Event_run/I");
    tree_out->Branch("Event_luminosityBlock", &Event_luminosityBlock, "Event_luminosityBlock/I");


    ///////////////////////////////// BEAM SPOT BRANCHES ////////////////////////////////

    tree_out->Branch("BeamSpot_x0", &BeamSpot_x0, "BeamSpot_x0/F");
    tree_out->Branch("BeamSpot_y0", &BeamSpot_y0, "BeamSpot_y0/F");
    tree_out->Branch("BeamSpot_z0", &BeamSpot_z0, "BeamSpot_z0/F");
    tree_out->Branch("BeamSpot_BeamWidthX", &BeamSpot_BeamWidthX, "BeamSpot_BeamWidthX/F");
    tree_out->Branch("BeamSpot_BeamWidthY", &BeamSpot_BeamWidthY, "BeamSpot_BeamWidthY/F");


    /////////////////////////////// GENERAL TRACK BRANCHES //////////////////////////////

    tree_out->Branch("nTrack", &nTrack, "nTrack/I");
    tree_out->Branch("TrackSel_pt", TrackSel_pt, "TrackSel_pt[nTrack]/F");
    tree_out->Branch("TrackSel_eta", TrackSel_eta, "TrackSel_eta[nTrack]/F");
    tree_out->Branch("TrackSel_phi", TrackSel_phi, "TrackSel_phi[nTrack]/F");
    tree_out->Branch("TrackSel_dxy", TrackSel_dxy, "TrackSel_dxy[nTrack]/F");
    tree_out->Branch("TrackSel_dxyError", TrackSel_dxyError, "TrackSel_dxyError[nTrack]/F");
    tree_out->Branch("TrackSel_d0", TrackSel_d0, "TrackSel_d0[nTrack]/F");
    tree_out->Branch("TrackSel_d0Error", TrackSel_d0Error, "TrackSel_d0Error[nTrack]/F");
    tree_out->Branch("TrackSel_dz", TrackSel_dz, "TrackSel_dz[nTrack]/F");
    tree_out->Branch("TrackSel_dzError", TrackSel_dzError, "TrackSel_dzError[nTrack]/F");
    tree_out->Branch("TrackSel_vx", TrackSel_vx, "TrackSel_vx[nTrack]/F");
    tree_out->Branch("TrackSel_vy", TrackSel_vy, "TrackSel_vy[nTrack]/F");
    tree_out->Branch("TrackSel_vz", TrackSel_vz, "TrackSel_vz[nTrack]/F");
    tree_out->Branch("TrackSel_isDisplaced", TrackSel_isDisplaced, "TrackSel_isDisplaced[nTrack]/I");
    tree_out->Branch("TrackSel_trackIndex", TrackSel_trackIndex, "TrackSel_trackIndex[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidTrackerHits", TrackSel_numberOfValidTrackerHits, "TrackSel_numberOfValidTrackerHits[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidPixelHits", TrackSel_numberOfValidPixelHits, "TrackSel_numberOfValidPixelHits[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidPixelBarrelHits", TrackSel_numberOfValidPixelBarrelHits, "TrackSel_numberOfValidPixelBarrelHits[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidPixelEndcapHits", TrackSel_numberOfValidPixelEndcapHits, "TrackSel_numberOfValidPixelEndcapHits[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidStripHits", TrackSel_numberOfValidStripHits, "TrackSel_numberOfValidStripHits[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidStripTIBHits", TrackSel_numberOfValidStripTIBHits, "TrackSel_numberOfValidStripTIBHits[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidStripTIDHits", TrackSel_numberOfValidStripTIDHits, "TrackSel_numberOfValidStripTIDHits[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidStripTOBHits", TrackSel_numberOfValidStripTOBHits, "TrackSel_numberOfValidStripTOBHits[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidStripTECHits", TrackSel_numberOfValidStripTECHits, "TrackSel_numberOfValidStripTECHits[nTrack]/I");



    ////////////////////////////// PRIMARY VERTEX BRANCHES //////////////////////////////

    tree_out->Branch("nPV", &nPV, "nPV/I");
    tree_out->Branch("PV_vx", &PV_vx, "PV_vx/F");
    tree_out->Branch("PV_vy", &PV_vy, "PV_vy/F");
    tree_out->Branch("PV_vz", &PV_vz, "PV_vz/F");
    tree_out->Branch("PV_nTracks", &PV_nTracks, "PV_nTracks/I");
    tree_out->Branch("PV_fittingtracks", &PV_fittingtracks, "PV_fittingtracks/I");
    tree_out->Branch("PV_tracksSize", &PV_tracksSize, "PV_tracksSize/I");


    ////////////////////////////////// PHOTON BRANCHES //////////////////////////////////

    tree_out->Branch("nPhoton", &nPhoton, "nPhoton/I");
    tree_out->Branch("nPhotonOriginal", &nPhotonOriginal, "nPhotonOirginal/I");
    tree_out->Branch("PhotonSel_et", PhotonSel_et, "PhotonSel_et[nPhoton]/F");
    tree_out->Branch("PhotonSel_eta", PhotonSel_eta, "PhotonSel_eta[nPhoton]/F");
    tree_out->Branch("PhotonSel_phi", PhotonSel_phi, "PhotonSel_phi[nPhoton]/F");
    tree_out->Branch("PhotonSel_hadronicOverEm", PhotonSel_hadronicOverEm, "PhotonSel_hadronicOverEm[nPhoton]/F");
    tree_out->Branch("PhotonSel_full5x5_sigmaIetaIeta", PhotonSel_full5x5_sigmaIetaIeta, "PhotonSel_full5x5_sigmaIetaIeta[nPhoton]/F");
    tree_out->Branch("PhotonSel_isEB", PhotonSel_isEB, "PhotonSel_isEB[nPhoton]/I");
    tree_out->Branch("PhotonSel_isEE", PhotonSel_isEE, "PhotonSel_isEE[nPhoton]/I");


    //////////////////////////////// GENPARTICLE BRANCHES ///////////////////////////////

    tree_out->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
    tree_out->Branch("GenParticleSel_pt", GenParticleSel_pt, "GenParticleSel_pt[nGenParticle]/F");
    tree_out->Branch("GenParticleSel_eta", GenParticleSel_eta, "GenParticleSel_eta[nGenParticle]/F");
    tree_out->Branch("GenParticleSel_phi", GenParticleSel_phi, "GenParticleSel_phi[nGenParticle]/F");
    tree_out->Branch("GenParticleSel_pdgId", GenParticleSel_pdgId, "GenParticleSel_pdgId[nGenParticle]/I"); 


    //////////////////////////// MUON TRIGGER OBJECT BRANCHES ///////////////////////////

    tree_out->Branch("nMuonTriggerObject", &nMuonTriggerObject, "nMuonTriggerObject/I");
    tree_out->Branch("MuonTriggerObjectSel_pt", MuonTriggerObjectSel_pt, "MuonTriggerObjectSel_pt[nMuonTriggerObject]/F");
    tree_out->Branch("MuonTriggerObjectSel_eta", MuonTriggerObjectSel_eta, "MuonTriggerObjectSel_eta[nMuonTriggerObject]/F");
    tree_out->Branch("MuonTriggerObjectSel_phi", MuonTriggerObjectSel_phi, "MuonTriggerObjectSel_phi[nMuonTriggerObject]/F");



    ////////////////////////////////// TRIGGER INFORMATION //////////////////////////////

    tree_out->Branch("HLT_DoubleMu43NoFiltersNoVtx_v3", &HLT_DoubleMu43NoFiltersNoVtx_v3, "HLT_DoubleMu43NoFiltersNoVtx_v3/O");


    //////////////////////////// ELECTRON CANDIDATE BRANCHES ////////////////////////////

    tree_out->Branch("nElectronCandidate", &nElectronCandidate, "nElectronCandidate/I");
    tree_out->Branch("ElectronCandidate_pt", ElectronCandidate_pt, "ElectronCandidate_pt[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_et", ElectronCandidate_et, "ElectronCandidate_et[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_eta", ElectronCandidate_eta, "ElectronCandidate_eta[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_phi", ElectronCandidate_phi, "ElectronCandidate_phi[nElectronCandidate]/F");
    tree_out->Branch("ElectronCandidate_photonIdx", ElectronCandidate_photonIdx, "ElectronCandidate_photonIdx[nElectronCandidate]/I");
    tree_out->Branch("ElectronCandidate_trackIdx", ElectronCandidate_trackIdx, "ElectronCandidate_trackIdx[nElectronCandidate]/I");


    ////////////////////////////// MUON CANDIDATE BRANCHES /////////////////////////////

    tree_out->Branch("nMuonCandidate", &nMuonCandidate, "nMuonCandidate/I");
    tree_out->Branch("MuonCandidate_pt", MuonCandidate_pt, "MuonCandidate_pt[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_eta", MuonCandidate_eta, "MuonCandidate_eta[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_phi", MuonCandidate_phi, "MuonCandidate_phi[nMuonCandidate]/F");
    tree_out->Branch("MuonCandidate_muonTriggerObjectIdx", MuonCandidate_muonTriggerObjectIdx, "MuonCandidate_muonTriggerObjectIdx[nMuonCandidate]/I");
    tree_out->Branch("MuonCandidate_trackIdx", MuonCandidate_trackIdx, "MuonCandidate_trackIdx[nMuonCandidate]/I");



    /////////////////////////// REFITTED PRIMARY VERTEX BRANCHES ////////////////////////

    tree_out->Branch("RefittedPV_vx", &RefittedPV_vx, "RefittedPV_vx/F");
    tree_out->Branch("RefittedPV_vy", &RefittedPV_vy, "RefittedPV_vy/F");
    tree_out->Branch("RefittedPV_vz", &RefittedPV_vz, "RefittedPV_vz/F");






}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void AODAnalysis::endJob() 
{


    std::cout << "The event is writen" << std::endl;
    file_out->cd();
    tree_out->Write();
    file_out->Close();

}
//=======================================================================================================================================================================================================================//

/*
void AODAnalysis::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,"HLT",changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
    // The HLT config has actually changed wrt the previous Run, hence rebook your
    //  histograms or do anything else dependent on the revised HLT config

    const std::vector<std::string>& filters = hltConfig_.moduleLabels("HLT_DoubleMu43NoFiltersNoVtx_v3");

    std::cout << filters[0] << std::endl;

    }
  } else {
    //LogError("MyAnalyzer") << " HLT config extraction failure with process name " << processName_;

  }

}
*/

//=======================================================================================================================================================================================================================//
void AODAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//=======================================================================================================================================================================================================================//





DEFINE_FWK_MODULE(AODAnalysis);
