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
//#include "DataFormats/PatCandidates/interface/Photon.h"
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


//-> GENERAL TRACKS
const Int_t nTrackMax = 1000;
Int_t nTrack;
Int_t nTrackOriginal;
Float_t TrackSel_pt[nTrackMax];
Float_t TrackSel_eta[nTrackMax];
Float_t TrackSel_phi[nTrackMax];
Float_t TrackSel_dxy[nTrackMax];
Float_t TrackSel_dxyError[nTrackMax];
Float_t TrackSel_dz[nTrackMax];
Float_t TrackSel_dzError[nTrackMax];
Float_t TrackSel_vx[nTrackMax];
Float_t TrackSel_vy[nTrackMax];
Float_t TrackSel_vz[nTrackMax];
Int_t TrackSel_numberOfValidTrackerHits[nTrackMax];
Int_t TrackSel_numberOfValidPixelHits[nTrackMax];


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
   edm::Handle<reco::BeamSpot> beamSpot;
   edm::Handle<edm::View<reco::Photon> > photons;
   edm::Handle<edm::View<reco::GenParticle> > genParticles;

   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<trigger::TriggerEvent> trigSummary;


   //-> Get by Token from the file
   iEvent.getByToken(thePrimaryVertexCollection, primaryVertices);
   iEvent.getByToken(theGeneralTrackCollection, generalTracks);
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


   ////////////////////////////////// GENERAL TRACKS //////////////////////////////////

   std::vector<int> iT; // track selection indexes

   // Loop over general tracks
   for (size_t i = 0; i < generalTracks->size(); i++){

       const reco::Track & track = (*generalTracks)[i];

       // Check if the track fulfils the preselection requirements
       if (generalTrackPreselection(track)) { iT.push_back(i); }

   }


   // Sort preselected track indexes by pt
   std::sort( std::begin(iT), std::end(iT), [&](int i1, int i2){ return generalTracks->at(i1).pt() < generalTracks->at(i2).pt(); });


   // Loop over the preselected tracks
   for (size_t i = 0; i < iT.size(); i++){

       const reco::Track & track = (*generalTracks)[iT.at(i)];
       
       // Get the variables
       TrackSel_pt[i] = track.pt();
       TrackSel_eta[i] = track.eta();
       TrackSel_phi[i] = track.phi();
       TrackSel_dxy[i] = track.dxy();
       TrackSel_dxyError[i] = track.dxyError();
       TrackSel_dz[i] = track.dz();
       TrackSel_dzError[i] = track.dzError();
       TrackSel_vx[i] = track.vx();
       TrackSel_vy[i] = track.vy();
       TrackSel_vz[i] = track.vz(); 

       // Hit info
       const reco::HitPattern &hits = track.hitPattern();

       TrackSel_numberOfValidTrackerHits[i] = hits.numberOfValidTrackerHits();
       TrackSel_numberOfValidPixelHits[i] = hits.numberOfValidPixelHits();

   }


   nTrack = iT.size(); // number of selected tracks
   nTrackOriginal = generalTracks->size(); // original number of general tracks in AOD


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

       trigger::TriggerObject foundObject = triggerObjects[i];

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


   std::cout << "Original PV: " << PV_vx << "\t" << "Refitted_PV: " << RefittedPV_vx << std::endl;
   
   
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
    tree_out->Branch("nTrackOriginal", &nTrackOriginal, "nTrackOriginal/I");
    tree_out->Branch("TrackSel_pt", TrackSel_pt, "TrackSel_pt[nTrack]/F");
    tree_out->Branch("TrackSel_eta", TrackSel_eta, "TrackSel_eta[nTrack]/F");
    tree_out->Branch("TrackSel_phi", TrackSel_phi, "TrackSel_phi[nTrack]/F");
    tree_out->Branch("TrackSel_dxy", TrackSel_dxy, "TrackSel_dxy[nTrack]/F");
    tree_out->Branch("TrackSel_dxyError", TrackSel_dxyError, "TrackSel_dxyError[nTrack]/F");
    tree_out->Branch("TrackSel_dz", TrackSel_dz, "TrackSel_dz[nTrack]/F");
    tree_out->Branch("TrackSel_dzError", TrackSel_dzError, "TrackSel_dzError[nTrack]/F");
    tree_out->Branch("TrackSel_vx", TrackSel_vx, "TrackSel_vx[nTrack]/F");
    tree_out->Branch("TrackSel_vy", TrackSel_vy, "TrackSel_vy[nTrack]/F");
    tree_out->Branch("TrackSel_vz", TrackSel_vz, "TrackSel_vz[nTrack]/F");
    tree_out->Branch("TrackSel_numberOfValidTrackerHits", TrackSel_numberOfValidTrackerHits, "TrackSel_numberOfValidTrackerHits[nTrack]/I");
    tree_out->Branch("TrackSel_numberOfValidPixelHits", TrackSel_numberOfValidPixelHits, "TrackSel_numberOfValidPixelHits[nTrack]/I");



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
