//=======================================================================================================================================================================================================================// 
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
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Candidate/interface/Candidate.h"


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
Int_t nOriginalTrack;
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



   //-> Get by Token from the file
   iEvent.getByToken(thePrimaryVertexCollection, primaryVertices);
   iEvent.getByToken(theGeneralTrackCollection, generalTracks);
   iEvent.getByToken(theBeamSpot, beamSpot);
   iEvent.getByToken(thePhotonCollection, photons);


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
   int itrack = 0; // track counter

   // Loop over general tracks
   for (size_t i = 0; i < generalTracks->size(); i++){

       const reco::Track & track = (*generalTracks)[i];

       // Check if the track fulfils the preselection requirements
       if (!generalTrackPreselection(track)) { continue; }

       iT.push_back(i);

       
       // Get the variables
       TrackSel_pt[itrack] = track.pt();
       TrackSel_eta[itrack] = track.eta();
       TrackSel_phi[itrack] = track.phi();
       TrackSel_dxy[itrack] = track.dxy();
       TrackSel_dxyError[itrack] = track.dxyError();
       TrackSel_dz[itrack] = track.dz();
       TrackSel_dzError[itrack] = track.dzError();
       TrackSel_vx[itrack] = track.vx();
       TrackSel_vy[itrack] = track.vy();
       TrackSel_vz[itrack] = track.vz(); 

       // Hit info
       const reco::HitPattern &hits = track.hitPattern();

       TrackSel_numberOfValidTrackerHits[i] = hits.numberOfValidTrackerHits();
       TrackSel_numberOfValidPixelHits[i] = hits.numberOfValidPixelHits();

       // Update the counter
       itrack++;

   }


   nTrack = iT.size(); // number of selected tracks
   nOriginalTrack = generalTracks->size(); // original number of general tracks in AOD


   ////////////////////////////// PRIMARY VERTEX FEATURES //////////////////////////////


   nPV = primaryVertices->size();
   const reco::Vertex *vertex = &((*primaryVertices)[0]);

   PV_vx = vertex->position().x();
   PV_vy = vertex->position().y();
   PV_vz = vertex->position().z();


   PV_nTracks = vertex->nTracks(); // to check
   PV_tracksSize= vertex->tracksSize(); // to check

 
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
    tree_out->Branch("nOriginalTrack", &nOriginalTrack, "nOriginalTrack/I");
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
