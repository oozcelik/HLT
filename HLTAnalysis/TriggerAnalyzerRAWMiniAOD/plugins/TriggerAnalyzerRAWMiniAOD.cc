// -*- C++ -*-
//
// Package:    HLTAnalysis/TriggerAnalyzerRAWMiniAOD
// Class:      TriggerAnalyzerRAWMiniAOD
// 
/**\class TriggerAnalyzerRAWMiniAOD TriggerAnalyzerRAWMiniAOD.cc HLTAnalysis/TriggerAnalyzerRAWMiniAOD/plugins/TriggerAnalyzerRAWMiniAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Laurent Thomas
//         Created:  Fri, 24 Mar 2017 04:09:55 GMT
//
//


// system include files
#include <memory>
#include "TEfficiency.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"


#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TriggerAnalyzerRAWMiniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet&);
      ~TriggerAnalyzerRAWMiniAOD();
  

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      bool PassOfflineMuonSelection(const reco::Muon *mu, reco::Vertex::Point PV);
      bool RecoHLTMatching(const edm::Event&,double recoeta, double recophi, std::string filtername, double dRmatching = 0.3);
      unsigned int runNum, evtNum, lumiNum;
      vector<float> *JMass;
      TTree* X_One_Tree_;
      // ----------member data ---------------------------


  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectsMINIAODToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsORIGToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsRAWToken_;
  edm::EDGetTokenT<edm::TriggerResults>  trgresultsHLT2Token_;

  edm::EDGetTokenT<std::vector<pat::Jet> > jet_token;
  edm::EDGetTokenT<std::vector<reco::Muon> > muon_token;
  edm::EDGetTokenT<std::vector<pat::Electron> > electron_token;
  edm::EDGetTokenT<reco::VertexCollection>  PV_token;

  edm::EDGetTokenT<std::vector<reco::BeamSpot>> beam_token_;

  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> et_Filter_Token_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> showershape_Filter_Token_;
  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> dphi_Filter_Token_;

  edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > showershape_Var_Token_;
  edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > hovere_Var_Token_;
  edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > trackiso_Var_Token_;

  edm::Service<TFileService> fs;

  TH1F* h_mu3pfjet200csv1p5_vs_leadbjetpt_den;
  TH1F* h_mu3pfjet200csv1p5_vs_leadbjetpt_num;
  TH1F* h_mu3pfjet200csv1p5_vs_leadbjetpt_numl1;
  TH1F* h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_den;
  TH1F* h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_num;
  TH1F* h_mu3pfjet200csv1p5_vs_nbjetspt200_den;
  TH1F* h_mu3pfjet200csv1p5_vs_nbjetspt200_num;
  TH1F* h_DoubleMu43_vs_leadingmuonpt_den;
  TH1F* h_DoubleMu43_vs_leadingmuonpt_num;
  TH1F* h_ele35wptight_lastfilter_den;
  TH1F* h_ele35wptight_lastfilter_num;
  TH1F* h_sietaieta_HLT;
  TH1F* h_hoe_HLT;
  TH1F* h_trackiso_HLT;
  TEfficiency *pEff;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerAnalyzerRAWMiniAOD::TriggerAnalyzerRAWMiniAOD(const edm::ParameterSet& iConfig)

{
  trigobjectsMINIAODToken_ = consumes<pat::TriggerObjectStandAloneCollection>( edm::InputTag("patTrigger")); //selectedPatTrigger
  trigobjectsRAWToken_=consumes<trigger::TriggerEvent>(edm::InputTag("hltTriggerSummaryAOD", "","MYHLT"));  
  trgresultsORIGToken_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults","","HLT") );
  trgresultsHLT2Token_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults","","MYHLT") );
  muon_token = consumes<std::vector<reco::Muon >>(edm::InputTag("muons") );
  PV_token = consumes<reco::VertexCollection>  (edm::InputTag("offlinePrimaryVertices"));
  beam_token_ = consumes<std::vector<reco::BeamSpot> > (edm::InputTag("offlineBeamSpot"));

  //now do what ever initialization is needed
  //   usesResource("TFileService");

  h_mu3pfjet200csv1p5_vs_leadbjetpt_den= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_leadbjetpt_den","",50,0,500);
  h_mu3pfjet200csv1p5_vs_leadbjetpt_num= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_leadbjetpt_num","",50,0,500);
  h_mu3pfjet200csv1p5_vs_leadbjetpt_numl1= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_leadbjetpt_numl1","",50,0,500);
  h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_den= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_den","",101,0,1.01);
  h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_num= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_num","",101,0,1.01);
  h_mu3pfjet200csv1p5_vs_nbjetspt200_den= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_nbjetspt200_den","",5,0,5);
  h_mu3pfjet200csv1p5_vs_nbjetspt200_num= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_nbjetspt200_num","",5,0,5);
  h_DoubleMu43_vs_leadingmuonpt_den= fs->make<TH1F>("h_DoubleMu43_vs_leadingmuonpt_den","",100,0,100);
  h_DoubleMu43_vs_leadingmuonpt_num= fs->make<TH1F>("h_DoubleMu43_vs_leadingmuonpt_num","",100,0,100);
  h_ele35wptight_lastfilter_den= fs->make<TH1F>("h_ele35wptight_lastfilter_den","",20,0,100);
  h_ele35wptight_lastfilter_num= fs->make<TH1F>("h_ele35wptight_lastfilter_num","",20,0,100);

  h_sietaieta_HLT= fs->make<TH1F>("h_sietaieta_HLT","",100,0,0.05);
  h_hoe_HLT= fs->make<TH1F>("h_hoe_HLT","",100,0,0.2);
  h_trackiso_HLT= fs->make<TH1F>("h_trackiso_HLT","",100,0,0.5);
  pEff = fs->make<TEfficiency >("eff","my efficiency;x;#epsilon",100,0,100);

}


TriggerAnalyzerRAWMiniAOD::~TriggerAnalyzerRAWMiniAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerAnalyzerRAWMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace reco;
   using namespace std;


  runNum = iEvent.id().run();
  evtNum = iEvent.id().event();
  lumiNum = iEvent.id().luminosityBlock();


//    ****************Part 1. Accessing some trigger information ************* 
   bool passHLT_DoubleMu43(false);

   //Accessing trigger bits:
   //This works in both RAW, AOD or MINIAOD 
   //Here we access the decision provided by the HLT (i.e. original trigger step). 
   edm::Handle<edm::TriggerResults> trigResults;
   iEvent.getByToken(trgresultsHLT2Token_, trigResults);

   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
   
   ESHandle < MagneticField > bFieldHandle;
   iSetup.get < IdealMagneticFieldRecord > ().get(bFieldHandle);

   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByToken(beam_token_, beamSpotHandle);

   reco::BeamSpot beamSpot;
   Vertex theBeamSpotV;

   if (beamSpotHandle.isValid())
    {
        beamSpot = *beamSpotHandle;
        theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }

     
if( !trigResults.failedToGet() ) {
     int N_Triggers = trigResults->size();
     const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
     for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {

       if (trigResults.product()->accept(i_Trig)) {
	 TString TrigPath =trigName.triggerName(i_Trig);
	 if(TrigPath.Index("HLT_DoubleMu4_3_Jpsi") >=0)  passHLT_DoubleMu43 = true;
       }
     }
    }
 
   //Accessing the trigger objects in RAW/AOD
   edm::Handle<trigger::TriggerEvent> triggerObjectsSummary;
   iEvent.getByToken(trigobjectsRAWToken_ ,triggerObjectsSummary);
   trigger::TriggerObjectCollection selectedObjects;
   bool passedFilter = false;
   if (triggerObjectsSummary.isValid()) {
     size_t filterIndex = (*triggerObjectsSummary).filterIndex( edm::InputTag("hltmumuFilterDoubleMu43Jpsi","","MYHLT") );
     trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();
     if (filterIndex < (*triggerObjectsSummary).sizeFilters()) { 
       const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
       for (size_t j = 0; j < keys.size(); j++) {
	 trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];	
	 passedFilter=true;	
         h_mu3pfjet200csv1p5_vs_leadbjetpt_den->Fill(foundObject.pt());
	// cout <<"object found, printing pt, eta, phi, id: " <<filterIndex << " " << foundObject.pt()<<", "<<foundObject.eta()<<", "<< foundObject.phi() << ", " << foundObject.id() << endl;
       }
     }
   }
   

   //Offline muons and PVs for AODSIM
   edm::Handle< std::vector<reco::Muon> > muons;
   iEvent.getByToken(muon_token,muons );
   Handle<reco::VertexCollection> theVertices;
   iEvent.getByToken(PV_token,theVertices);
   int nvertex = theVertices->size();
   Vertex::Point PV(0,0,0);
   Vertex thePrimaryV;
   if( nvertex>0) {
	   PV = theVertices->begin()->position();
	   thePrimaryV = Vertex(*(theVertices->begin()));
   }
   //Count the nb of offline muons with pt >3
   int nmuonspt3 =0;
   double leadingmuonpt(-10);
   int npassed(0);
   bool passed(false);
   for( std::vector<reco::Muon>::const_iterator mu1 = (muons)->begin(); mu1 != (muons)->end(); mu1++ ) {
           if(!(muon::isSoftMuon(*mu1, thePrimaryV))) continue;
	   if( (mu1->pt() < 3. || mu1->pt() >4.) || fabs(mu1->eta()) >2.5 ) continue;

	   for(std::vector<reco::Muon>::const_iterator mu2 = mu1 + 1; mu2 != (muons)->end(); mu2++){
		   if(!( mu1->charge() * mu2->charge() < 0 )) continue; 
		   if( mu2 == mu1 ) continue;	
		   if(!(muon::isSoftMuon(*mu2, thePrimaryV))) continue;
		   if( (mu2->pt() < 3. || mu2->pt() >4.) || fabs(mu2->eta()) >2.5 ) continue;

		   /// Fit the mu1 & mu2 to a common vertex and make J/Psi
		   std::vector<reco::TransientTrack> j_tks;
		   j_tks.push_back((*theB).build(mu1->track().get()));
		   j_tks.push_back((*theB).build(mu2->track().get()));


		   KalmanVertexFitter jkvf;
		   TransientVertex jtv = jkvf.vertex(j_tks);
		   if (!jtv.isValid()) continue;

		   reco::Vertex jpsivertex = jtv;
                   double dimuonCL(0);
		   if( (jpsivertex.chi2()>=0.0) && (jpsivertex.ndof()>0) ) dimuonCL = TMath::Prob(jpsivertex.chi2(), jpsivertex.ndof() );
		   if(dimuonCL<0.1) continue;
		   math::XYZVector jpperp(mu1->px() + mu2->px() , mu1->py() + mu2->py(), 0.);
		   GlobalPoint jVertex = jtv.position();
		   GlobalError jerr    = jtv.positionError();

		   GlobalPoint displacementFromBeamspot (-1*((beamSpot.x0() - jVertex.x()) + (jVertex.z() - beamSpot.z0()) * beamSpot.dxdz()), -1*((beamSpot.y0() - jVertex.y()) + (jVertex.z() - beamSpot.z0()) * beamSpot.dydz()), 0);
     
                  reco::Vertex::Point vperp(displacementFromBeamspot.x(), displacementFromBeamspot.y(),  0.);
                  float cosAlpha = vperp.Dot(jpperp) / (vperp.R() * jpperp.R());
                  if (cosAlpha<0.9) continue;

		  reco::Particle::LorentzVector MuMuP4_1, MuMuP4_2, JpsiP4;
		  ParticleMass muon_mass = 0.10565837;

		  MuMuP4_1 = reco::Particle::LorentzVector(mu1->px() , mu1->py() , mu1->pz() , muon_mass );
		  MuMuP4_2 = reco::Particle::LorentzVector(mu2->px() , mu2->py() , mu2->pz() , muon_mass );

		  JpsiP4 = MuMuP4_1 + MuMuP4_2;

		  if(JpsiP4.M()<2.9 || JpsiP4.M()>3.3) continue;
		  if(JpsiP4.Pt()<6.9) continue;
                  leadingmuonpt = JpsiP4.Pt();
		  npassed++;
		  passed = true;
		  if(passed) JMass->push_back(JpsiP4.M());
	   }

   }    

   //Effcy vs leading muon pt:
   if(passed) h_DoubleMu43_vs_leadingmuonpt_den->Fill(leadingmuonpt); ///all events that pass the offline selections.
   if(passHLT_DoubleMu43 && passedFilter) h_DoubleMu43_vs_leadingmuonpt_num->Fill(leadingmuonpt); //only events pass my HLT path.
   pEff->Fill(passedFilter, leadingmuonpt);
}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::beginJob()
{
	  X_One_Tree_ = fs->make < TTree > ("HLT", "Efficiency");
	  
	  X_One_Tree_->Branch("evtNum", &evtNum, "evtNum/i");
	  X_One_Tree_->Branch("runNum", &runNum, "runNum/i");
	  X_One_Tree_->Branch("lumiNum", &lumiNum, "lumiNum/i");
	  X_One_Tree_->Branch("JMass", &JMass);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::endJob() 
{
	X_One_Tree_->GetDirectory()->cd();
	X_One_Tree_->Write();
}
//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzerRAWMiniAOD);
