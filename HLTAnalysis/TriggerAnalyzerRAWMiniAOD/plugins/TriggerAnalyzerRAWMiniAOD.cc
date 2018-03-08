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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"
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

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  bool PassOfflineMuonSelection(const pat::Muon *mu, reco::Vertex::Point PV);
  bool PassOfflineElectronSelection(const pat::Electron * ele, reco::Vertex::Point PV);
  bool RecoHLTMatching(const edm::Event&,double recoeta, double recophi, std::string filtername, double dRmatching = 0.3);
  double VarStudied( const edm::Event& iEvent, double recoeta, double recophi,edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > varToken_,  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> candToken_,   bool  dividebyE, bool dividebyEt, double dRmatching =0.3);

      // ----------member data ---------------------------


  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectsMINIAODToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsORIGToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsRAWToken_;
  edm::EDGetTokenT<edm::TriggerResults>  trgresultsHLT2Token_;

  edm::EDGetTokenT<std::vector<pat::Jet> > jet_token;
  edm::EDGetTokenT<std::vector<pat::Muon> > muon_token;
  edm::EDGetTokenT<std::vector<pat::Electron> > electron_token;
  edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;

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
  TH1F* h_mu3pfjet200csv1p5_vs_leadingmuonpt_den;
  TH1F* h_mu3pfjet200csv1p5_vs_leadingmuonpt_num;
  TH1F* h_ele35wptight_lastfilter_den;
  TH1F* h_ele35wptight_lastfilter_num;
  TH1F* h_sietaieta_HLT;
  TH1F* h_hoe_HLT;
  TH1F* h_trackiso_HLT;
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
  trigobjectsRAWToken_=consumes<trigger::TriggerEvent>(edm::InputTag("hltTriggerSummaryAOD::MYHLT"));  

  trgresultsORIGToken_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );
  trgresultsHLT2Token_= consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::MYHLT") );


  showershape_Var_Token_  = consumes<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > ( edm::InputTag("hltEgammaClusterShape","sigmaIEtaIEta5x5","HLT2") );
  hovere_Var_Token_  = consumes<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > ( edm::InputTag("hltEgammaHoverE","","HLT2") );
  trackiso_Var_Token_  = consumes<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > ( edm::InputTag("hltEgammaEleGsfTrackIso","","HLT2")  );

  et_Filter_Token_ = consumes<trigger::TriggerFilterObjectWithRefs> ( edm::InputTag("hltEG35L1SingleEGOrEtFilter","","HLT2") ) ;
  showershape_Filter_Token_ = consumes<trigger::TriggerFilterObjectWithRefs> ( edm::InputTag("hltEle35noerWPTightClusterShapeFilter","","HLT2") );
  dphi_Filter_Token_ = consumes<trigger::TriggerFilterObjectWithRefs> ( edm::InputTag("hltEle35noerWPTightGsfDphiFilter","","HLT2") );
  

  jet_token = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJets") );
  muon_token = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons") );
  electron_token = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons") );
  PV_token = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));
  

  //now do what ever initialization is needed
  //   usesResource("TFileService");

  h_mu3pfjet200csv1p5_vs_leadbjetpt_den= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_leadbjetpt_den","",50,0,500);
  h_mu3pfjet200csv1p5_vs_leadbjetpt_num= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_leadbjetpt_num","",50,0,500);
  h_mu3pfjet200csv1p5_vs_leadbjetpt_numl1= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_leadbjetpt_numl1","",50,0,500);
  h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_den= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_den","",101,0,1.01);
  h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_num= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_highestcsv_jetpt250_num","",101,0,1.01);
  h_mu3pfjet200csv1p5_vs_nbjetspt200_den= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_nbjetspt200_den","",5,0,5);
  h_mu3pfjet200csv1p5_vs_nbjetspt200_num= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_nbjetspt200_num","",5,0,5);
  h_mu3pfjet200csv1p5_vs_leadingmuonpt_den= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_leadingmuonpt_den","",100,0,100);
  h_mu3pfjet200csv1p5_vs_leadingmuonpt_num= fs->make<TH1F>("h_mu3pfjet200csv1p5_vs_leadingmuonpt_num","",100,0,100);
  h_ele35wptight_lastfilter_den= fs->make<TH1F>("h_ele35wptight_lastfilter_den","",20,0,100);
  h_ele35wptight_lastfilter_num= fs->make<TH1F>("h_ele35wptight_lastfilter_num","",20,0,100);

  h_sietaieta_HLT= fs->make<TH1F>("h_sietaieta_HLT","",100,0,0.05);
  h_hoe_HLT= fs->make<TH1F>("h_hoe_HLT","",100,0,0.2);
  h_trackiso_HLT= fs->make<TH1F>("h_trackiso_HLT","",100,0,0.5);


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



   // ****************Part 1. Accessing some trigger information ************* 
   bool passHLT_IsoMu27(false);
   bool passHLT_Mu3_PFJet200CSV_1p5(false), passHLT_Mu3_L1SingleJet180(false), passHLT_PFJet200CSV_1p5(false);   

   //Accessing trigger bits:
   //This works in both RAW, AOD or MINIAOD 
   //Here we access the decision provided by the HLT (i.e. original trigger step). 
   edm::Handle<edm::TriggerResults> trigResults;
   iEvent.getByToken(trgresultsORIGToken_, trigResults);
   if( !trigResults.failedToGet() ) {
     int N_Triggers = trigResults->size();
     const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);

     for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
       if (trigResults.product()->accept(i_Trig)) {
	 TString TrigPath =trigName.triggerName(i_Trig);
	 	 cout << "Passed path: " << TrigPath<<endl;
//	 if(TrigPath.Index("HLT_DoubleMu4_3_JPsi_Displaced_v") >=0) passHLT_IsoMu27=true;          
	 //Notice the special syntax: since the path version can change during data taking one only looks for the string "HLT_IsoMu27_v"
       }
     }
   }
   //Exercise 1: 
   //Clone and *then* modify the code above in order to save the decision of your customized HLT menu in the booleans passHLT_Mu3_PFJet200CSV_1p5, passHLT_Mu3_L1SingleJet180, passHLT_PFJet200CSV_1p5
   //Do not directly edit the code above as you will also need the use the original HLT_IsoMu27 decision later on.

   



   //Accessing the trigger objects in MINIAOD
   //This recipe works for MINIAOD only


   //Exercise 2: uncomment the lines above to print all the trigger objects and their corresponding pt, eta, phi. 
   

   //Accessing the trigger objects in RAW/AOD
   //Printing here all trigger objects corresponding to the filter hltL3MuFiltered3
   edm::Handle<trigger::TriggerEvent> triggerObjectsSummary;
   iEvent.getByToken(trigobjectsRAWToken_ ,triggerObjectsSummary);
   trigger::TriggerObjectCollection selectedObjects;
   if (triggerObjectsSummary.isValid()) {
     size_t filterIndex = (*triggerObjectsSummary).filterIndex( edm::InputTag("hltmumuFilterDoubleMu43Jpsi","","MYHLT") );
     trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();
     if (filterIndex < (*triggerObjectsSummary).sizeFilters()) { 
       const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
       for (size_t j = 0; j < keys.size(); j++) {
	 trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	 cout <<"object found, printing pt, eta, phi: " <<foundObject.pt()<<", "<<foundObject.eta()<<", "<< foundObject.phi() <<endl;
       }
     }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnalyzerRAWMiniAOD::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerAnalyzerRAWMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool TriggerAnalyzerRAWMiniAOD::PassOfflineMuonSelection(const pat::Muon *mu, reco::Vertex::Point PV){
  if ( !(mu->isGlobalMuon() || mu->isTrackerMuon() )) return false;
  if ( !(mu->isPFMuon()) ) return false;
  const reco::TrackRef innerTrack = mu->innerTrack();
  if( innerTrack.isNull() )return false;
  
  bool goodGlb =  mu->isGlobalMuon() &&  mu->globalTrack()->normalizedChi2() < 3
    &&  mu->combinedQuality().chi2LocalPosition < 12  && mu->combinedQuality().trkKink < 20;
  bool good =  mu->innerTrack()->validFraction() >= 0.8 &&  mu->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451)  ;
  
  if(!good) return false;
  if(TMath::Abs(innerTrack->dxy(PV)) >0.1 ) return false;
  if(TMath::Abs(innerTrack->dz(PV)) >0.1 ) return false;  
  
  double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
  double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
  double photonIso = mu->pfIsolationR03().sumPhotonEt;
  
  double beta = mu->pfIsolationR03().sumPUPt;
  double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
  
  if(pfRelIsoMu >0.4) return false;
  return true;
}









bool TriggerAnalyzerRAWMiniAOD::PassOfflineElectronSelection(const pat::Electron * ele, reco::Vertex::Point PV){
  const reco::GsfTrackRef gsfTrack = ele->gsfTrack();
  if (!gsfTrack.isNonnull()) return false;
  if( TMath::Abs(gsfTrack->dxy(PV)) > 0.05  )  return false;
  if( TMath::Abs(gsfTrack->dz(PV)) > 0.2  )  return false;
  if(TMath::Abs(ele->superCluster()->eta()) >2.5) return false;
  else if( TMath::Abs(ele->superCluster()->eta()) < 1.479  ) {
    if( TMath::Abs(ele->full5x5_sigmaIetaIeta()) > 0.0103 ) return  false;
    if( TMath::Abs(ele->deltaEtaSuperClusterTrackAtVtx()) > 0.0105 ) return  false;
    if( TMath::Abs(ele->deltaPhiSuperClusterTrackAtVtx()) > 0.115  ) return  false;
    if( TMath::Abs(ele->hadronicOverEm())  > 0.104  ) return  false;
    if( TMath::Abs(1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy() )>0.102 ) return  false;
    if( TMath::Abs(gsfTrack->dxy(PV)) > 0.0261) return false;
    if( TMath::Abs(gsfTrack->dz(PV)) > 0.41) return false;
  }
  else {
    if( TMath::Abs(ele->full5x5_sigmaIetaIeta()) > 0.0301 ) return  false;
    if( TMath::Abs(ele->deltaEtaSuperClusterTrackAtVtx()) > 0.00814 ) return  false;
    if( TMath::Abs(ele->deltaPhiSuperClusterTrackAtVtx()) > 0.182  ) return  false;
    if( TMath::Abs(ele->hadronicOverEm())  > 0.0897  ) return  false;
    if( TMath::Abs(1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy() )>0.126 ) return  false;
    if( TMath::Abs(gsfTrack->dxy(PV)) > 0.0118) return false;
    if( TMath::Abs(gsfTrack->dz(PV)) > 0.822) return false;


  }

  double iso = (ele->pfIsolationVariables().sumChargedHadronPt 
		+ TMath::Max(0.0, ele->pfIsolationVariables().sumNeutralHadronEt + ele->pfIsolationVariables().sumPhotonEt - 0.5*ele->pfIsolationVariables().sumPUPt ) 
		) /ele->pt() ; 
  if(iso>0.2) return false; 
  return true;



}


bool TriggerAnalyzerRAWMiniAOD::RecoHLTMatching(const edm::Event& iEvent, double recoeta, double recophi, std::string filtername, double dRmatching){
  //In the next few lines one loops over all the trigger objects (corresponding to a given filter) and check whether one of them matches the reco object under study                                       
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsORIGToken_, trigResults);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(trigobjectsMINIAODToken_, triggerObjects);

  const edm::TriggerNames &names = iEvent.triggerNames(*trigResults);
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackFilterLabels(iEvent,*trigResults);
    obj.unpackPathNames(names);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      std::string myfillabl=obj.filterLabels()[h];
      if( myfillabl.find(filtername)!=std::string::npos   && deltaR(recoeta,recophi, obj.eta(),obj.phi())<dRmatching ) return true;
    }
  }

  return false;
}




double TriggerAnalyzerRAWMiniAOD::VarStudied( const edm::Event& iEvent, double recoeta, double recophi,
					      edm::EDGetTokenT<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > varToken_,  edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> candToken_,   bool  dividebyE, bool dividebyEt, double dRmatching ){

  double thevar = 0.;

  //Inspired from http://cmslxr.fnal.gov/source/HLTrigger/Egamma/src/HLTGenericFilter.cc        
  edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;
  iEvent.getByToken (candToken_, PrevFilterOutput);

  edm::Handle<edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > > > depMap;
  iEvent.getByToken (varToken_,depMap);

  std::vector<edm::Ref<std::vector<reco::RecoEcalCandidate> > > recoCands;


  if(PrevFilterOutput.isValid()&&  depMap.isValid() ){

    PrevFilterOutput->getObjects(trigger::TriggerCluster, recoCands);
    if(recoCands.empty())PrevFilterOutput->getObjects(trigger::TriggerPhoton, recoCands);

    double dRmin = dRmatching;
    for (unsigned int i=0; i<recoCands.size(); i++) {
      edm::Ref<std::vector<reco::RecoEcalCandidate> > ref = recoCands[i];
      typename edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>, float > >::const_iterator mapi = (*depMap).find( ref );
      float vali = mapi->val;
      float EtaSC = ref->eta();
      float PhiSC = ref->phi();

      if(deltaR(recoeta,recophi,EtaSC,PhiSC ) > dRmin )continue;
      dRmin = deltaR(recoeta,recophi,EtaSC,PhiSC ) ;
      float energy = ref->superCluster()->energy();
      float et = ref->superCluster()->energy() * sin (2*atan(exp(-ref->eta())));
      thevar = (double) vali;
      if(dividebyE)thevar = (double)vali/energy;
      if(dividebyEt)thevar =(double) vali/et;

    }
  }
  return thevar;
}



//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnalyzerRAWMiniAOD);
