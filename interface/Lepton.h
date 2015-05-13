#ifndef LEPTON_H 
#define LEPTON_H 

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h"

class Lepton{
	public:
		Lepton(){}
		Lepton( LepInfoBranches& LepInfo, int idx ){
			Fill( LepInfo, idx );
		}
		void Fill( LepInfoBranches& LepInfo, int idx ){
			Index = LepInfo.Index[idx];
			isEcalDriven = LepInfo.isEcalDriven[idx];
			isTrackerDriven = LepInfo.isTrackerDriven[idx];
			LeptonType = LepInfo.LeptonType[idx];
			Charge = LepInfo.Charge[idx];
			ChargeGsf = LepInfo.ChargeGsf[idx];
			ChargeCtf = LepInfo.ChargeCtf[idx];
			ChargeScPix = LepInfo.ChargeScPix[idx];
			Pt = LepInfo.Pt[idx];
			Et = LepInfo.Et[idx];
			Eta = LepInfo.Eta[idx];
			caloEta = LepInfo.caloEta[idx];
			Phi = LepInfo.Phi[idx];
			TrackIso = LepInfo.TrackIso[idx];
			EcalIso = LepInfo.EcalIso[idx];
			HcalIso = LepInfo.HcalIso[idx];
			HcalDepth1Iso = LepInfo.HcalDepth1Iso[idx];
			HcalDepth2Iso = LepInfo.HcalDepth2Iso[idx];
			ChargedHadronIso = LepInfo.ChargedHadronIso[idx];
			NeutralHadronIso = LepInfo.NeutralHadronIso[idx];
			PhotonIso = LepInfo.PhotonIso[idx];

			ChargedHadronIsoR03 = LepInfo.ChargedHadronIsoR03[idx];
			NeutralHadronIsoR03 = LepInfo.NeutralHadronIsoR03[idx];
			PhotonIsoR03 = LepInfo.PhotonIsoR03[idx];
			sumPUPtR03 = LepInfo.sumPUPtR03[idx];
			IsoRhoCorrR03 = LepInfo.IsoRhoCorrR03[idx];
			ChargedHadronIsoR04 = LepInfo.ChargedHadronIsoR04[idx];
			NeutralHadronIsoR04 = LepInfo.NeutralHadronIsoR04[idx];
			PhotonIsoR04 = LepInfo.PhotonIsoR04[idx];
			sumPUPtR04 = LepInfo.sumPUPtR04[idx];
			IsoRhoCorrR04 = LepInfo.IsoRhoCorrR04[idx];

			Ip3dPV = LepInfo.Ip3dPV[idx];
			Ip3dPVErr = LepInfo.Ip3dPVErr[idx];
			Ip3dPVSignificance = LepInfo.Ip3dPVSignificance[idx];

			MuontimenDof = LepInfo.MuontimenDof[idx];
			MuontimeAtIpInOut = LepInfo.MuontimeAtIpInOut[idx];
			MuontimeAtIpOutIn = LepInfo.MuontimeAtIpOutIn[idx];

			Muondirection = LepInfo.Muondirection[idx];

			CaloEnergy = LepInfo.CaloEnergy[idx];
			e1x5 = LepInfo.e1x5[idx];
			e2x5Max = LepInfo.e2x5Max[idx];
			e5x5 = LepInfo.e5x5[idx];

			Px = LepInfo.Px[idx]; 
			Py = LepInfo.Py[idx]; 
			Pz = LepInfo.Pz[idx]; 
			Energy = LepInfo.Energy[idx]; 
			isGoodMuonTMOneStationTight = LepInfo.isGoodMuonTMOneStationTight[idx]; 
			innerTracknormalizedChi2 = LepInfo.innerTracknormalizedChi2[idx];   

			vertexZ = LepInfo.vertexZ[idx]; 

			isPFMuon = LepInfo.isPFMuon[idx];
			MuIDGlobalMuonPromptTight = LepInfo.MuIDGlobalMuonPromptTight[idx];

			MuInnerPtError = LepInfo.MuInnerPtError[idx];  
			MuGlobalPtError = LepInfo.MuGlobalPtError[idx];  
			MuInnerTrackDz = LepInfo.MuInnerTrackDz[idx];  
			MuInnerTrackD0 = LepInfo.MuInnerTrackD0[idx];  
			MuInnerTrackDxy_BS = LepInfo.MuInnerTrackDxy_BS[idx];  
			MuInnerTrackDxy_PV = LepInfo.MuInnerTrackDxy_PV[idx];  
			MuInnerTrackDxy_PVBS = LepInfo.MuInnerTrackDxy_PVBS[idx];  
			MuInnerTrackNHits = LepInfo.MuInnerTrackNHits[idx];
			MuNTrackerHits = LepInfo.MuNTrackerHits[idx];

			MuGlobalNormalizedChi2 = LepInfo.MuGlobalNormalizedChi2[idx]; 

			MuCaloCompat = LepInfo.MuCaloCompat[idx];
			MuNChambers = LepInfo.MuNChambers[idx];
			MuNChambersMatchesSegment = LepInfo.MuNChambersMatchesSegment[idx];
			MuNMatchedStations = LepInfo.MuNMatchedStations[idx];
			MuNPixelLayers = LepInfo.MuNPixelLayers[idx];
			MuNPixelLayersWMeasurement = LepInfo.MuNPixelLayersWMeasurement[idx]; 
			MuNTrackLayersWMeasurement = LepInfo.MuNTrackLayersWMeasurement[idx];
			MuNLostInnerHits = LepInfo.MuNLostInnerHits[idx];
			MuNLostOuterHits = LepInfo.MuNLostOuterHits[idx];
			MuNMuonhits = LepInfo.MuNMuonhits[idx];
			MuDThits = LepInfo.MuDThits[idx];
			MuCSChits = LepInfo.MuCSChits[idx];
			MuRPChits = LepInfo.MuRPChits[idx];
			MuType = LepInfo.MuType[idx]; 
			// GlobalMuon 0x02 -> 0000:0010-> 2
			// TrackerMuon 0x04 -> 0000:0100 -> 4
			// StandAloneMuon 0x08 -> 0000:1000 -> 8
			// CaloMuon 0x16 -> 0001:0110 -> 16+6 = 22

			EgammaMVANonTrig = LepInfo.EgammaMVANonTrig[idx]; 
			EgammaMVATrig = LepInfo.EgammaMVATrig[idx]; 
			EgammaCutBasedEleIdTRIGGERTIGHT = LepInfo.EgammaCutBasedEleIdTRIGGERTIGHT[idx]; 
			EgammaCutBasedEleIdTRIGGERWP70 = LepInfo.EgammaCutBasedEleIdTRIGGERWP70[idx]; 
			EgammaCutBasedEleIdVETO = LepInfo.EgammaCutBasedEleIdVETO[idx]; 
			EgammaCutBasedEleIdLOOSE = LepInfo.EgammaCutBasedEleIdLOOSE[idx]; 
			EgammaCutBasedEleIdMEDIUM = LepInfo.EgammaCutBasedEleIdMEDIUM[idx]; 
			EgammaCutBasedEleIdTIGHT = LepInfo.EgammaCutBasedEleIdTIGHT[idx]; 

			Eldr03HcalDepth1TowerSumEtBc = LepInfo.Eldr03HcalDepth1TowerSumEtBc[idx];
			Eldr03HcalDepth2TowerSumEtBc = LepInfo.Eldr03HcalDepth2TowerSumEtBc[idx];
			Eldr04HcalDepth1TowerSumEtBc = LepInfo.Eldr04HcalDepth1TowerSumEtBc[idx];
			Eldr04HcalDepth2TowerSumEtBc = LepInfo.Eldr04HcalDepth2TowerSumEtBc[idx];
			ElhcalOverEcalBc = LepInfo.ElhcalOverEcalBc[idx];
			ElEcalE = LepInfo.ElEcalE[idx];
			ElEoverP = LepInfo.ElEoverP[idx];
			EldeltaEta = LepInfo.EldeltaEta[idx];
			EldeltaPhi = LepInfo.EldeltaPhi[idx]; 
			ElHadoverEm = LepInfo.ElHadoverEm[idx];
			ElsigmaIetaIeta = LepInfo.ElsigmaIetaIeta[idx];	
			ElscSigmaIetaIeta = LepInfo.ElscSigmaIetaIeta[idx];	
			ElEnergyErr = LepInfo.ElEnergyErr[idx];
			ElMomentumErr = LepInfo.ElMomentumErr[idx];
			ElTrackNHits = LepInfo.ElTrackNHits[idx]; 
			ElSharedHitsFraction = LepInfo.ElSharedHitsFraction[idx]; 
			dR_gsf_ctfTrack = LepInfo.dR_gsf_ctfTrack[idx]; 
			dPt_gsf_ctfTrack = LepInfo.dPt_gsf_ctfTrack[idx]; 
			ElhasConv = LepInfo.ElhasConv[idx];

			ElTrackNLostHits = LepInfo.ElTrackNLostHits[idx];  
			ElTrackDz = LepInfo.ElTrackDz[idx];  
			ElTrackDz_BS = LepInfo.ElTrackDz_BS[idx];  
			ElTrackD0 = LepInfo.ElTrackD0[idx];  
			ElTrackDxy_BS = LepInfo.ElTrackDxy_BS[idx];  
			ElTrackDxy_PV = LepInfo.ElTrackDxy_PV[idx];  
			ElTrackDxy_PVBS = LepInfo.ElTrackDxy_PVBS[idx]; 
			ElNClusters = LepInfo.ElNClusters[idx];
			ElClassification = LepInfo.ElClassification[idx];
			ElFBrem = LepInfo.ElFBrem[idx];
			NumberOfExpectedInnerHits = LepInfo.NumberOfExpectedInnerHits[idx]; 
			Eldist = LepInfo.Eldist[idx]; 
			Eldcot = LepInfo.Eldcot[idx]; 
			Elconvradius = LepInfo.Elconvradius[idx]; 
			ElConvPoint_x = LepInfo.ElConvPoint_x[idx]; 
			ElConvPoint_y = LepInfo.ElConvPoint_y[idx]; 
			ElConvPoint_z = LepInfo.ElConvPoint_z[idx]; 

			dcotdist = LepInfo.dcotdist[idx];
			ElseedEoverP = LepInfo.ElseedEoverP[idx];
			ElEcalIso04 = LepInfo.ElEcalIso04[idx];
			ElHcalIso04 = LepInfo.ElHcalIso04[idx];

			ElNumberOfBrems = LepInfo.ElNumberOfBrems[idx];
			GenPt = LepInfo.GenPt[idx];
			GenEta = LepInfo.GenEta[idx];
			GenPhi = LepInfo.GenPhi[idx];
			GenPdgID = LepInfo.GenPdgID[idx];
			GenMCTag = LepInfo.GenMCTag[idx]; 	

			TrgPt = LepInfo.TrgPt[idx];
			TrgEta = LepInfo.TrgEta[idx];
			TrgPhi = LepInfo.TrgPhi[idx];
			TrgID = LepInfo.TrgID[idx];

			isPFTau = LepInfo.isPFTau[idx];   
			decayModeFinding = LepInfo.decayModeFinding[idx];
			byVLooseCombinedIsolationDeltaBetaCorr = LepInfo.byVLooseCombinedIsolationDeltaBetaCorr[idx];
			byLooseCombinedIsolationDeltaBetaCorr = LepInfo.byLooseCombinedIsolationDeltaBetaCorr[idx]; 
			byMediumCombinedIsolationDeltaBetaCorr = LepInfo.byMediumCombinedIsolationDeltaBetaCorr[idx]; 
			byTightCombinedIsolationDeltaBetaCorr = LepInfo.byTightCombinedIsolationDeltaBetaCorr[idx]; 
			againstElectronLoose = LepInfo.againstElectronLoose[idx]; 
			againstElectronMedium = LepInfo.againstElectronMedium[idx]; 
			againstElectronTight = LepInfo.againstElectronTight[idx]; 
			againstElectronMVA = LepInfo.againstElectronMVA[idx]; 
			againstMuonLoose = LepInfo.againstMuonLoose[idx]; 
			againstMuonMedium = LepInfo.againstMuonMedium[idx]; 
			againstMuonTight = LepInfo.againstMuonTight[idx];

			P3.SetXYZ( Px, Py, Pz );
			P4.SetPxPyPzE( Px, Py, Pz, Energy );
		}

		TVector3 P3;
		TLorentzVector P4;	
		
		int	Index;
		int	isEcalDriven;
		int	isTrackerDriven;
		int	LeptonType;
		int	Charge;
		int	ChargeGsf;
		int	ChargeCtf;
		int	ChargeScPix;
		float Pt;
		float Et;
		float Eta;
		float caloEta;
		float Phi;
		float TrackIso;
		float EcalIso;
		float HcalIso;
		float HcalDepth1Iso;
		float HcalDepth2Iso;
		float ChargedHadronIso;
		float NeutralHadronIso;
		float PhotonIso;

		float ChargedHadronIsoR03;
		float NeutralHadronIsoR03;
		float PhotonIsoR03;
		float sumPUPtR03;
		float IsoRhoCorrR03;
		float ChargedHadronIsoR04;
		float NeutralHadronIsoR04;
		float PhotonIsoR04;
		float sumPUPtR04;
		float IsoRhoCorrR04;

		float Ip3dPV;
		float Ip3dPVErr;
		float Ip3dPVSignificance;

		int   MuontimenDof;
		float MuontimeAtIpInOut;
		float MuontimeAtIpOutIn;
		// enum Direction { OutsideIn = -1, Undefined = 0, InsideOut = 1 };
		int   Muondirection;

		float CaloEnergy;
		float e1x5;
		float e2x5Max;
		float e5x5;

		float Px; 
		float Py; 
		float Pz; 
		float Energy; 
		bool isGoodMuonTMOneStationTight; 
		float innerTracknormalizedChi2;   

		float vertexZ; 

		bool  isPFMuon;
		bool  MuIDGlobalMuonPromptTight;

		float MuInnerPtError;  
		float MuGlobalPtError;  
		float MuInnerTrackDz;  
		float MuInnerTrackD0;  
		float MuInnerTrackDxy_BS;  
		float MuInnerTrackDxy_PV;  
		float MuInnerTrackDxy_PVBS;  
		int   MuInnerTrackNHits;
		int   MuNTrackerHits;

		float MuGlobalNormalizedChi2; 

		float MuCaloCompat;
		int   MuNChambers;
		int   MuNChambersMatchesSegment;
		int   MuNMatchedStations;
		int   MuNPixelLayers;
		int   MuNPixelLayersWMeasurement; 
		int   MuNTrackLayersWMeasurement;
		int   MuNLostInnerHits;
		int   MuNLostOuterHits;
		int   MuNMuonhits;
		int   MuDThits;
		int   MuCSChits;
		int   MuRPChits;
		int   MuType;

		float EgammaMVANonTrig; 
		float EgammaMVATrig; 
		bool EgammaCutBasedEleIdTRIGGERTIGHT; 
		bool EgammaCutBasedEleIdTRIGGERWP70; 
		bool EgammaCutBasedEleIdVETO; 
		bool EgammaCutBasedEleIdLOOSE; 
		bool EgammaCutBasedEleIdMEDIUM; 
		bool EgammaCutBasedEleIdTIGHT; 

		float Eldr03HcalDepth1TowerSumEtBc;
		float Eldr03HcalDepth2TowerSumEtBc;
		float Eldr04HcalDepth1TowerSumEtBc;
		float Eldr04HcalDepth2TowerSumEtBc;
		float ElhcalOverEcalBc;
		float ElEcalE;
		float ElEoverP;
		float EldeltaEta;
		float EldeltaPhi; 
		float ElHadoverEm;
		float ElsigmaIetaIeta;	
		float ElscSigmaIetaIeta;	
		float ElEnergyErr;
		float ElMomentumErr;
		int	ElTrackNHits; 
		float ElSharedHitsFraction; 
		float dR_gsf_ctfTrack; 
		float dPt_gsf_ctfTrack; 
		bool	ElhasConv;

		float ElTrackNLostHits;  
		float ElTrackDz;  
		float ElTrackDz_BS;  
		float ElTrackD0;  
		float ElTrackDxy_BS;  
		float ElTrackDxy_PV;  
		float ElTrackDxy_PVBS; 
		int	ElNClusters;
		int	ElClassification;
		float	ElFBrem;
		int NumberOfExpectedInnerHits; 
		float Eldist; 
		float Eldcot; 
		float Elconvradius; 
		float ElConvPoint_x; 
		float ElConvPoint_y; 
		float ElConvPoint_z; 

		// For CIC	
		float dcotdist;
		float ElseedEoverP;
		float ElEcalIso04;
		float ElHcalIso04;

		int	ElNumberOfBrems;
		float GenPt;
		float GenEta;
		float GenPhi;
		int	GenPdgID;
		int	GenMCTag; 	// 0: unknown, 1: decay from W, 2: decay from Z, 
		// 3: from b, 4: from c, 5: match to a parton (q or g), 6: match to a photon
		// (+10) from b'
		// (+20) from t'
		float TrgPt;
		float TrgEta;
		float TrgPhi;
		int TrgID;

		//tau : hpsPFTau ID  
		int isPFTau;   
		float decayModeFinding;
		float byVLooseCombinedIsolationDeltaBetaCorr;
		float byLooseCombinedIsolationDeltaBetaCorr; 
		float byMediumCombinedIsolationDeltaBetaCorr; 
		float byTightCombinedIsolationDeltaBetaCorr; 
		float againstElectronLoose; 
		float againstElectronMedium; 
		float againstElectronTight; 
		float againstElectronMVA; 
		float againstMuonLoose; 
		float againstMuonMedium; 
		float againstMuonTight; 
};

#endif
