#ifndef REREGISTLEP_H
#define REREGISTLEP_H

#include "TTBarCPV/TTBarCPVAnalysisRun1/interface/format.h"

void reRegistLep( LepInfoBranches& OldLep, LepInfoBranches& NewLep)
{
    int size=0;
    for( int i=0; i<OldLep.Size; i++){
        size++;
        NewLep.Index[i] = OldLep.Index[i];
        NewLep.isEcalDriven[i] = OldLep.isEcalDriven[i];
        NewLep.isTrackerDriven[i] = OldLep.isTrackerDriven[i];
        NewLep.LeptonType[i] = OldLep.LeptonType[i];
        NewLep.Charge[i] = OldLep.Charge[i];
        NewLep.ChargeGsf[i] = OldLep.ChargeGsf[i];
        NewLep.ChargeCtf[i] = OldLep.ChargeCtf[i];
        NewLep.ChargeScPix[i] = OldLep.ChargeScPix[i];
        NewLep.Pt[i] = OldLep.Pt[i];
        NewLep.Et[i] = OldLep.Et[i];
        NewLep.Eta[i] = OldLep.Eta[i];
        NewLep.caloEta[i] = OldLep.caloEta[i];
        NewLep.Phi[i] = OldLep.Phi[i];
        NewLep.TrackIso[i] = OldLep.TrackIso[i];
        NewLep.EcalIso[i] = OldLep.EcalIso[i];
        NewLep.HcalIso[i] = OldLep.HcalIso[i];
        NewLep.HcalDepth1Iso[i] = OldLep.HcalDepth1Iso[i];
        NewLep.HcalDepth2Iso[i] = OldLep.HcalDepth2Iso[i];
        NewLep.ChargedHadronIso[i] = OldLep.ChargedHadronIso[i];
        NewLep.NeutralHadronIso[i] = OldLep.NeutralHadronIso[i];
        NewLep.PhotonIso[i] = OldLep.PhotonIso[i];

        NewLep.ChargedHadronIsoR03[i] = OldLep.ChargedHadronIsoR03[i];
        NewLep.NeutralHadronIsoR03[i] = OldLep.NeutralHadronIsoR03[i];
        NewLep.PhotonIsoR03[i] = OldLep.PhotonIsoR03[i];
        NewLep.sumPUPtR03[i] = OldLep.sumPUPtR03[i];
        NewLep.IsoRhoCorrR03[i] = OldLep.IsoRhoCorrR03[i];
        NewLep.ChargedHadronIsoR04[i] = OldLep.ChargedHadronIsoR04[i];
        NewLep.NeutralHadronIsoR04[i] = OldLep.NeutralHadronIsoR04[i];
        NewLep.PhotonIsoR04[i] = OldLep.PhotonIsoR04[i];
        NewLep.sumPUPtR04[i] = OldLep.sumPUPtR04[i];
        NewLep.IsoRhoCorrR04[i] = OldLep.IsoRhoCorrR04[i];

        NewLep.Ip3dPV[i] = OldLep.Ip3dPV[i];
        NewLep.Ip3dPVErr[i] = OldLep.Ip3dPVErr[i];
        NewLep.Ip3dPVSignificance[i] = OldLep.Ip3dPVSignificance[i];

        NewLep.MuontimenDof[i] = OldLep.MuontimenDof[i];
        NewLep.MuontimeAtIpInOut[i] = OldLep.MuontimeAtIpInOut[i];
        NewLep.MuontimeAtIpOutIn[i] = OldLep.MuontimeAtIpOutIn[i];
        NewLep.Muondirection[i] = OldLep.Muondirection[i];

        NewLep.CaloEnergy[i] = OldLep.CaloEnergy[i];
        NewLep.e1x5[i] = OldLep.e1x5[i];
        NewLep.e2x5Max[i] = OldLep.e2x5Max[i];
        NewLep.e5x5[i] = OldLep.e5x5[i];

        NewLep.Px[i] = OldLep.Px[i];
        NewLep.Py[i] = OldLep.Py[i];
        NewLep.Pz[i] = OldLep.Pz[i];
        NewLep.Energy[i] = OldLep.Energy[i];
        NewLep.sGoodMuonTMOneStationTight[i] = OldLep.sGoodMuonTMOneStationTight[i];
        NewLep.innerTracknormalizedChi2[i] = OldLep.innerTracknormalizedChi2[i];

        NewLep.vertexZ[i] = OldLep.vertexZ[i];

        NewLep.isPFMuon[i] = OldLep.isPFMuon[i];
        NewLep.MuIDGlobalMuonPromptTight[i] = OldLep.MuIDGlobalMuonPromptTight[i];

        NewLep.MuInnerPtError[i] = OldLep.MuInnerPtError[i];
        NewLep.MuGlobalPtError[i] = OldLep.MuGlobalPtError[i];
        NewLep.MuInnerTrackDz[i] = OldLep.MuInnerTrackDz[i];
        NewLep.MuInnerTrackD0[i] = OldLep.MuInnerTrackD0[i];
        NewLep.MuInnerTrackDxy_BS[i] = OldLep.MuInnerTrackDxy_BS[i];
        NewLep.MuInnerTrackDxy_PV[i] = OldLep.MuInnerTrackDxy_PV[i];
        NewLep.MuInnerTrackDxy_PVBS[i] = OldLep.MuInnerTrackDxy_PVBS[i];
        NewLep.MuInnerTrackNHits[i] = OldLep.MuInnerTrackNHits[i];
        NewLep.MuNTrackerHits[i] = OldLep.MuNTrackerHits[i];

        NewLep.MuGlobalNormalizedChi2[i] = OldLep.MuGlobalNormalizedChi2[i];

        NewLep.MuCaloCompat[i] = OldLep.MuCaloCompat[i];
        NewLep.MuNChambers[i] = OldLep.MuNChambers[i];
        NewLep.MuNChambersMatchesSegment[i] = OldLep.MuNChambersMatchesSegment[i];
        NewLep.MuNMatchedStations[i] = OldLep.MuNMatchedStations[i];
        NewLep.MuNPixelLayers[i] = OldLep.MuNPixelLayers[i];
        NewLep.MuNPixelLayersWMeasurement[i] = OldLep.MuNPixelLayersWMeasurement[i];
        NewLep.MuNTrackLayersWMeasurement[i] = OldLep.MuNTrackLayersWMeasurement[i];
        NewLep.MuNLostInnerHits[i] = OldLep.MuNLostInnerHits[i];
        NewLep.MuNLostOuterHits[i] = OldLep.MuNLostOuterHits[i];
        NewLep.MuNMuonhits[i] = OldLep.MuNMuonhits[i];
        NewLep.MuDThits[i] = OldLep.MuDThits[i];
        NewLep.MuCSChits[i] = OldLep.MuCSChits[i];
        NewLep.MuRPChits[i] = OldLep.MuRPChits[i];
        NewLep.MuType[i] = OldLep.MuType[i];

        NewLep.EgammaMVANonTrig[i] = OldLep.EgammaMVANonTrig[i];
        NewLep.EgammaMVATrig[i] = OldLep.EgammaMVATrig[i];
        NewLep.EgammaCutBasedEleIdTRIGGERTIGHT[i] = OldLep.EgammaCutBasedEleIdTRIGGERTIGHT[i];
        NewLep.EgammaCutBasedEleIdTRIGGERWP70[i] = OldLep.EgammaCutBasedEleIdTRIGGERWP70[i];
        NewLep.EgammaCutBasedEleIdVETO[i] = OldLep.EgammaCutBasedEleIdVETO[i];
        NewLep.EgammaCutBasedEleIdLOOSE[i] = OldLep.EgammaCutBasedEleIdLOOSE[i];
        NewLep.EgammaCutBasedEleIdMEDIUM[i] = OldLep.EgammaCutBasedEleIdMEDIUM[i];
        NewLep.EgammaCutBasedEleIdTIGHT[i] = OldLep.EgammaCutBasedEleIdTIGHT[i];

        NewLep.Eldr03HcalDepth1TowerSumEtBc[i] = OldLep.Eldr03HcalDepth1TowerSumEtBc[i];
        NewLep.Eldr03HcalDepth2TowerSumEtBc[i] = OldLep.Eldr03HcalDepth2TowerSumEtBc[i];
        NewLep.Eldr04HcalDepth1TowerSumEtBc[i] = OldLep.Eldr04HcalDepth1TowerSumEtBc[i];
        NewLep.Eldr04HcalDepth2TowerSumEtBc[i] = OldLep.Eldr04HcalDepth2TowerSumEtBc[i];
        NewLep.ElhcalOverEcalBc[i] = OldLep.ElhcalOverEcalBc[i];
        NewLep.ElEcalE[i] = OldLep.ElEcalE[i];
        NewLep.ElEoverP[i] = OldLep.ElEoverP[i];
        NewLep.EldeltaEta[i] = OldLep.EldeltaEta[i];
        NewLep.EldeltaPhi[i] = OldLep.EldeltaPhi[i];
        NewLep.ElHadoverEm[i] = OldLep.ElHadoverEm[i];
        NewLep.ElsigmaIetaIeta[i] = OldLep.ElsigmaIetaIeta[i];
        NewLep.ElscSigmaIetaIeta[i] = OldLep.ElscSigmaIetaIeta[i];
        NewLep.ElEnergyErr[i] = OldLep.ElEnergyErr[i];
        NewLep.ElMomentumErr[i] = OldLep.ElMomentumErr[i];
        NewLep.ElTrackNHits[i] = OldLep.ElTrackNHits[i];
        NewLep.ElSharedHitsFraction[i] = OldLep.ElSharedHitsFraction[i];
        NewLep.dR_gsf_ctfTrack[i] = OldLep.dR_gsf_ctfTrack[i];
        NewLep.dPt_gsf_ctfTrack[i] = OldLep.dPt_gsf_ctfTrack[i];
        NewLep.ElhasConv[i] = OldLep.ElhasConv[i];

        NewLep.ElTrackNLostHits[i] = OldLep.ElTrackNLostHits[i];
        NewLep.ElTrackDz[i] = OldLep.ElTrackDz[i];
        NewLep.ElTrackDz_BS[i] = OldLep.ElTrackDz_BS[i];
        NewLep.ElTrackD0[i] = OldLep.ElTrackD0[i];
        NewLep.ElTrackDxy_BS[i] = OldLep.ElTrackDxy_BS[i];
        NewLep.ElTrackDxy_PV[i] = OldLep.ElTrackDxy_PV[i];
        NewLep.ElTrackDxy_PVBS[i] = OldLep.ElTrackDxy_PVBS[i];
        NewLep.ElNClusters[i] = OldLep.ElNClusters[i];
        NewLep.ElClassification[i] = OldLep.ElClassification[i];
        NewLep.ElFBrem[i] = OldLep.ElFBrem[i];
        NewLep.NumberOfExpectedInnerHits[i] = OldLep.NumberOfExpectedInnerHits[i];
        NewLep.Eldist[i] = OldLep.Eldist[i];
        NewLep.Eldcot[i] = OldLep.Eldcot[i];
        NewLep.Elconvradius[i] = OldLep.Elconvradius[i];
        NewLep.ElConvPoint_x[i] = OldLep.ElConvPoint_x[i];
        NewLep.ElConvPoint_y[i] = OldLep.ElConvPoint_y[i];
        NewLep.ElConvPoint_z[i] = OldLep.ElConvPoint_z[i];

        NewLep.dcotdist[i] = OldLep.dcotdist[i];
        NewLep.ElseedEoverP[i] = OldLep.ElseedEoverP[i];
        NewLep.ElEcalIso04[i] = OldLep.ElEcalIso04[i];
        NewLep.ElHcalIso04[i] = OldLep.ElHcalIso04[i];

        NewLep.ElNumberOfBrems[i] = OldLep.ElNumberOfBrems[i];
        NewLep.GenPt[i] = OldLep.GenPt[i];
        NewLep.GenEta[i] = OldLep.GenEta[i];
        NewLep.GenPhi[i] = OldLep.GenPhi[i];
        NewLep.GenPdgID[i] = OldLep.GenPdgID[i];
        NewLep.GenMCTag[i] = OldLep.GenMCTag[i];
        NewLep.TrgPt[i] = OldLep.TrgPt[i];
        NewLep.TrgEta[i] = OldLep.TrgEta[i];
        NewLep.TrgPhi[i] = OldLep.TrgPhi[i];
        NewLep.TrgID[i] = OldLep.TrgID[i];

        //tau : hpsPFTau ID  
        NewLep.isPFTau[i] = OldLep.isPFTau[i];
        NewLep.decayModeFinding[i] = OldLep.decayModeFinding[i];
        NewLep.byVLooseCombinedIsolationDeltaBetaCorr[i] = OldLep.byVLooseCombinedIsolationDeltaBetaCorr[i];
        NewLep.byLooseCombinedIsolationDeltaBetaCorr[i] = OldLep.byLooseCombinedIsolationDeltaBetaCorr[i];
        NewLep.byMediumCombinedIsolationDeltaBetaCorr[i] = OldLep.byMediumCombinedIsolationDeltaBetaCorr[i];
        NewLep.byTightCombinedIsolationDeltaBetaCorr[i] = OldLep.byTightCombinedIsolationDeltaBetaCorr[i];
        NewLep.againstElectronLoose[i] = OldLep.againstElectronLoose[i];
        NewLep.againstElectronMedium[i] = OldLep.againstElectronMedium[i];
        NewLep.againstElectronTight[i] = OldLep.againstElectronTight[i];
        NewLep.againstElectronMVA[i] = OldLep.againstElectronMVA[i];
        NewLep.againstMuonLoose[i] = OldLep.againstMuonLoose[i];
        NewLep.againstMuonMedium[i] = OldLep.againstMuonMedium[i];
        NewLep.againstMuonTight[i] = OldLep.againstMuonTight[i];
    }
}

#endif
