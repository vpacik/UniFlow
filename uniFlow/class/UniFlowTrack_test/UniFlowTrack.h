class AliVTrack;
// class AliAODTrack;

class UniFlowTrack
{
  public:

    UniFlowTrack() : fPt{-1.0}, fPhi{0.0}, fEta{0.0}, fWeight{0.0}, fMass{0.0}, fSpecies{0} {};
    UniFlowTrack(Double_t pt, Double_t phi, Double_t eta, Double_t w, Double_t mass = 0.0, Int_t spec = 0) : fPt{pt}, fPhi{phi}, fEta{eta}, fWeight{w}, fMass{mass}, fSpecies{spec} {};
    UniFlowTrack(const AliVTrack* track) : UniFlowTrack{} { if(track) { fPt = track->Pt(); } };
    ~UniFlowTrack() {};

    // bool operator < (const UniFlowTrack& track) const { return (fPt < track.fPt); }

    Double_t fPt;
    Double_t fPhi;
    Double_t fEta;
    Double_t fWeight;
    Double_t fMass;
    Int_t fSpecies;

  private:
};
