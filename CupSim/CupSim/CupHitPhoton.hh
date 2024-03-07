
#ifndef __CupHitPhoton_hh__
#define __CupHitPhoton_hh__

#include <iostream>

#include "G4ThreeVector.hh"
using namespace CLHEP;

class CupHitPhoton {
  public:
    CupHitPhoton() {}

    void SetPMTID(int id) { fPMTID = id; }
    void SetTime(double t) { fTime = t; }
    void SetKineticEnergy(double KE);
    void SetWavelength(double wl);
    void SetPosition(double x, double y, double z);
    void SetMomentum(double x, double y, double z);
    void SetPolarization(double x, double y, double z);
    void SetCount(int count) { fCount = count; }
    void AddCount(int dcount) { fCount += dcount; }
    void SetProcessTag(int proctag) { fProctag = proctag; } // EJ: 2007-11-06

    int GetPMTID() const { return fPMTID; }
    double GetTime() const { return fTime; }
    double GetKineticEnergy() const;
    double GetWavelength() const;
    template <class T>
    inline void GetPosition(T &x, T &y, T &z) const;
    template <class T>
    inline void GetMomentum(T &x, T &y, T &z) const;
    template <class T>
    inline void GetPolarization(T &x, T &y, T &z) const;
    int GetCount() const { return fCount; }
    int GetProcessTag() const { return fProctag; } // EJ: 2007-11-06
    void Print(std::ostream &) const;

  private:
    double fTime;           /// time of hit
    int fPMTID;             /// ID number of PMT the HitPhoton hit
    float fKE;              /// kinetic energy
    float fPosition[3];     /// x,y,z components of position
    float fMomentum[3];     /// x,y,z components of momentum (normalized?)
    float fPolarization[3]; /// x,y,z components of polarization
    int fCount;             /// count of photons, often 1
    int fProctag; /// process tag: 1=Cerenkov, 2=Scintillation, 3=Reemission: EJ:2007-11-06
};

template <class T>
inline void CupHitPhoton::GetPosition(T &x, T &y, T &z) const {
    x = fPosition[0];
    y = fPosition[1];
    z = fPosition[2];
}

template <class T>
inline void CupHitPhoton::GetMomentum(T &x, T &y, T &z) const {
    x = fMomentum[0];
    y = fMomentum[1];
    z = fMomentum[2];
}

template <class T>
inline void CupHitPhoton::GetPolarization(T &x, T &y, T &z) const {
    x = fPolarization[0];
    y = fPolarization[1];
    z = fPolarization[2];
}

/** comparison function for sorting CupHitPhoton pointers
 */
inline bool Compare_HitPhotonPtr_TimeAscending(const CupHitPhoton *a, const CupHitPhoton *b) {
    return a->GetTime() < b->GetTime();
}

#endif // __CupHitPhoton_hh__
