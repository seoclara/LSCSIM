#include "CupSim/CupPMTOpticalModel.hh"
#include "CupSim/CupPMTSD.hh"

#include "G4Version.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpticalSurface.hh"
#include "G4TransportationManager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"

#include "Randomize.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#include <complex.h>

#include "G4VProcess.hh" 

#include "G4GeometryTolerance.hh" // for kCarTolerance

G4UIdirectory *CupPMTOpticalModel::fgCmdDir = NULL;

// constructor -- also handles all initialization
// CupPMTOpticalModel::CupPMTOpticalModel (G4String modelName,
//				      G4VPhysicalVolume* envelope_phys)
//  : G4VFastSimulationModel(modelName, envelope_phys->GetLogicalVolume())
// EJ: start
CupPMTOpticalModel::CupPMTOpticalModel(G4String modelName, G4VPhysicalVolume *envelope_phys)
    : G4VFastSimulationModel(
          modelName,
          envelope_phys->GetLogicalVolume()->GetRegion()) // EJ: argument is changed from
                                                          // G4LogicalVolume to G4Region EJ: end
{
    //  _verbosity= 3;
    //  _luxlevel= 1;
    // get material properties vectors
    // ... material properties of glass
    G4LogicalVolume *envelope_log       = envelope_phys->GetLogicalVolume();
    G4MaterialPropertiesTable *glass_pt = envelope_log->GetMaterial()->GetMaterialPropertiesTable();
    if (glass_pt == NULL)
        // G4Exception("CupSim/CupPMTOpticalModel: glass lacks a properties table!");
        G4Exception(" ", " ", JustWarning,
                    "CupSim/CupPMTOpticalModel: glass lacks a properties table!");

    _rindex_glass = glass_pt->GetProperty("RINDEX");
    if (_rindex_glass == NULL)
        // G4Exception("CupSim/CupPMTOpticalModel: glass does not have RINDEX!");
        G4Exception(" ", " ", JustWarning,
                    "CupSim/CupPMTOpticalModel: glass does not have RINDEX!");

    // ... material properties of photocathode (first get photocathode surface)
    // here we assume that the first daughter volume is the "inner1" volume
    _inner1_phys = envelope_log->GetDaughter(0);
    //G4cout << "EJ: envelope= " << envelope_phys->GetName()
    //       << ", daughter= " << _inner1_phys->GetName() << G4endl;
    G4LogicalBorderSurface *pc_log_surface =
        G4LogicalBorderSurface::GetSurface(envelope_phys, _inner1_phys);
    if (pc_log_surface == NULL)
        // G4Exception("CupSim/CupPMTOpticalModel: no photocathode logical surface!?!");
        G4Exception(" ", " ", JustWarning,
                    "CupSim/CupPMTOpticalModel: no photocathode logical surface!?!");

#if (G4VERSION_NUMBER < 460)
    // it used to be so easy...
    G4OpticalSurface *pc_opsurf = pc_log_surface->GetOpticalSurface();
#else
    // G4LogicalSurface::GetOpticalSurface() function was eliminated from Geant4
    // version 6.  G4OpticalSurface now inherits from G4SurfaceProperty.
    // We have to trust that the G4SurfaceProperty returned by
    // pc_log_surface->GetSurfaceProperty() will be a G4OpticalSurface,
    // and cast the pointer accordingly by faith alone!
    G4OpticalSurface *pc_opsurf = (G4OpticalSurface *)(pc_log_surface->GetSurfaceProperty());
#endif

    if (pc_opsurf == nullptr) {
        // G4Exception("CupSim/CupPMTOpticalModel: no photocathode optical surface!?!");
        G4Exception(" ", " ", JustWarning,
                    "CupSim/CupPMTOpticalModel: no photocathode optical surface!?!");
    } else {
        G4MaterialPropertiesTable *pc_pt = pc_opsurf->GetMaterialPropertiesTable();
        if (pc_pt == nullptr) {
            // G4Exception("CupSim/CupPMTOpticalModel: photocathode lacks a properties table!");
            G4Exception(" ", " ", JustWarning,
                        "CupSim/CupPMTOpticalModel: photocathode lacks a properties table!");
        } else {
            _rindex_photocathode = pc_pt->GetProperty("RINDEX");
            if (_rindex_photocathode == nullptr)
                G4Exception(" ", " ", JustWarning,
                            "CupSim/CupPMTOpticalModel: photocathode does not have RINDEX!");
            _kindex_photocathode = pc_pt->GetProperty("KINDEX");
            if (_kindex_photocathode == nullptr)
                G4Exception(" ", " ", JustWarning,
                            "CupSim/CupPMTOpticalModel: photocathode does not have KINDEX!");

            _thickness_photocathode = pc_pt->GetProperty("THICKNESS");
            if (_thickness_photocathode == nullptr)
                G4Exception(" ", " ", JustWarning,
                            "CupSim/CupPMTOpticalModel: photocathode does not have THICKNESS!");

            _efficiency_photocathode = pc_pt->GetProperty("EFFICIENCY");
            if (_efficiency_photocathode == nullptr)
                G4Exception(" ", " ", JustWarning,
                            "CupSim/CupPMTOpticalModel: photocathode does not have EFFICIENCY!");
        }
    }
    _inner1_solid = _inner1_phys->GetLogicalVolume()->GetSolid();

    // initialize _photon_energy to a nonsense value to indicate that the other
    // values are not initialized
    _photon_energy = -1.0;

    // add UI commands
    if (fgCmdDir == NULL) {
        fgCmdDir = new G4UIdirectory("/PMTOpticalModel/");
        fgCmdDir->SetGuidance("PMT optical model control.");
        G4UIcommand *cmd;

        cmd = new G4UIcommand("/PMTOpticalModel/verbose", this);
        cmd->SetGuidance("Set verbose level\n"
                         " 0 == quiet\n"
                         " 1 == minimal entrance/exit info\n"
                         " 2 == +print verbose tracking info\n"
                         " >= 10  +lots of info on thin photocathode calcs\n");
        cmd->SetParameter(new G4UIparameter("level", 'i', false));

        cmd = new G4UIcommand("/PMTOpticalModel/luxlevel", this);
        cmd->SetGuidance("Set \"luxury level\" for PMT Optical Model\n"
                         " 0 == standard \"black bucket\": photons stop in PC, maybe make pe, \n"
                         " 1 == shiny translucent brown film: photons only stop if they make a PE, "
                         "otherwise 50/50 chance of reflecting/transmitting\n"
                         " 2 or greater == full model\n"
                         "The default value is 3.");
        cmd->SetParameter(new G4UIparameter("level", 'i', false));
    }
}

// destructor
CupPMTOpticalModel::~CupPMTOpticalModel() {
    // nothing to delete
    // Note: The "MaterialPropertyVector"s are owned by the material, not us.
}

// IsApplicable() method overriding virtual function of G4VFastSimulationModel
// returns true if model is applicable to given particle.
// -- see also Geant4 docs
G4bool CupPMTOpticalModel::IsApplicable(const G4ParticleDefinition &particleType) {
    return (&particleType == G4OpticalPhoton::OpticalPhotonDefinition());
}

// ModelTrigger() method overriding virtual function of G4VFastSimulationModel
// returns true if model should take over this specific track.
// -- see also Geant4 docs
G4bool CupPMTOpticalModel::ModelTrigger(const G4FastTrack &fastTrack) {
    // we trigger if the track position is above the equator
    // or if it is on the equator and heading up
    // if ( fastTrack.GetPrimaryTrackLocalPosition().z() >  kCarTolerance )
    if (fastTrack.GetPrimaryTrackLocalPosition().z() >
        G4GeometryTolerance::GetInstance()->GetSurfaceTolerance())
        return true;
    // if ( fastTrack.GetPrimaryTrackLocalPosition().z() > -kCarTolerance
    if (fastTrack.GetPrimaryTrackLocalPosition().z() >
            -(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) &&
        fastTrack.GetPrimaryTrackLocalDirection().z() > 0.0)
        return true;
    return false;
}

// DoIt() method overriding virtual function of G4VFastSimulationModel
// does the fast simulation for this track.  It is basically a faster but
// complete tracking code for the two-volume case.  It is a monster.
// -- see also Geant4 docs and comments below
void CupPMTOpticalModel::DoIt(const G4FastTrack &fastTrack, G4FastStep &fastStep) {
    // Logic summary:
    //  1) If track is outside the "inner1" vacuum, then track
    //     is advanced to either the inner1/body interface or to the body/outside
    //     interface, whichever comes first.
    //   a) If track hits the body/outside interface, then we just update
    //      its position and time and return;
    //   b) else, we do the thin-layer reflection/transmission/absorption thing
    //      with n1=n_glass and n3=1.0:
    //    i) Make a binary random decision on whether to absorb a photon.
    //       Probability is equal to "weight" of track times absorption coeff.
    //       (If P > 1.0, weight hit by floor(P+G4UniformRand()).)
    //       If photon is absorbed, then we are done tracking.
    //    ii) A binary random decision is made on whether to reflect or
    //       refract (transmit) the remaining track, if any.
    //  2) If track is in the "inner1" vacuum (either on entry to DoIt or
    //     after transmission in step 1.b.ii), then we advance to the
    //     inner1/body interface (surface of inner1).  Then we do the
    //     thin-layer reflection/transmission/absoprtion thing exactly as
    //     in step 1.b, but with n1=1.0 and n3=n_glass.
    //  3) Steps 1 and 2 are repeated until track hits the body outer surface,
    //     descends below the equator, or is completely absorbed on the
    //     photocathode.  Note it is best to check if step 1 applies even
    //     after doing step 1, to avoid making assumptions about the geometry.
    //     (It could happen someday that we have a photocathode with a negative
    //     curvature someplace.)
    //
    // Tracking elsewhere in the PMT (in the vacuum or glass below the equator,
    // for example) is handled by the usual Geant4 tracking.  If we ever
    // need to handle more than two volumes or more than one surface, then
    // we should use a general Geant4 "navigator" as in ExN05EMShowerModel,
    // but for this simple case, we can be more efficient with this custom
    // coding.  -GHS.

    G4double dist, dist1;
    G4ThreeVector pos;
    G4ThreeVector dir;
    G4ThreeVector pol;
    G4ThreeVector norm;
    G4double time;
    G4int weight;
    G4double energy;
    G4double n_glass;
    G4VSolid *envelope_solid = fastTrack.GetEnvelopeSolid();
    //  G4VSensitiveDetector* detector=
    //    fastTrack.GetEnvelope()->GetSensitiveDetector();
    // EJ: start
    G4VSensitiveDetector *detector =
        fastTrack.GetEnvelopeLogicalVolume()
            ->GetSensitiveDetector(); // EJ: argument is changed from G4LogicalVolume to G4Region
    // EJ: end
    enum EWhereAmI { kInGlass, kInVacuum } whereAmI;
    int ipmt = -1;

    // find which pmt we are in
    // The following doesn't work anymore (due to new geometry optimization?)
    //  ipmt=fastTrack.GetEnvelopePhysicalVolume()->GetMother()->GetCopyNo();
    // so we do this:
    {
        const G4VTouchable *touch = fastTrack.GetPrimaryTrack()->GetTouchable();
        int nd                    = touch->GetHistoryDepth();
        int id;
        for (id = 0; id < nd; id++) {
            if (touch->GetVolume(id) == fastTrack.GetEnvelopePhysicalVolume()) {
                ipmt = touch->GetReplicaNumber(id + 1);
                break;
            }
        }
        if (ipmt < 0) {
            // G4Exception("CupSim/CupPMTOpticalModel: could not find envelope -- where am I !?!");
            G4Exception(" ", " ", JustWarning,
                        "CupSim/CupPMTOpticalModel: could not find envelope -- where am I !?!");
        }
    }

    // processTag: EJ: 2007-11-06 (2016-0914)
    const G4VProcess *tProc = fastTrack.GetPrimaryTrack()->GetCreatorProcess();
    G4int processTag        = -1;
    if (tProc != NULL) {
        G4String procname = tProc->GetProcessName();
        if ((!strcmp(procname, "Cerenkov"))) processTag = 1;
        if ((!strcmp(procname, "Scintillation"))) processTag = 2;
        if ((!strcmp(procname, "Reemission"))) processTag = 3;
    }

    // get position and direction in local coordinates
    pos = fastTrack.GetPrimaryTrackLocalPosition();
    dir = fastTrack.GetPrimaryTrackLocalDirection();
    pol = fastTrack.GetPrimaryTrackLocalPolarization();

    // get weight and time
    time   = fastTrack.GetPrimaryTrack()->GetGlobalTime(); // "global" is correct
    weight = (G4int)(fastTrack.GetPrimaryTrack()->GetWeight());

    // get n_glass, _n2, _k2, etc., for this wavelength
    energy = fastTrack.GetPrimaryTrack()->GetKineticEnergy();
    if (energy == _photon_energy) // equal to last energy?
    {
        // use old values
        if (_n1 == 1.0)
            n_glass = _n3;
        else
            n_glass = _n1;
    } else {
        _photon_energy = energy;
        _wavelength    = twopi * hbarc / energy;
        // n_glass= _rindex_glass->GetProperty( energy );
        n_glass = _rindex_glass->Value(energy);
        _n1     = n_glass; // just in case we exit before setting _n1
        //_n2= _rindex_photocathode->GetProperty( energy );
        _n2 = _rindex_photocathode->Value(energy);
        //_k2= _kindex_photocathode->GetProperty( energy );
        _k2 = _kindex_photocathode->Value(energy);
        _n3 = 1.0; // just in case we exit before setting _n3
        //_efficiency= _efficiency_photocathode->GetProperty( energy );
        _efficiency = _efficiency_photocathode->Value(energy);
    }

    // initialize "whereAmI"
    if (fastTrack.GetPrimaryTrack()->GetVolume() == _inner1_phys)
        whereAmI = kInVacuum;
    else
        whereAmI = kInGlass;

    // print verbose info
    if (_verbosity > 0) {
        G4cout.flush();
        G4cout << "> Enter CupPMTOpticalModel, ipmt=" << ipmt
               << (whereAmI == kInVacuum ? " vacuum" : " glass") << ", pos=" << pos
               << ", dir=" << dir << ", weight=" << weight << ", pol=" << pol
               << ", energy=" << _photon_energy << ", wavelength=" << _wavelength
               << ", (n1,n3,n2,k2,efficiency)=(" << _n1 << "," << _n3 << "," << _n2 << "," << _k2
               << "," << _efficiency << ")\n";
    }

    G4int iloop;
    G4int max_iloop = 100;
    for (iloop = 0; iloop < max_iloop; iloop++) {

        if (whereAmI == kInGlass) { // in the glass
            // advance to next interface
            dist1 = envelope_solid->DistanceToOut(pos, dir);
            dist  = _inner1_solid->DistanceToIn(pos, dir);
            if (dist1 < dist) {
                // we hit the envelope outer surface, not the inner surface
                dist = dist1;
                if (dir.z() < 0.0)              // headed towards equator?
                {                               // make sure we don't cross the equator
                    dist1 = -pos.z() / dir.z(); // distance to equator
                    if (dist1 < dist) dist = dist1;
                }
                pos += dist * dir;
                time += dist * n_glass / c_light;
                break;
            }
            pos += dist * dir;
            time += dist * n_glass / c_light;
            _n1 = n_glass;
            _n3 = 1.0;
        } else { // in the "inner1" vacuum
            // advance to next interface
            dist = _inner1_solid->DistanceToOut(pos, dir);
            if (dist < 0.0) {
                G4cerr << "CupSim/CupPMTOpticalModel::DoIt(): "
                       << "Warning, strangeness detected! inner1->DistanceToOut()=" << dist
                       << G4endl;
                dist = 0.0;
            }
            pos += dist * dir;
            time += dist / c_light;
            // if ( pos.z() < kCarTolerance ) // we're passing through the equator
            if (pos.z() < G4GeometryTolerance::GetInstance()
                              ->GetSurfaceTolerance()) // we're passing through the equator
                break;
            _n1 = 1.0;
            _n3 = n_glass;
        }

        if (_verbosity >= 2) {
            G4cout << " " << iloop << " dist=" << dist << " newpos=" << pos << G4endl;
        }

        // get outward-pointing normal in local coordinates at this position
        norm = _inner1_solid->SurfaceNormal(pos);
        // reverse sign if incident from glass, so normal points into region "3"
        if (whereAmI == kInGlass)
            norm *= -1.0; // in principle, this is more efficient than norm= -norm;

        // set _thickness and _cos_theta1
        //_thickness=  _thickness_photocathode->GetProperty( pos.z() );
        _thickness  = _thickness_photocathode->Value(pos.z());
        _cos_theta1 = dir * norm;
        if (_cos_theta1 < 0.0) {
            G4cerr << "CupSim/CupPMTOpticalModel::DoIt(): "
                   << " The normal points the wrong way!\n"
                   << "  norm: " << norm << G4endl << "  dir:  " << dir << G4endl
                   << "  _cos_theta1:  " << _cos_theta1 << G4endl << "  pos:  " << pos << G4endl
                   << "  whereAmI:  " << (int)(whereAmI) << G4endl << " Reversing normal!"
                   << G4endl;
            _cos_theta1 = -_cos_theta1;
            norm        = -norm;
        }

        // Now calculate coefficients
        CalculateCoefficients();

        // Calculate Transmission, Reflection, and Absorption coefficients
        G4double T, R, A, An, collection_eff;
        G4double E_s2;
        if (_sin_theta1 > 0.0) {
            E_s2 = (pol * dir.cross(norm)) / _sin_theta1;
            E_s2 *= E_s2;
        } else
            E_s2 = 0.0;
        T              = fT_s * E_s2 + fT_p * (1.0 - E_s2);
        R              = fR_s * E_s2 + fR_p * (1.0 - E_s2);
        A              = 1.0 - (T + R);
        An             = 1.0 - (fT_n + fR_n); // The absorption at normal incidence
        collection_eff = _efficiency / An;    // net QE = _efficiency for normal inc.

#ifdef G4DEBUG
        if (A < 0.0 || A > 1.0 || collection_eff < 0.0 || collection_eff > 1.0) {
            G4cerr << "CupSim/CupPMTOpticalModel::DoIt(): Strange coefficients!\n";
            G4cout << "T, R, A, An, weight: " << T << " " << R << " " << A << " " << An << " "
                   << weight << G4endl;
            G4cout << "collection eff, std QE: " << collection_eff << " " << _efficiency << G4endl;
            G4cout << "=========================================================" << G4endl;
            A = collection_eff = 0.5; // safe values???
        }
#endif

        // Now decide how many pe we make.
        // When weight == 1, probability of a pe is A*collection_eff.
        // There is a certain correlation between "a pe is made" and
        // "the track is absorbed", which is implemented correctly below for
        // the weight == 1 case, and as good as can be done for weight>1 case.
        if (_verbosity > 0) {
            G4cout << "EJ: weight= " << weight << ", A= " << A
                   << ", collection_eff= " << collection_eff << ", qefficiency= " << _efficiency
                   << G4endl;
        }
/*
        //    G4double mean_N_pe= weight*A*collection_eff;
        G4double ranno_absorb = G4UniformRand();
        //    G4int N_pe= (G4int)( mean_N_pe + (1.0-ranno_absorb) );
        //    G4int N_pe= (G4int)( mean_N_pe + (1.0-ranno_absorb) );
        G4int N_pe;
        if (ranno_absorb < _efficiency) N_pe = 1;
*/
        G4double mean_N_pe= weight*A*collection_eff;
        G4double ranno_absorb= G4UniformRand();
        G4int N_pe= (G4int)( mean_N_pe + (1.0-ranno_absorb) );
        if (N_pe > 0) {
            if (detector != NULL && detector->isActive())
                ((CupPMTSD *)detector)
                    ->SimpleHit(ipmt, time, energy, pos, dir, pol, N_pe,
                                processTag); // EJ
            if (_verbosity >= 2) {
                G4cout << "CupSim/CupPMTOpticalModel made " << N_pe << " pe\n";
            }
        }

        // Now maybe absorb the track.
        // The probability is independent of weight, and the entire
        // track is either absorbed or not as a whole.
        // This is consistent with how the track is treated in other
        // processes (absorption, G4OpBoundary, etc.), and the statistics
        // for the final number of pe detected overall is made consistent
        // by the poissonian statistics of number of tracks implemented in
        // CupScint.  (Actually, for weights > 1, the correlation between
        // number absorbed and number transmitted/reflected and subsequently
        // absorbed cannot be made as big as it should be, so the width of
        // the distribution of the total number of pe per event will be slightly
        // too large.  Only the weight=1 case is "guaranteed" to get everything
        // correct, assuming there are no bugs in the code.)
        if (ranno_absorb < A) {
            weight = 0;
            if (_verbosity >= 2) G4cout << "CupSim/CupPMTOpticalModel absorbed track\n";
            break;
        }

        // reflect or refract the unabsorbed track
        if (G4UniformRand() < R / (R + T)) { // reflect
            Reflect(dir, pol, norm);
            if (_verbosity >= 2) G4cout << "CupSim/CupPMTOpticalModel reflects track\n";
        } else { // transmit
            Refract(dir, pol, norm);
            if (whereAmI == kInGlass)
                whereAmI = kInVacuum;
            else
                whereAmI = kInGlass;
            if (_verbosity >= 2)
                G4cout << "CupSim/CupPMTOpticalModel transmits track, now in "
                       << (whereAmI == kInVacuum ? " vacuum" : " glass") << G4endl;
        }
    }

    fastStep.SetPrimaryTrackFinalPosition(pos);
    fastStep.SetPrimaryTrackFinalTime(time);
    fastStep.SetPrimaryTrackFinalMomentum(dir);
    fastStep.SetPrimaryTrackFinalPolarization(pol);
    // fastStep.SetPrimaryTrackPathLength( trackLength ); // does anyone care?
    if (weight <= 0) {
        fastStep.ProposeTrackStatus(fStopAndKill);
        if (weight < 0) {
            G4cerr << "CupSim/CupPMTOpticalModel::DoIt(): Logic error, weight = " << weight
                   << G4endl;
        }
    } else {
        // in case multiphoton has been partly absorbed and partly reflected
        fastStep.SetPrimaryTrackFinalEventBiasingWeight(weight);
    }
    if (iloop >= max_iloop) {
        G4cerr << "CupSim/CupPMTOpticalModel::DoIt(): Too many loops, particle trapped!"
               << " Killing it." << G4endl;
        fastStep.ProposeTrackStatus(fStopAndKill);
    }

    if (_verbosity > 0) {
        G4cout.flush();
        G4cout << "> Exit CupPMTOpticalModel, ipmt=" << ipmt
               << (whereAmI == kInVacuum ? " vacuum" : " glass") << ", pos=" << pos
               << ", dir=" << dir << ", weight=" << weight << ", pol=" << pol << ", iloop=" << iloop
               << "\n";
    }

    return;
}

// CalculateCoefficients() method used by DoIt() above.
// *** THE PHYSICS, AT LAST!!! :-) ***
// Correct formalism implemented by Dario Motta (CEA-Saclay) 23 Feb 2005

void CupPMTOpticalModel::CalculateCoefficients()
// calculate and set fR_s, etc.
{
    if (_luxlevel <= 0) {
        // no reflection or transmission, just a black "light bucket"
        // 100% absorption, and QE will be renormalized later
        fR_s = fR_p = 0.0;
        fT_s = fT_p = 0.0;
        fR_n        = 0.0;
        fT_n        = 0.0;
        return;
    } else if (_luxlevel == 1) {
        // this is what was calculated before, when we had no good defaults
        // for cathode thickness and complex rindex
        // set normal incidence coefficients: 50/50 refl/transm if not absorb.
        fR_n = fT_n = 0.5 * (1.0 - _efficiency);
        // set sines and cosines
        _sin_theta1 = sqrt(1.0 - _cos_theta1 * _cos_theta1);
        _sin_theta3 = _n1 / _n3 * _sin_theta1;
        if (_sin_theta3 > 1.0) {
            // total non-transmission -- what to do?
            // total reflection or absorption
            _cos_theta3 = 0.0;
            fR_s = fR_p = 1.0 - _efficiency;
            fT_s = fT_p = 0.0;
            return;
        }
        _cos_theta3 = sqrt(1.0 - _sin_theta3 * _sin_theta3);
        fR_s = fR_p = fR_n;
        fT_s = fT_p = fT_n;
        return;
    }
    // else...

    // declare the prototypes of some useful functions
    G4complex carcsin(G4complex theta); // complex sin^-1
    G4complex gfunc(G4complex ni, G4complex nj, G4complex ti, G4complex tj);
    G4complex rfunc(G4complex ni, G4complex nj, G4complex ti, G4complex tj);
    G4complex trfunc(G4complex ni, G4complex nj, G4complex ti, G4complex tj, G4complex tk);

    // declare some useful constants
    G4complex _n2comp(_n2, -_k2); // complex photocathode refractive index
    G4complex eta = twopi * _n2comp * _thickness / _wavelength;
    G4complex zi(0., 1.); // imaginary unit

    // declare local variables

    G4complex theta1, theta2, theta3, delta; // geometric parameters
    G4complex r12, r23, t12, t21, t23;       // reflection- and transmission-related terms
    G4complex ampr, ampt;                    // relfection and transmission amplitudes

    // first set sines and cosines
    _sin_theta1 = sqrt(1.0 - _cos_theta1 * _cos_theta1);
    _sin_theta3 = _n1 / _n3 * _sin_theta1;
    if (_sin_theta3 > 1.0) {
        // total non-transmission -- what to do???
        // these variables only used to decide refracted track direction,
        // so doing the following should be okay:
        _sin_theta3 = 1.0;
    }
    _cos_theta3 = sqrt(1.0 - _sin_theta3 * _sin_theta3);

    // Determine all angles
    theta1 = asin(_sin_theta1);                      // incidence angle
    theta2 = carcsin((_n1 / _n2comp) * _sin_theta1); // complex angle in the photocathode
    theta3 = carcsin((_n2comp / _n3) * sin(theta2)); // angle of refraction into vacuum
    if (imag(theta3) < 0.) theta3 = conj(theta3);    // needed! (sign ambiguity arcsin)

    delta = eta * cos(theta2);

    // Calculation for the s-polarization

    r12 = rfunc(_n1, _n2comp, theta1, theta2);
    r23 = rfunc(_n2comp, _n3, theta2, theta3);
    t12 = trfunc(_n1, _n2comp, theta1, theta1, theta2);
    t21 = trfunc(_n2comp, _n1, theta2, theta2, theta1);
    t23 = trfunc(_n2comp, _n3, theta2, theta2, theta3);

    ampr =
        r12 + (t12 * t21 * r23 * exp(-2. * zi * delta)) / (1. + r12 * r23 * exp(-2. * zi * delta));
    ampt = (t12 * t23 * exp(-zi * delta)) / (1. + r12 * r23 * exp(-2. * zi * delta));

    // And finally...!
    fR_s = real(ampr * conj(ampr));
    fT_s = real(gfunc(_n3, _n1, theta3, theta1) * ampt * conj(ampt));

    // Calculation for the p-polarization

    r12 = rfunc(_n1, _n2comp, theta2, theta1);
    r23 = rfunc(_n2comp, _n3, theta3, theta2);
    t12 = trfunc(_n1, _n2comp, theta1, theta2, theta1);
    t21 = trfunc(_n2comp, _n1, theta2, theta1, theta2);
    t23 = trfunc(_n2comp, _n3, theta2, theta3, theta2);

    ampr =
        r12 + (t12 * t21 * r23 * exp(-2. * zi * delta)) / (1. + r12 * r23 * exp(-2. * zi * delta));
    ampt = (t12 * t23 * exp(-zi * delta)) / (1. + r12 * r23 * exp(-2. * zi * delta));

    // And finally...!
    fR_p = real(ampr * conj(ampr));
    fT_p = real(gfunc(_n3, _n1, theta3, theta1) * ampt * conj(ampt));

    // Now calculate the reference values at normal incidence (to scale QE)

    delta = eta;
    // Calculation for both polarization (the same at normal incidence)
    r12 = rfunc(_n1, _n2comp, 0., 0.);
    r23 = rfunc(_n2comp, _n3, 0., 0.);
    t12 = trfunc(_n1, _n2comp, 0., 0., 0.);
    t21 = trfunc(_n2comp, _n1, 0., 0., 0.);
    t23 = trfunc(_n2comp, _n3, 0., 0., 0.);

    ampr =
        r12 + (t12 * t21 * r23 * exp(-2. * zi * delta)) / (1. + r12 * r23 * exp(-2. * zi * delta));
    ampt = (t12 * t23 * exp(-zi * delta)) / (1. + r12 * r23 * exp(-2. * zi * delta));

    // And finally...!
    fR_n = real(ampr * conj(ampr));
    fT_n = real(gfunc(_n3, _n1, 0., 0.) * ampt * conj(ampt));

#ifdef G4DEBUG
    if (_verbosity >= 10) {
        G4cout << "=> lam, n1, n2: " << _wavelength / nm << " " << _n1 << " " << _n2comp << G4endl;
        G4cout << "=> Angles: " << real(theta1) / degree << " " << theta2 / degree << " "
               << theta3 / degree << G4endl;
        G4cout << "Rper, Rpar, Tper, Tpar: " << fR_s << " " << fR_p << " " << fT_s << " " << fT_p;
        G4cout << "\nRn, Tn : " << fR_n << " " << fT_n;
        G4cout << "\n-------------------------------------------------------" << G4endl;
    }
#endif
}

G4complex carcsin(G4complex theta) // complex sin^-1
{
    G4complex zi(0., 1.);
    G4complex value = (1. / zi) * (log(zi * theta + sqrt(1. - theta * theta)));
    return value;
}

G4complex gfunc(G4complex ni, G4complex nj, G4complex ti, G4complex tj) {
    G4complex value = (ni * cos(ti)) / (nj * cos(tj));
    return value;
}

G4complex rfunc(G4complex ni, G4complex nj, G4complex ti, G4complex tj) {
    G4complex value = (ni * cos(ti) - nj * cos(tj)) / (ni * cos(ti) + nj * cos(tj));
    return value;
}

G4complex trfunc(G4complex ni, G4complex nj, G4complex ti, G4complex tj, G4complex tk) {
    G4complex value = 2. * (ni * cos(ti)) / (ni * cos(tj) + nj * cos(tk));
    return value;
}

// Reflect() method, used by DoIt()
void CupPMTOpticalModel::Reflect(G4ThreeVector &dir, G4ThreeVector &pol, G4ThreeVector &norm) {
    dir -= 2. * (dir * norm) * norm;
    pol -= 2. * (pol * norm) * norm;
}

// Refract() method, used by DoIt()
void CupPMTOpticalModel::Refract(G4ThreeVector &dir, G4ThreeVector &pol, G4ThreeVector &norm) {
    dir = (_cos_theta3 - _cos_theta1 * _n1 / _n3) * norm + (_n1 / _n3) * dir;
    pol = (pol - (pol * dir) * dir).unit();
}

// user command handling functions
void CupPMTOpticalModel::SetNewValue(G4UIcommand *command, G4String newValues) {
    G4String commandName = command->GetCommandName();
    if (commandName == "verbose") {
        _verbosity = strtol((const char *)newValues, NULL, 0);
    } else if (commandName == "luxlevel") {
        _luxlevel = strtol((const char *)newValues, NULL, 0);
    } else {
        G4cerr << "No PMTOpticalModel command named " << commandName << G4endl;
    }
    return;
}

G4String CupPMTOpticalModel::GetCurrentValue(G4UIcommand *command) {
    G4String commandName = command->GetCommandName();
    if (commandName == "verbose") {
        char outbuff[64];
        sprintf(outbuff, "%d", _verbosity);
        return G4String(outbuff);
    } else if (commandName == "luxlevel") {
        char outbuff[64];
        sprintf(outbuff, "%d", _luxlevel);
        return G4String(outbuff);
    } else {
        return (commandName + " is not a valid PMTOpticalModel command");
    }
}
