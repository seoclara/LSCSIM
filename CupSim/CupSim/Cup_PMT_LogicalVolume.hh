
#ifndef __Cup_PMT_LogicalVolume_hh__
#define __Cup_PMT_LogicalVolume_hh__

#include "Cup_PMT_LogicalVolume.hh" // EJ: for Cup_PMT_LogicalVolume.hh
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"

class G4Material;
class G4OpticalSurface;

using namespace CLHEP;

typedef enum { kPmtStyle_TorusStack, kPmtStyle_Ellipsoid } ePmtStyle;

// Cup_PMT_LogicalVolume is parent class for specific tube classes below
class Cup_PMT_LogicalVolume : public G4LogicalVolume {
  public:
    Cup_PMT_LogicalVolume(const G4String &plabel, // label -- subvolume names are derived from this
                          G4double r_bound,       // radius of bounding cylinder
                          G4double hh_bound,      // half height of bounding cylinder
                          G4Material *ExteriorMat); // material which fills the bounding cylinder

    G4double GetZEquator() { return z_equator; }

  protected:
    G4double z_equator; // Z location of equator of tube

    static G4OpticalSurface *our_Mirror_opsurf;

    void ConstructPMT_UsingTorusStack(
        const G4int n_edge, const G4double outer_z_edge[], const G4double outer_rho_edge[],
        const G4double outer_z_o[],
        G4double r_dynode,                     // radius of dynode stack
        G4double z_dynode,                     // z coordinate of top of dynode stack, equator=0.
        G4double d_wall,                       // thickness of the walls
        G4Material *Exterior,                  // material outside tube
        G4Material *Glass,                     // glass material
        G4OpticalSurface *Photocathode_opsurf, // photocathode surface
        G4Material *PMT_Vac,                   // tube interior
        G4Material *Dynode_mat,                // dynode stack metal
        G4VSensitiveDetector *detector         // sensitive detector hook
    );

    void ConstructPMT_UsingEllipsoid(
        G4double r_equat,                      // equatorial radius of phototube
        G4double h_face,                       // height of face above eq.
        G4double h_back,                       // depth of back below eq.
        G4double r_stem,                       // radius of stem
        G4double l_stem,                       // length of stem
        G4double r_dynode,                     // radius of dynode stack
        G4double z_dynode,                     // z coordinate of top of dynode stack, equator=0.
        G4double d_wall,                       // thickness of the walls
        G4Material *Exterior,                  // material outside tube
        G4Material *Glass,                     // glass material
        G4OpticalSurface *Photocathode_opsurf, // photocathode surface
        G4Material *PMT_Vac,                   // tube interior
        G4Material *Dynode_mat,                // dynode stack metal
        G4VSensitiveDetector *detector         // sensitive detector hook
    );
};

// Hamamatsu R7250 PMT (so-called "17-inch") with optional face mask:
// note that the R7250 has the same outside shape and size as the
// R3600 "20-inch" PMT, but has less usable photocathode area.
class Cup_17inch_LogicalVolume : public Cup_PMT_LogicalVolume {
  public:
    Cup_17inch_LogicalVolume(
        const G4String &plabel,                // label -- subvolume names are derived from this
        G4Material *ExteriorMat,               // material which fills the bounding cylinder
        G4Material *GlassMat,                  // glass material
        G4OpticalSurface *Photocathode_opsurf, // photocathode surface
        G4Material *PMT_Vacuum,                // vacuum inside tube
        G4Material *DynodeMat,                 // dynode material
        G4Material *MaskMat,                   // material for photocathode mask (e.g, blk acryl)
                                               // OK to set MaskMat == NULL for no mask
        G4VSensitiveDetector *detector,        // sensitive detector hook
        ePmtStyle PmtStyle                     // PMT approx. style (TorusStack or Ellipsoid)
    );
};

// Hamamatsu R3600 ("20-inch") PMT with optional face mask
class Cup_20inch_LogicalVolume : public Cup_PMT_LogicalVolume {
  public:
    Cup_20inch_LogicalVolume(
        const G4String &plabel,                // label -- subvolume names are derived from this
        G4Material *ExteriorMat,               // material which fills the bounding cylinder
        G4Material *GlassMat,                  // glass material
        G4OpticalSurface *Photocathode_opsurf, // photocathode surface
        G4Material *PMT_Vacuum,                // vacuum inside tube
        G4Material *DynodeMat,                 // dynode material
        G4Material *MaskMat,                   // material for photocathode mask (e.g, blk acryl)
                                               // OK to set MaskMat == NULL for no mask
        G4VSensitiveDetector *detector,        // sensitive detector hook
        ePmtStyle PmtStyle                     // PMT approx. style (TorusStack or Ellipsoid)
    );
};

class Cup_10inch_LogicalVolume : public Cup_PMT_LogicalVolume {
  public:
    Cup_10inch_LogicalVolume(
        const G4String &plabel,                // label -- subvolume names are derived from this
        G4Material *ExteriorMat,               // material which fills the bounding cylinder
        G4Material *GlassMat,                  // glass material
        G4OpticalSurface *Photocathode_opsurf, // photocathode surface
        G4Material *PMT_Vacuum,                // vacuum inside tube
        G4Material *DynodeMat,                 // dynode material
        G4Material *MaskMat,                   // material for photocathode mask (e.g, blk acryl)
        // OK to set MaskMat == NULL for no mask
        G4VSensitiveDetector *detector, // sensitive detector hook
        ePmtStyle PmtStyle              // PMT approx. style (TorusStack or Ellipsoid)
    );
};

// Hamamatsu R1408 ("8-inch") PMT with optional face mask
//  (note R5912 and R4558 use same glass envelope, but dynode may be different)
class Cup_8inch_LogicalVolume : public Cup_PMT_LogicalVolume {
  public:
    Cup_8inch_LogicalVolume(
        const G4String &plabel,                // label -- subvolume names are derived from this
        G4Material *ExteriorMat,               // material which fills the bounding cylinder
        G4Material *GlassMat,                  // glass material
        G4OpticalSurface *Photocathode_opsurf, // photocathode surface
        G4Material *PMT_Vacuum,                // vacuum inside tube
        G4Material *DynodeMat,                 // dynode material
        G4Material *MaskMat,                   // material for photocathode mask (e.g, blk acryl)
                                               // OK to set MaskMat == NULL for no mask
        G4VSensitiveDetector *detector,        // sensitive detector hook
        ePmtStyle PmtStyle                     // PMT approx. style (TorusStack or Ellipsoid)
    );
};

class Cup_5inch_LogicalVolume : public Cup_PMT_LogicalVolume {
  public:
    Cup_5inch_LogicalVolume(
        const G4String &plabel,                // label -- subvolume names are derived from this
        G4Material *ExteriorMat,               // material which fills the bounding cylinder
        G4Material *GlassMat,                  // glass material
        G4OpticalSurface *Photocathode_opsurf, // photocathode surface
        G4Material *PMT_Vacuum,                // vacuum inside tube
        G4Material *DynodeMat,                 // dynode material
        G4Material *MaskMat,                   // material for photocathode mask (e.g, blk acryl)
                                               // OK to set MaskMat == NULL for no mask
        G4VSensitiveDetector *detector,        // sensitive detector hook
        ePmtStyle PmtStyle                     // PMT approx. style (TorusStack or Ellipsoid)
    );
};

// Hamamatsu H7195 ("2-inch") PMT with optional face mask
// EJ added 2inch PMT
class Cup_2inch_LogicalVolume : public Cup_PMT_LogicalVolume {
  public:
    Cup_2inch_LogicalVolume(
        const G4String &plabel,                // label -- subvolume names are derived from this
        G4Material *ExteriorMat,               // material which fills the bounding cylinder
        G4Material *GlassMat,                  // glass material
        G4OpticalSurface *Photocathode_opsurf, // photocathode surface
        G4Material *PMT_Vacuum,                // vacuum inside tube
        G4Material *DynodeMat,                 // dynode material
        G4Material *MaskMat,                   // material for photocathode mask (e.g, blk acryl)
                                               // OK to set MaskMat == NULL for no mask
        G4VSensitiveDetector *detector,        // sensitive detector hook
        ePmtStyle PmtStyle                     // PMT approx. style (TorusStack or Ellipsoid)
    );
};

#endif /*__Cup_PMT_LogicalVolume_hh__*/
