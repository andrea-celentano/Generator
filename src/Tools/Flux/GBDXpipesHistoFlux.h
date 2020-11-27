//____________________________________________________________________________
/*!

 \class   genie::flux::GBDXpipesHistoFlux

 \brief   A generic GENIE flux driver.
 Generates a 'cylindrical' neutrino beam along the input direction,
 with the input transverse radius and centered at the input position.
 The energies are generated from the input energy spectrum (TH1D).
 Multiple neutrino species can be generated (you will need to supply
 an energy spectrum for each).

 \author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
 University of Liverpool & STFC Rutherford Appleton Lab

 \created July 4, 2005

 \cpright  Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE
 */
//____________________________________________________________________________
#ifndef _G_BDX_HISTO_FLUX_H
#define _G_BDX_HISTO_FLUX_H

#include <string>
#include <vector>

#include <TLorentzVector.h>

#include "Framework/EventGen/GFluxI.h"
#include "Tools/Flux/GFluxExposureI.h"
#include "Tools/Flux/GFluxFileConfigI.h"
#include "Framework/ParticleData/PDGUtils.h"

class TH1D;
class TH2D;
class TH3D;
class TVector3;

using std::string;
using std::vector;

namespace genie {
namespace flux {

class GBDXpipesHistoFlux: public GFluxI, public genie::flux::GFluxExposureI {
//, public genie::flux::GFluxFileConfigI {

public:
	GBDXpipesHistoFlux();
	~GBDXpipesHistoFlux();

	// methods specific to this flux object
	void SetBeamSpot(const TVector3 & spot);
	void SetBeamSpotZ(double Zspot);
	void SetTransverseRadius(double Rt);
	void SetEnergyRange(double Emin, double Emax);
	void SetEmin(double Emin);
	void SetEmax(double Emax);
	void AddEnergySpectrum(int nu_pdgc, TH1D * spectrum1D, TH3D * spectrum3D);

	// methods implementing the GENIE GFluxI interface
	const PDGCodeList & FluxParticles(void) {
		return *fPdgCList;
	}
	double MaxEnergy(void) {
		return fMaxEv;
	}
	double MinEnergy(void) {
		return fMinEv;
	}
	bool GenerateNext(void);
	int PdgCode(void) {
		return fgPdgC;
	}
	double Weight(void) {
		return 1.0;
	}
	const TLorentzVector & Momentum(void) {
		return fgP4;
	}
	const TLorentzVector & Position(void) {
		return fgX4;
	}
	bool End(void) {
		return false;
	}
	long int Index(void) {
		return -1;
	}
	void Clear(Option_t * opt);
	void GenerateWeighted(bool gen_weighted);

	virtual double GetTotalExposure() const;  ///< GFluxExposureI interface
	virtual long int NFluxNeutrinos() const;    ///< # of rays generated

	double UsedEOTs(void) const;       ///< # of protons-on-target used

private:

	// private methods
	void Initialize(void);
	void CleanUp(void);
	void ResetSelection(void);
	void AddAllFluxes(void);
	int SelectNeutrino(double Ev);
	double GeneratePhi(void) const;

	// private data members
	double fMaxEv;       ///< maximum energy
	double fMinEv;      ///<  minimum energy
	PDGCodeList * fPdgCList;    ///< list of neutrino pdg-codes
	int fgPdgC;       ///< running generated nu pdg-code
	TLorentzVector fgP4;         ///< running generated nu 4-momentum
	TLorentzVector fgX4;         ///< running generated nu 4-position
	vector<TH1D *> fSpectrum;    ///< flux = f(Ev), 1/neutrino species
	vector<TH3D *> fSpectrum3D;  ///distrubition of f(R,2pi(1-z),E)

	vector<TH2D **> fSpectrum3D_proj;


	TH1D * fTotSpectrum; ///< combined flux = f(Ev)
	double fTotNorm;

	TVector3 * fBeamSpot;    ///< beam spot position
	double fRt;          ///< transverse size of neutrino beam

	int fMaxTrials;

	//exposure
	double    fEffEOTsPerNu;        ///< what a entry is worth ...
	double    fAccumEOTs;           ///< POTs used so far
	long int  fNNeutrinos;          ///< number of flux neutrinos thrown so far

};

} // flux namespace
} // genie namespace

#endif // _G_TH1_CYLICDRICAL_FLUX_H_
