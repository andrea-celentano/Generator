//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
 University of Liverpool & STFC Rutherford Appleton Lab - July 04, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Feb 22, 2011 - JD
 Implemented dummy versions of the new GFluxI::Clear, GFluxI::Index and
 GFluxI::GenerateWeighted methods needed for pre-generation of flux
 interaction probabilities in GMCJDriver.

 */
//____________________________________________________________________________
#include <cassert>
#include <algorithm>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TVector3.h>

#include "Framework/Conventions/Constants.h"
#include "Tools/Flux/GBDXpipesHistoFlux.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/PrintUtils.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie, flux, GBDXpipesHistoFlux, genie::flux::GBDXpipesHistoFlux)

using namespace genie;
using namespace genie::constants;
using namespace genie::flux;

//____________________________________________________________________________
GBDXpipesHistoFlux::GBDXpipesHistoFlux() :
		GFluxExposureI(genie::flux::kPOTs) {

	this->Initialize();
}
//___________________________________________________________________________
GBDXpipesHistoFlux::~GBDXpipesHistoFlux() {
	this->CleanUp();
}
//___________________________________________________________________________
double GBDXpipesHistoFlux::GetTotalExposure() const {
	// complete the GFluxExposureI interface
	return UsedEOTs();
}

//___________________________________________________________________________
double GBDXpipesHistoFlux::UsedEOTs(void) const {
// Compute current number of flux EOTs
	return fAccumEOTs;
}

//___________________________________________________________________________
long int GBDXpipesHistoFlux::NFluxNeutrinos(void) const {
	///< number of flux neutrinos looped so far
	return fNNeutrinos;
}
//___________________________________________________________________________
bool GBDXpipesHistoFlux::GenerateNext(void) {
	//-- Reset previously generated neutrino code / 4-p / 4-x
	this->ResetSelection();

	double R, ctheta, stheta, phi, yy;
	double x, y, z;
	//-- Generate an energy from the 'combined' spectrum histogram
	double Ev = fTotSpectrum->GetRandom();

	//check the range
	int trials = 0;
	while (1) {
		if ((Ev > fMinEv) && (Ev < fMaxEv)) break;
		if (trials > fMaxTrials){
			LOG ("Flux", pWARN) << "Warn GenerateNext() had to break in energy loop";
			break;
		}
		Ev = fTotSpectrum->GetRandom();
		trials++;
	}

	//-- Select a neutrino species from the flux fractions at the
	//   selected energy
	int inu = this->SelectNeutrino(Ev);
	fgPdgC = (*fPdgCList)[inu];

	//-- Extract R,2*pi*(1-z) from the projected 3D spectrum.
	//-- In the projected 3D spectrum: x-axis is the ANGLE, y-axis is the RADIUS
	TH3D *hSpectrum3D = fSpectrum3D[inu];
	int theBin = hSpectrum3D->GetZaxis()->FindBin(Ev);
	fSpectrum3D_proj[inu][theBin - 1]->GetRandom2(yy, R);
	trials = 0;

	while (1) {
		if ((R < fRt)||(fRt <= 0)) break;
		if (trials > fMaxTrials){
			LOG ("Flux", pWARN) << "Warn GenerateNext() had to break in angle/radius loop";
			break;
		}
		fSpectrum3D_proj[inu][theBin - 1]->GetRandom2(yy, R);
		trials++;
	}

	//ctheta is here yy=2*pi*(1-ctheta). Convert back to ctheta = 1-yy/2pi
	ctheta = 1 - yy / (2. * kPi);
	stheta = sqrt(1 - ctheta * ctheta);

	//compute random phi -> both the XY phi and the Z phi
	phi = this->GeneratePhi();	// rndm angle [0,2pi]


	//compute neutrino 4-momentum
	fgP4.SetXYZT(Ev * stheta * cos(phi), Ev * stheta * sin(phi), Ev * ctheta, Ev);


	//Compute neutrino position
	x = R * cos(phi);
	y = R * sin(phi);
	z = 0;


	//Add the "center"
	if (fBeamSpot) {
		x = x + fBeamSpot->X();
		y = y + fBeamSpot->Y();
		z = z + fBeamSpot->Z();
	}

	//Important: the units in this driver are cm, but GENIE expects the output in S.I. units
	double fLengthScale = units::centimeter/units::meter;

	x=x*fLengthScale;
	y=y*fLengthScale;
	z=z*fLengthScale;


	fgX4.SetXYZT(x, y, z, 0.);


	LOG("Flux", pINFO) << "Generated neutrino pdg-code: " << fgPdgC;
	LOG("Flux", pINFO) << "Generated neutrino p4: " << utils::print::P4AsShortString(&fgP4);
	LOG("Flux", pINFO) << "Generated neutrino x4: " << utils::print::X4AsString(&fgX4);

	// update the # EOTs and number of neutrinos
	fAccumEOTs += fEffEOTsPerNu;
	fNNeutrinos++;

	return true;
}
//___________________________________________________________________________
void GBDXpipesHistoFlux::Clear(Option_t * opt) {
// Dummy clear method needed to conform to GFluxI interface 
//
	LOG("Flux", pERROR) <<
	"No Clear(Option_t * opt) method implemented for opt: " << opt;
}
//___________________________________________________________________________
void GBDXpipesHistoFlux::GenerateWeighted(bool gen_weighted) {
// Dummy implementation needed to conform to GFluxI interface
//
	LOG("Flux", pERROR) <<
	"No GenerateWeighted(bool gen_weighted) method implemented for " << "gen_weighted: " << gen_weighted;
}
//___________________________________________________________________________
void GBDXpipesHistoFlux::Initialize(void) {
	LOG("Flux", pNOTICE) << "Initializing GBDXpipesHistoFlux driver";

	fMaxEv = -1;
	fMinEv = -1;
	fPdgCList = new PDGCodeList;
	fTotSpectrum = 0;
	fTotNorm = 0;
	fBeamSpot = 0;
	fRt = -1;

	fMaxTrials = 100000;
	fEffEOTsPerNu = 1;
	fAccumEOTs = 0;
	fNNeutrinos = 0;

	this->ResetSelection();
}
//___________________________________________________________________________
void GBDXpipesHistoFlux::ResetSelection(void) {
// initializing running neutrino pdg-code, 4-position, 4-momentum
	fgPdgC = 0;
	fgP4.SetPxPyPzE(0., 0., 0., 0.);
	fgX4.SetXYZT(0., 0., 0., 0.);
}
//___________________________________________________________________________
void GBDXpipesHistoFlux::CleanUp(void) {
	LOG("Flux", pNOTICE) << "Cleaning up...";

	if (fBeamSpot) delete fBeamSpot;
	if (fPdgCList) delete fPdgCList;
	if (fTotSpectrum) delete fTotSpectrum;

	unsigned int nspectra = fSpectrum.size();
	for (unsigned int i = 0; i < nspectra; i++) {
		TH1D * spectrum = fSpectrum[i];
		delete spectrum;
		spectrum = 0;

		TH3D * spectrum3D = fSpectrum3D[i];
		int nbins = spectrum3D->GetNbinsZ();
		for (int ibin = 1; ibin <= nbins; ibin++) {
			delete fSpectrum3D_proj[i][ibin - 1];
		}
		delete[] fSpectrum3D_proj[i];

		delete spectrum3D;
		spectrum3D = 0;

	}
}

//___________________________________________________________________________
void GBDXpipesHistoFlux::SetBeamSpot(const TVector3 & spot) {
	if (fBeamSpot) delete fBeamSpot;
	fBeamSpot = new TVector3(spot);
}

void GBDXpipesHistoFlux::SetBeamSpotZ(double Zspot) {
	if (fBeamSpot) delete fBeamSpot;
	fBeamSpot = new TVector3(0, 0, Zspot);
}

//___________________________________________________________________________
void GBDXpipesHistoFlux::SetTransverseRadius(double Rt) {
	LOG ("Flux", pNOTICE) << "Setting R[transverse] = " << Rt;
	fRt = Rt;

}

//___________________________________________________________________________
void GBDXpipesHistoFlux::SetEmin(double Emin) {

	//Spectrum not yet set, do checks when setting it.
	if (fTotSpectrum == 0) {
		fMinEv = Emin;
		return;
	} else {
		double axisMax = fTotSpectrum->GetXaxis()->GetXmax();
		double axisMin = fTotSpectrum->GetXaxis()->GetXmin();
		if (Emin > axisMax) {
			LOG ("Flux", pWARN) << " fTotSpectrum axis max is: " << axisMax << " and the minimum energy was requested larger than this: " << Emin << " IGNORE THIS REQUEST ";
			fMinEv = axisMin;
		} else {
			LOG ("Flux", pNOTICE) << "Setting E[min] = " << Emin;
			fMinEv = Emin;
		}

		//need to recompute the normalization
		fTotNorm = fTotSpectrum->Integral(fTotSpectrum->FindBin(fMinEv), fTotSpectrum->FindBin(fMaxEv));
		fEffEOTsPerNu = 1. / fTotNorm;
		LOG ("Flux", pNOTICE) << "SetEmin, fixing norm between: " << fMinEv << " and " << fMaxEv << " : " << fTotNorm;
		LOG ("Flux", pNOTICE) << "EOT per nu: " << fEffEOTsPerNu;
	}
}

//___________________________________________________________________________
void GBDXpipesHistoFlux::SetEmax(double Emax) {

	//Spectrum not yet set, do checks when setting it.
	if (fTotSpectrum == 0) {
		fMaxEv = Emax;
		return;
	} else {
		double axisMax = fTotSpectrum->GetXaxis()->GetXmax();
		double axisMin = fTotSpectrum->GetXaxis()->GetXmin();
		if (Emax < axisMin) {
			LOG ("Flux", pWARN) << " fTotSpectrum axis min is: " << axisMin << " and the maximum energy was requested smaller than this: " << Emax << " IGNORE THIS REQUEST ";
			fMaxEv = axisMax;
		} else {
			LOG ("Flux", pNOTICE) << "Setting E[max] = " << Emax;
			fMaxEv = Emax;
		}

		//need to recompute the normalization
		fTotNorm = fTotSpectrum->Integral(fTotSpectrum->FindBin(fMinEv), fTotSpectrum->FindBin(fMaxEv));
		fEffEOTsPerNu = 1. / fTotNorm;
		LOG ("Flux", pNOTICE) << "SetEmax, fixing norm between: " << fMinEv << " and " << fMaxEv << " : " << fTotNorm;
		LOG ("Flux", pNOTICE) << "EOT per nu: " << fEffEOTsPerNu;
	}
}

//___________________________________________________________________________
void GBDXpipesHistoFlux::SetEnergyRange(double Emin, double Emax) {
	this->SetEmin(Emin);
	this->SetEmax(Emax);
}

//___________________________________________________________________________
//3D spectrum X: R (in cm), Y: 2*pi*(1-cos(theta)), Z: Energy (GeV)
void GBDXpipesHistoFlux::AddEnergySpectrum(int nu_pdgc, TH1D * spectrum, TH3D * spectrum3D) {
	LOG("Flux", pNOTICE) << "Adding flux spectrum for pdg = " << nu_pdgc;

	/*VERY IMPORTANT. ROOT uses the inversion method to extract random numberss from histograms, by pre-caching the integral
	 * u = int_0^x n(y) dy, where n(y) is the PDF.
	 *
	 *
	 *
	 * However, in doing so it IGNORES the bin width. Therefore if the bin width is non-uniform, the result is wrong.
	 *
	 * u_i = SUM_{j=1}^{j=1} BIN_CONTENT(j)
	 *
	 * AND NOT:
	 *
	 * u_i = SUM_{j=1}^{j=1} BIN_CONTENT(j) * BIN_WIDTH(j)
	 *
	 * Trick: change the spectrum!
	 *
	 */
	fPdgCList->push_back(nu_pdgc);

	bool accepted = (count(fPdgCList->begin(), fPdgCList->end(), nu_pdgc) == 1);
	if (!accepted) {
		LOG ("Flux", pWARN) << "The pdg-code isn't recognized and the spectrum was ignored";
	} else {

		if (spectrum->GetXaxis()->GetXmin() != spectrum3D->GetZaxis()->GetXmin()) {
			LOG ("Flux", pFATAL) << "histo has min x-axis range: " << spectrum->GetXaxis()->GetXmin() << " while histo3d has min z-axis range: " << spectrum3D->GetZaxis()->GetXmin() << " FATAL ERROR";
			exit(1);
		}
		if (spectrum->GetXaxis()->GetXmax() != spectrum3D->GetZaxis()->GetXmax()) {
			LOG ("Flux", pFATAL) << "histo has max x-axis range: " << spectrum->GetXaxis()->GetXmin() << " while histo3d has max z-axis range: " << spectrum3D->GetZaxis()->GetXmax() << " FATAL ERROR";
			exit(1);
		}

		fSpectrum.push_back(spectrum);
		fSpectrum3D.push_back(spectrum3D);

		//here do the trick on the spectrum (the integral was cached properly before this)
		int nbins = spectrum->GetNbinsX();
		for (int ii = 1; ii <= nbins; ii++) {
			double data = spectrum->GetBinContent(ii);
			double width = spectrum->GetBinWidth(ii);
			spectrum->SetBinContent(ii, data * width);
		}

		//prepare the projections - one projections for each z-bin of the 3D histogram
		TH2D **hArrayProjs = new TH2D*[nbins];
		TH2D *hProj;

		nbins = spectrum3D->GetNbinsZ();
		for (int ii = 1; ii <= nbins; ii++) {
			spectrum3D->GetZaxis()->SetRange(ii, ii);
			hProj = (TH2D*) spectrum3D->Project3D("XY");
			hArrayProjs[ii - 1] = (TH2D*) hProj->Clone();
			hArrayProjs[ii - 1]->SetName(Form("%s_%i", spectrum3D->GetName(), ii));
		}

		fSpectrum3D_proj.push_back(hArrayProjs);
		this->AddAllFluxes(); // update combined flux

		//fix fRt - the RADIUS is in the X axis
		double axisMax = spectrum3D->GetXaxis()->GetXmax();
		if (fRt < 0) { //rmax not yet set
			fRt = axisMax;
		}
	}
}

//___________________________________________________________________________
void GBDXpipesHistoFlux::AddAllFluxes(void) {
	LOG("Flux", pNOTICE) << "Computing combined flux";

	if (fTotSpectrum) delete fTotSpectrum;

	vector<TH1D *>::const_iterator spectrum_iter;

	unsigned int inu = 0;
	fTotNorm = 0;
	for (spectrum_iter = fSpectrum.begin(); spectrum_iter != fSpectrum.end(); ++spectrum_iter) {
		TH1D * spectrum = *spectrum_iter;

		if (inu == 0) {
			fTotSpectrum = new TH1D(*spectrum);
		} else {
			fTotSpectrum->Add(spectrum);
		}

		inu++;
	}

	/*Range checks*/
	double axisMin = fTotSpectrum->GetXaxis()->GetXmin();
	double axisMax = fTotSpectrum->GetXaxis()->GetXmax();

	if (fMaxEv < 0) {	//MAX - not yet set
		fMaxEv = axisMax;
	} else {	// MAX was set. Check it is reasonable.
		if (fMaxEv < axisMin) {
			LOG ("Flux", pWARN) << " fTotSpectrum axis min is: " << axisMin << " and the maximum energy was requested smaller than this: " << fMaxEv << " - RESET Emax to max histo range " << axisMax;
			fMaxEv = axisMax;
		}
	}

	if (fMinEv < 0) {	//MIN - not yet set
		fMinEv = axisMin;
	} else {	// MIN was set. Check it is reasonable.
		if (fMinEv > axisMax) {
			LOG ("Flux", pWARN) << " fTotSpectrum axis max is: " << axisMax << " and the minimum energy was requested bigger than this: " << fMinEv << " - RESET Emin to min histo range " << axisMin;
			fMinEv = axisMin;
		}
	}

	fTotNorm = fTotSpectrum->Integral(fTotSpectrum->FindBin(fMinEv), fTotSpectrum->FindBin(fMaxEv)); //very important: do not use the width, since spectrum was fixed before.
	fEffEOTsPerNu = 1. / fTotNorm;
	if (fTotNorm <= 0) {
		LOG ("Flux", pERROR) << "fTotNorm is <0: " << fTotNorm;
		exit(1);
	} else {
		LOG("Flux",pINFO) << "fTotNorm is: " << fTotNorm;
		LOG ("Flux", pNOTICE) << "EOT per nu: " << fEffEOTsPerNu;
	}

}
//___________________________________________________________________________
int GBDXpipesHistoFlux::SelectNeutrino(double Ev) {
	const unsigned int n = fPdgCList->size();
	double fraction[n];

	vector<TH1D *>::const_iterator spectrum_iter;

	unsigned int inu = 0;
	for (spectrum_iter = fSpectrum.begin(); spectrum_iter != fSpectrum.end(); ++spectrum_iter) {
		TH1D * spectrum = *spectrum_iter;
		fraction[inu++] = spectrum->GetBinContent(spectrum->FindBin(Ev));
	}

	double sum = 0;
	for (inu = 0; inu < n; inu++) {
		sum += fraction[inu];
		fraction[inu] = sum;
		LOG("Flux", pDEBUG) << "SUM-FRACTION(0->" << inu << ") = " << sum;
	}

	RandomGen * rnd = RandomGen::Instance();
	double R = sum * rnd->RndFlux().Rndm();

	LOG("Flux", pDEBUG) << "R e [0,SUM] = " << R;

	for (inu = 0; inu < n; inu++) {
		if (R < fraction[inu]) return inu;
	}

	LOG("Flux", pERROR) << "Could not select a neutrino species";
	assert(false);

	return -1;
}
//___________________________________________________________________________
double GBDXpipesHistoFlux::GeneratePhi(void) const {
	RandomGen * rnd = RandomGen::Instance();
	double phi = 2. * kPi * rnd->RndFlux().Rndm(); // [0,2pi]
	return phi;
}
//___________________________________________________________________________
