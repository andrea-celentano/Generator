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
#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TVector3.h>

#include "Framework/Conventions/Constants.h"
#include "Tools/Flux/GBDXpipesHistoFlux_DM.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/PrintUtils.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie, flux, GBDXpipesHistoFlux_DM, genie::flux::GBDXpipesHistoFlux_DM)

using namespace genie;
using namespace genie::constants;
using namespace genie::flux;

//____________________________________________________________________________
GBDXpipesHistoFlux_DM::GBDXpipesHistoFlux_DM() :
		GFluxExposureI(genie::flux::kPOTs) {

	this->Initialize();
}
//___________________________________________________________________________
GBDXpipesHistoFlux_DM::~GBDXpipesHistoFlux_DM() {
	this->CleanUp();
}
//___________________________________________________________________________
double GBDXpipesHistoFlux_DM::GetTotalExposure() const {
	// complete the GFluxExposureI interface
	return UsedEOTs();
}

//___________________________________________________________________________
double GBDXpipesHistoFlux_DM::UsedEOTs(void) const {
// Compute current number of flux EOTs
	return fAccumEOTs;
}

//___________________________________________________________________________
long int GBDXpipesHistoFlux_DM::NFluxNeutrinos(void) const {
	///< number of flux neutrinos looped so far
	return fNDM;
}
//___________________________________________________________________________
bool GBDXpipesHistoFlux_DM::GenerateNext(void) {
	//-- Reset previously generated neutrino code / 4-p / 4-x
	this->ResetSelection();

	double ctheta, stheta, phi, yy;
	double x, y, z;
	//-- Generate an energy from the 'combined' spectrum histogram
	double E = fTotSpectrum->GetRandom();

	//check the range
	int trials = 0;
	while (1) {
		if ((E > fMinE) && (E < fMaxE)) break;
		if (trials > fMaxTrials) {
			LOG ("Flux", pWARN) << "Warn GenerateNext() had to break in energy loop";
			break;
		}
		E = fTotSpectrum->GetRandom();
		trials++;
	}

	//-- Select a DM species from the flux fractions at the
	//   selected energy
	int iDM = this->SelectDM(E);
	fgPdgC = (*fPdgCList)[iDM];

	//-- Extract z=sin(theta) from the projected 2D spectrum.
	//-- In the projected 2D spectrum: x-axis is the angle
	TH2D *hSpectrum2D = fSpectrum2D[iDM];
	int theBin = hSpectrum2D->GetYaxis()->FindBin(E);
	yy = fSpectrum2D_proj[iDM][theBin - 1]->GetRandom();


	//stheta here is yy
	stheta = yy;
	ctheta = sqrt(1 - stheta * stheta);

	//compute random phi -> both the XY phi and the Z phi
	phi = this->GeneratePhi();	// rndm angle [0,2pi]




	//compute DM 4-momentum
	fgP4.SetXYZT(E * stheta * cos(phi), E * stheta * sin(phi), E * ctheta, E);

	//Compute DM position --> assume prompt position at (0,0,0)
	x = 0;
	y = 0;
	z = 0;

	//Add the "center"
	if (fBeamSpot) {
		x = x + fBeamSpot->X();
		y = y + fBeamSpot->Y();
		z = z + fBeamSpot->Z();
	}

	//Important: the units in this driver are cm, but GENIE expects the output in S.I. units
	double fLengthScale = units::centimeter / units::meter;

	x = x * fLengthScale;
	y = y * fLengthScale;
	z = z * fLengthScale;

	fgX4.SetXYZT(x, y, z, 0.);

	LOG("Flux", pINFO) << "Generated DM code: " << fgPdgC;
	LOG("Flux", pINFO) << "Generated DM p4: " << utils::print::P4AsShortString(&fgP4);
	LOG("Flux", pINFO) << "Generated DM x4: " << utils::print::X4AsString(&fgX4);

	// update the # EOTs and number of neutrinos
	fAccumEOTs += fEffEOTsPerDM;
	fNDM++;

	std::cout<<"Generated neutrino DM-code: " << fgPdgC<<std::endl;
	std::cout<<"Generated neutrino p4: " << utils::print::P4AsShortString(&fgP4)<<std::endl;
	std::cout<<"Generated neutrino x4: " << utils::print::P4AsShortString(&fgX4)<<std::endl;
	return true;
}
//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::Clear(Option_t *opt) {
// Dummy clear method needed to conform to GFluxI interface 
//


	LOG("Flux", pERROR) <<
	"Clear(Option_t * opt) called for opt: " << opt;



	 fNDM = 0;
	 fAccumEOTs  = 0;
}
//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::GenerateWeighted(bool gen_weighted) {
// Dummy implementation needed to conform to GFluxI interface
//
	LOG("Flux", pERROR) <<
	"No GenerateWeighted(bool gen_weighted) method implemented for " << "gen_weighted: " << gen_weighted;
}
//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::Initialize(void) {
	LOG("Flux", pNOTICE) << "Initializing GBDXpipesHistoFlux_DM driver";

	fMaxE = -1;
	fMinE = -1;
	fPdgCList = new PDGCodeList;
	fTotSpectrum = 0;
	fTotNorm = 0;
	fBeamSpot = 0;

	fMaxTrials = 100000;
	fEffEOTsPerDM = 1;
	fAccumEOTs = 0;
	fNDM = 0;

	this->ResetSelection();
}
//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::ResetSelection(void) {
// initializing running neutrino pdg-code, 4-position, 4-momentum
	fgPdgC = 0;
	fgP4.SetPxPyPzE(0., 0., 0., 0.);
	fgX4.SetXYZT(0., 0., 0., 0.);
}
//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::CleanUp(void) {
	LOG("Flux", pNOTICE) << "Cleaning up...";

	if (fBeamSpot) delete fBeamSpot;
	if (fPdgCList) delete fPdgCList;
	if (fTotSpectrum) delete fTotSpectrum;

	unsigned int nspectra = fSpectrum.size();
	for (unsigned int i = 0; i < nspectra; i++) {
		TH1D *spectrum = fSpectrum[i];
		delete spectrum;
		spectrum = 0;

		TH2D *spectrum2D = fSpectrum2D[i];
		int nbins = spectrum2D->GetNbinsY();
		for (int ibin = 1; ibin <= nbins; ibin++) {
			delete fSpectrum2D_proj[i][ibin - 1];
		}
		delete[] fSpectrum2D_proj[i];

		delete spectrum2D;
		spectrum2D = 0;

	}
}

//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::SetBeamSpot(const TVector3 &spot) {
	if (fBeamSpot) delete fBeamSpot;
	fBeamSpot = new TVector3(spot);
}

void GBDXpipesHistoFlux_DM::SetBeamSpotZ(double Zspot) {
	if (fBeamSpot) delete fBeamSpot;
	fBeamSpot = new TVector3(0, 0, Zspot);
}

//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::SetEmin(double Emin) {

	//Spectrum not yet set, do checks when setting it.
	if (fTotSpectrum == 0) {
		fMinE = Emin;
		return;
	} else {
		double axisMax = fTotSpectrum->GetXaxis()->GetXmax();
		double axisMin = fTotSpectrum->GetXaxis()->GetXmin();
		if (Emin > axisMax) {
			LOG ("Flux", pWARN) << " fTotSpectrum axis max is: " << axisMax << " and the minimum energy was requested larger than this: " << Emin << " IGNORE THIS REQUEST ";
			fMinE = axisMin;
		} else {
			LOG ("Flux", pNOTICE) << "Setting E[min] = " << Emin;
			fMinE = Emin;
		}

		//need to recompute the normalization
		fTotNorm = fTotSpectrum->Integral(fTotSpectrum->FindBin(fMinE), fTotSpectrum->FindBin(fMaxE));
		fEffEOTsPerDM = 1. / fTotNorm;
		LOG ("Flux", pNOTICE) << "SetEmin, fixing norm between: " << fMinE << " and " << fMaxE << " : " << fTotNorm;
		LOG ("Flux", pNOTICE) << "EOT per DM: " << fEffEOTsPerDM;
	}
}

//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::SetEmax(double Emax) {

	//Spectrum not yet set, do checks when setting it.
	if (fTotSpectrum == 0) {
		fMaxE = Emax;
		return;
	} else {
		double axisMax = fTotSpectrum->GetXaxis()->GetXmax();
		double axisMin = fTotSpectrum->GetXaxis()->GetXmin();
		if (Emax < axisMin) {
			LOG ("Flux", pWARN) << " fTotSpectrum axis min is: " << axisMin << " and the maximum energy was requested smaller than this: " << Emax << " IGNORE THIS REQUEST ";
			fMaxE = axisMax;
		} else {
			LOG ("Flux", pNOTICE) << "Setting E[max] = " << Emax;
			fMaxE = Emax;
		}

		//need to recompute the normalization
		fTotNorm = fTotSpectrum->Integral(fTotSpectrum->FindBin(fMinE), fTotSpectrum->FindBin(fMaxE));
		fEffEOTsPerDM = 1. / fTotNorm;
		LOG ("Flux", pNOTICE) << "SetEmax, fixing norm between: " << fMinE << " and " << fMaxE << " : " << fTotNorm;
		LOG ("Flux", pNOTICE) << "EOT per nu: " << fEffEOTsPerDM;
	}
}

//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::SetEnergyRange(double Emin, double Emax) {
	this->SetEmin(Emin);
	this->SetEmax(Emax);
}

//___________________________________________________________________________
//2D spectrum X: Y=SinTheta Y: Energy (GeV)
void GBDXpipesHistoFlux_DM::AddEnergySpectrum(int DM_pdgc, TH1D *spectrum, TH2D *spectrum2D) {
	LOG("Flux", pNOTICE) << "Adding flux spectrum for pdg = " << DM_pdgc;

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
	fPdgCList->push_back(DM_pdgc);

	bool accepted = (count(fPdgCList->begin(), fPdgCList->end(), DM_pdgc) == 1);
	if (!accepted) {
		LOG ("Flux", pWARN) << "The pdg-code isn't recognized and the spectrum was ignored";
	} else {

		if (spectrum->GetXaxis()->GetXmin() != spectrum2D->GetYaxis()->GetXmin()) {
			LOG ("Flux", pFATAL) << "histo has min x-axis range: " << spectrum->GetXaxis()->GetXmin() << " while histo2d has min y-axis range: " << spectrum2D->GetYaxis()->GetXmin() << " FATAL ERROR";
			exit(1);
		}
		if (spectrum->GetXaxis()->GetXmax() != spectrum2D->GetYaxis()->GetXmax()) {
			LOG ("Flux", pFATAL) << "histo has max x-axis range: " << spectrum->GetXaxis()->GetXmax() << " while histo2d has max y-axis range: " << spectrum2D->GetYaxis()->GetXmax() << " FATAL ERROR";
			exit(1);
		}

		fSpectrum.push_back(spectrum);
		fSpectrum2D.push_back(spectrum2D);

		//here do the trick on the spectrum (the integral was cached properly before this)
		int nbins = spectrum->GetNbinsX();
		for (int ii = 1; ii <= nbins; ii++) {
			double data = spectrum->GetBinContent(ii);
			double width = spectrum->GetBinWidth(ii);
			spectrum->SetBinContent(ii, data * width);
		}

		//prepare the projections - one projections for each y-bin of the 2D histogram
		TH1D **hArrayProjs = new TH1D*[nbins];
		TH1D *hProj;

		nbins = spectrum2D->GetNbinsY();
		for (int ii = 1; ii <= nbins; ii++) {
			hProj = (TH1D*) spectrum2D->ProjectionX(Form("%s_%i", spectrum2D->GetName(), ii), ii, ii);
			hArrayProjs[ii - 1] = (TH1D*) hProj->Clone();
			hArrayProjs[ii - 1]->SetName(Form("%s_%i", spectrum2D->GetName(), ii));
		}

		fSpectrum2D_proj.push_back(hArrayProjs);
		this->AddAllFluxes(); // update combined flux

	}
}

//___________________________________________________________________________
void GBDXpipesHistoFlux_DM::AddAllFluxes(void) {
	LOG("Flux", pNOTICE) << "Computing combined flux";

	if (fTotSpectrum) delete fTotSpectrum;

	vector<TH1D*>::const_iterator spectrum_iter;

	unsigned int inu = 0;
	fTotNorm = 0;
	for (spectrum_iter = fSpectrum.begin(); spectrum_iter != fSpectrum.end(); ++spectrum_iter) {
		TH1D *spectrum = *spectrum_iter;

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

	if (fMaxE < 0) {	//MAX - not yet set
		fMaxE = axisMax;
	} else {	// MAX was set. Check it is reasonable.
		if (fMaxE < axisMin) {
			LOG ("Flux", pWARN) << " fTotSpectrum axis min is: " << axisMin << " and the maximum energy was requested smaller than this: " << fMaxE << " - RESET Emax to max histo range " << axisMax;
			fMaxE = axisMax;
		}
	}

	if (fMinE < 0) {	//MIN - not yet set
		fMinE = axisMin;
	} else {	// MIN was set. Check it is reasonable.
		if (fMinE > axisMax) {
			LOG ("Flux", pWARN) << " fTotSpectrum axis max is: " << axisMax << " and the minimum energy was requested bigger than this: " << fMinE << " - RESET Emin to min histo range " << axisMin;
			fMinE = axisMin;
		}
	}

	fTotNorm = fTotSpectrum->Integral(fTotSpectrum->FindBin(fMinE), fTotSpectrum->FindBin(fMaxE)); //very important: do not use the width, since spectrum was fixed before.
	fEffEOTsPerDM = 1. / fTotNorm;
	if (fTotNorm <= 0) {
		LOG ("Flux", pERROR) << "fTotNorm is <0: " << fTotNorm;
		exit(1);
	} else {
		LOG("Flux",pINFO) << "fTotNorm is: " << fTotNorm;
		LOG ("Flux", pNOTICE) << "EOT per nu: " << fEffEOTsPerDM;
	}

}
//___________________________________________________________________________
int GBDXpipesHistoFlux_DM::SelectDM(double E) {
	const unsigned int n = fPdgCList->size();
	double fraction[n];

	vector<TH1D*>::const_iterator spectrum_iter;

	unsigned int iDM = 0;
	for (spectrum_iter = fSpectrum.begin(); spectrum_iter != fSpectrum.end(); ++spectrum_iter) {
		TH1D *spectrum = *spectrum_iter;
		fraction[iDM++] = spectrum->GetBinContent(spectrum->FindBin(E));
	}

	double sum = 0;
	for (iDM = 0; iDM < n; iDM++) {
		sum += fraction[iDM];
		fraction[iDM] = sum;
		LOG("Flux", pDEBUG) << "SUM-FRACTION(0->" << iDM << ") = " << sum;
	}

	RandomGen *rnd = RandomGen::Instance();
	double R = sum * rnd->RndFlux().Rndm();

	LOG("Flux", pDEBUG) << "R e [0,SUM] = " << R;

	for (iDM = 0; iDM < n; iDM++) {
		if (R < fraction[iDM]) return iDM;
	}

	LOG("Flux", pERROR) << "Could not select a DM species";
	assert(false);

	return -1;
}
//___________________________________________________________________________
double GBDXpipesHistoFlux_DM::GeneratePhi(void) const {
	RandomGen *rnd = RandomGen::Instance();
	double phi = 2. * kPi * rnd->RndFlux().Rndm(); // [0,2pi]
	return phi;
}
//___________________________________________________________________________
