////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Jun 10 09:21:35 EDT 2016
//             From a coffee shop in Brooklyn, NY
//
////////////////////////////////////////////////////////////////////

#include "TParticleBeamContainer.h"
#include "TRandomA.h"

// Random number defined elsewhere
extern TRandomA* gRandomA;


TParticleBeamContainer::TParticleBeamContainer ()
{
  // Constructor
}





TParticleBeamContainer::~TParticleBeamContainer ()
{
  // Destruction
}





TParticleBeam& TParticleBeamContainer::AddNewParticleBeam (std::string const& Type, std::string const& Name, TVector3D const& X0, TVector3D const& D0, double const E0, double const T0, double const Current, double const Weight, double const Charge, double const Mass)
{
  std::string const MyName = Name != "" ? Name : "_beam" + std::to_string(fParticleBeams.size());

  if (fParticleBeamMap.count(MyName) != 0) {
    std::cerr << "fParticleBeamMap.count(Name) != 0" << std::endl;
    throw std::invalid_argument("beam with this name already exists");
  }

  if (fParticleBeamWeightSums.size() == 0) {
    fParticleBeamWeightSums.push_back(Weight);
  } else {
    fParticleBeamWeightSums.push_back(fParticleBeamWeightSums.back() + Weight);
  }

  if (Type == "custom") {
    fParticleBeams.push_back( TParticleBeam(Type, MyName, X0, D0, E0, T0, Current, Charge, Mass, Weight) );
  } else {
    fParticleBeams.push_back( TParticleBeam(Type, MyName, X0, D0, E0, T0, Current, Weight) );
  }
  fParticleBeamMap[MyName] = fParticleBeams.size() - 1;

  return fParticleBeams[fParticleBeams.size() - 1];
}




TParticleBeam& TParticleBeamContainer::AddNewParticleBeam (std::string const& Beam, std::string const& Name, double const Weight)
{
  std::string const MyName = Name != "" ? Name : "_beam" + std::to_string(fParticleBeams.size());

  if (fParticleBeamMap.count(MyName) != 0) {
    std::cerr << "fParticleBeamMap.count(Name) != 0" << std::endl;
    throw std::invalid_argument("beam with this name already exists");
  }

  if (fParticleBeamWeightSums.size() == 0) {
    fParticleBeamWeightSums.push_back(Weight);
  } else {
    fParticleBeamWeightSums.push_back(fParticleBeamWeightSums.back() + Weight);
  }

  fParticleBeams.push_back( TParticleBeam(Beam, MyName, Weight) );

  fParticleBeamMap[MyName] = fParticleBeams.size() - 1;

  return fParticleBeams[fParticleBeams.size() - 1];
}




TParticleA TParticleBeamContainer::GetNewParticle ()
{
  // UPDATE: incomplete
  return fParticleBeams[ this->GetRandomBeamIndexByWeight() ].GetNewParticle();
}




TParticleBeam& TParticleBeamContainer::GetParticleBeam (size_t const i)
{
  // Return a reference to the particle beam given its name
  if (i >= fParticleBeams.size()) {
    throw std::length_error("beam index out of range");
  }

  return fParticleBeams[i];
}



TParticleBeam& TParticleBeamContainer::GetParticleBeam (std::string const& Name)
{
  // Return a reference to the particle beam given its name.  If "" is given returns
  // a random beam beased on weights

  if (Name == "") {
    return this->GetRandomBeam();
  }

  if (fParticleBeamMap.count(Name) == 0) {
    throw std::out_of_range("beam name not in map");
  }

  return this->GetParticleBeam(fParticleBeamMap[Name]);
}




TParticleBeam& TParticleBeamContainer::GetRandomBeam ()
{
  return this->GetParticleBeam(this->GetRandomBeamIndexByWeight());
}




size_t TParticleBeamContainer::GetRandomBeamIndexByWeight () const
{
  // UPDATE: Get a better random generator with seed setting elsewhere

  // Size of array
  size_t const N = fParticleBeamWeightSums.size();

  // If it's zero we don't really know what we are doing here..
  if (N == 0) {
    throw std::length_error("no beam defined");
  }

  // If we're 1, that's easy
  if (N == 1) {
    return 0;
  }

  // Get a random double [0, SumOfWeights)
  double const Random = gRandomA->Uniform() * fParticleBeamWeightSums[N - 1];

  // Not the fastest algorithm, but I guess you don't have thousands of different beams...
  // If you do, let's update this search...
  for (size_t i = 0; i != N; ++i) {
    if (Random < fParticleBeamWeightSums[i]) {
      return i;
    }
  }

  // Just in case you don't find it, something is seriously wrong..
  std::cerr << "ERROR: TParticleBeamContainer::GetRandomBeamIndexByWeight did not find a beam for this weight" << std::endl;
  throw std::out_of_range("random weight out of range.  SERIOUS ERROR");

  return 0;
}




size_t TParticleBeamContainer::GetNParticleBeams () const
{
  // Return the number of particle beams
  return fParticleBeams.size();
}





void TParticleBeamContainer::Clear ()
{
  // Clear the particle beam container contents

  fParticleBeamWeightSums.clear();
  fParticleBeams.clear();
  fParticleBeamMap.clear();

  return;
}






void TParticleBeamContainer::SetEmittance (std::string const& Beam,
                                           TVector2D const& Emittance)
{
  // Set the emittance for any beam matching exactly the name given
  // or to all beams if "" is given


  for (std::vector<TParticleBeam>::iterator it = fParticleBeams.begin(); it != fParticleBeams.end(); ++it) {
    if (Beam == "") {
      it->SetEmittance(Emittance);
    } else if (Beam == it->GetName()) {
      it->SetEmittance(Emittance);
    }
  }

  return;
}



void TParticleBeamContainer::SetTwissParameters (std::string const& Beam,
                                                 TVector2D const& Beta,
                                                 TVector2D const& Alpha,
                                                 TVector2D const& Gamma,
                                                 TVector3D const& Lattice_Reference,
                                                 bool const HasReferencePoint)
{
  // Set the twiss parameters for any beam matching exactly the name given
  // or to all beams if "" is given

  for (std::vector<TParticleBeam>::iterator it = fParticleBeams.begin(); it != fParticleBeams.end(); ++it) {
    if (Beam == "") {
      it->SetTwissParameters(Beta, Alpha, Gamma, Lattice_Reference, HasReferencePoint);
    } else if (Beam == it->GetName()) {
      it->SetTwissParameters(Beta, Alpha, Gamma, Lattice_Reference, HasReferencePoint);
    }
  }

  return;
}
