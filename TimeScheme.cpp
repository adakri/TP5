#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, FiniteVolume* adv) :
_df(data_file), _fin_vol(adv), _t(_df->Get_t0()), _sol(adv->Initial_condition())
{
}

EulerScheme::EulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
}

ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
  std::cout << "Build time scheme class." << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{
}

// Renvoie _sol (pratique pour vérifier la résolution)
const VectorXd & TimeScheme::Get_sol() const
{
  return _sol;
}

// Euler Explicite
void EulerScheme::Advance()
{
  //construction de A
  _fin_vol -> Build_flux_mat_and_rhs(_t);
  MatrixXd A=_fin_vol -> Get_flux_matrix();
  //construction de b
  VectorXd b=_fin_vol ->Get_BC_RHS();

  cout<<b<<endl;

  double dt= _df -> Get_dt();
  VectorXd S=_fin_vol ->Source_term(_t);
  cout<<S<<endl;

  _sol=+ dt*S -dt*(A*_sol+b);



}

// Euler Implicite
void ImplicitEulerScheme::Advance()
{
  //construction de A
  _fin_vol -> Build_flux_mat_and_rhs(_t);
  MatrixXd A=_fin_vol -> Get_flux_matrix();
  //construction de b
  VectorXd b=_fin_vol ->Get_BC_RHS();
;
  cout<<b<<endl;

  double dt= _df -> Get_dt();

  VectorXd S=_fin_vol ->Source_term(_t+dt);

  cout<<S<<endl;

  //construction de I
  SparseMatrix<double> I(_sol.size(),_sol.size());
  vector<Triplet<double>> triplets;	triplets.clear();
  for (int i=0; i<_sol.size(); i++)
  {
      triplets.push_back({i,i,1.});
  }
  I.setFromTriplets(triplets.begin(), triplets.end());


  //cout<<I<<endl;

  SparseMatrix<double> B= I + dt*A;


  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;

  solver.analyzePattern(B);
  // Compute the numerical factorization
  solver.factorize(B);
  //Use the factors to solve the linear system
  _sol = solver.solve(_sol - dt*b + dt*S);

}

#define _TIME_SCHEME_CPP
#endif
