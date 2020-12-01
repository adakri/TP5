#ifndef _FINITEVOLUME_CPP
#include "FiniteVolume.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur
FiniteVolume::FiniteVolume(Function* function, DataFile* data_file, Mesh2D* mesh) :
_fct(function), _df(data_file), _msh(mesh)
{
	std::cout << "Build finite volume class." << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
}

// Construit la matrice des flux
void FiniteVolume::Build_flux_mat_and_rhs(const double& t)
{
	// Matrix
	_mat_flux.resize(_msh->Get_triangles().size(),_msh->Get_triangles().size());
	// RHS
	_BC_RHS.resize(_msh->Get_triangles().size());
	_BC_RHS.setZero();
	vector<Triplet<double>> triplets;	triplets.clear();

	vector<Edge> edges;
	edges= _msh -> Get_edges();
	int n=edges.size();

	vector<Triangle> triangles;
	triangles = _msh -> Get_triangles();
	std::cout << "------------------------------------------------------" << std::endl;
	VectorXd area= _msh ->Get_triangles_area();
	VectorXd bw_length= _msh -> Get_triangles_length();
	VectorXd edges_ln= _msh -> Get_edges_length();


	std::cout << "------------------------------------------------------" << std::endl;



	for (int i = 0; i < _msh->Get_edges().size(); i++)
	{
		int t1 = _msh->Get_edges()[i].Get_T1();
		int t2 =	_msh->Get_edges()[i].Get_T2();

		if (t2 != -1)
		{
					cout<<t1<<" "<<t2<<endl;
					double area1=area[t1];
					double area2=area[t2];
					double edge_ln=edges_ln[t2]; 
		}
		else
		{

		}



	}
	_mat_flux.setFromTriplets(triplets.begin(), triplets.end());
}

const MatrixXd FiniteVolume::Get_Matrix_A()
{
	double t=0.;
	FiniteVolume::Build_flux_mat_and_rhs(t);
	return _mat_flux;
}


// --- Déjà implémenté ---
// Construit la condition initiale au centre des triangles
VectorXd FiniteVolume::Initial_condition()
{
	VectorXd sol0(_msh->Get_triangles().size());

	for (int i = 0; i < _msh->Get_triangles().size(); i++)
	sol0(i) = _fct->Initial_condition(_msh->Get_triangles_center()(i,0),
	_msh->Get_triangles_center()(i,1));

	return sol0;
}

// Terme source au centre des triangles
VectorXd FiniteVolume::Source_term(double t)
{
	VectorXd sourceterm(_msh->Get_triangles().size());

	for (int i = 0; i < _msh->Get_triangles().size(); i++)
	{
		sourceterm(i) = _fct->Source_term(_msh->Get_triangles_center()(i,0),
		_msh->Get_triangles_center()(i,1), t);
	}
	return sourceterm;
}

// Solution exacte au centre des triangles
VectorXd FiniteVolume::Exact_solution(const double t)
{
	VectorXd exactsol(_msh->Get_triangles().size());

	for (int i = 0; i < _msh->Get_triangles().size(); i++)
	exactsol(i) = _fct->Exact_solution(_msh->Get_triangles_center()(i,0),
	_msh->Get_triangles_center()(i,1), t);

	return exactsol;
}

// Sauvegarde la solution
void FiniteVolume::Save_sol(const Eigen::VectorXd& sol, int n, std::string st)
{
	double norm = 0;
	for (int i = 0; i < sol.rows(); i++)
	norm += sol(i)*sol(i)*_msh->Get_triangles_area()[i];
	norm = sqrt(norm);

	if (st == "solution")
	{
		cout << "Norme de u = " << norm << endl;
	}

	string name_file = _df->Get_results() + "/" + st + "_" + std::to_string(n) + ".vtk";
	int nb_vert = _msh->Get_vertices().size();
	assert((sol.size() == _msh->Get_triangles().size())
	&& "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((_msh->Get_vertices()[i]).Get_coor())[0] << " "
		<< ((_msh->Get_vertices()[i]).Get_coor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << _msh->Get_triangles().size() << " "
	<< _msh->Get_triangles().size()*4 << endl;
	for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
	{
		solution << 3 << " " << ((_msh->Get_triangles()[i]).Get_vertices())[0]
		<< " " << ((_msh->Get_triangles()[i]).Get_vertices())[1]
		<< " " << ((_msh->Get_triangles()[i]).Get_vertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << _msh->Get_triangles().size() << endl;
	for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << _msh->Get_triangles().size() << endl;
	solution << "SCALARS sol float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	double eps = 1.0e-10;
	for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,sol[i]) << endl;
	}
	solution << endl;

	//solution << "CELL_DATA " << _msh->Get_triangles().size() << endl;
	solution << "SCALARS CFL float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,_df->Get_dt()*fabs(sol[i])/_msh->Get_triangles_length()(i)) << endl;
	}
	solution << endl;

	if (_df->Get_mu() > 1e-10)
	{
		solution << "SCALARS Pe float 1" << endl;
		solution << "LOOKUP_TABLE default" << endl;
		// To avoid strange behaviour (which appear only with Apple)
		// with Paraview when we have very small data (e-35 for example)
		for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
		{
			solution << max(eps,_msh->Get_triangles_length()(i)*fabs(sol[i])/_df->Get_mu()) << endl;
		}
		solution << endl;
	}

	solution.close();
}

#define _FINITEVOLUME_CPP
#endif
