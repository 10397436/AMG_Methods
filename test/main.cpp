/**
* @file   main.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include  "common.h"
#include  "GetPot.h"
#include  "parameter_setup.h"
#include  "parameter_cycle.h"
#include "parameter_method.h"
#include  "setup.h"
#include  "setupDG.h"
#include  "cycle.h"
#include  "method.h"
#include  "output.h"

/**
* @brief The @b main function.
*/

int main(const int argc, char * argv[]) {

GetPot commandLine(argc, argv);
const string config_directory=commandLine.follow("../config.pot",2,"-f","--file");
GetPot config(config_directory.c_str());

/**
* Read input/output and files parameters.
*/

const string inputdirectory=config("inputdirectory","share/examples/");
const string outputdirectory=config("outputdirectory","share/results/");
const string inputA=config("inputA","error");
const string inputf=config("inputf","error");
const string testname=config("outputTest","Test1.txt");
const string print=config("print","screen");
const string savesol=config("savesol","N");

if(inputA=="error" || inputf=="error")
{
	throw invalid_argument("Received invalid argument: give name of file.");
}

SpMat A;
loadMarket(A,inputdirectory+inputA);
A.makeCompressed();
Vec f;
loadMarketVector(f,inputdirectory+inputf);

if(A.rows()!=f.size())
{
	throw runtime_error("Input matrix and vector of different sizes.");
}

SimplicialLLT<SpMat> lltOfA(A); // compute the Cholesky decomposition of A
if(lltOfA.info() == NumericalIssue)
{
	throw runtime_error("Possibly non symmetric semi-positive definitie matrix.");
} 

/**
* Read setup parameters.
*/

const Real theta=config("theta",0.25);
const int nlevel=config("nlevel",2)-1;

if(theta<=0 || theta>1 || nlevel<1)
{
	throw invalid_argument("Received invalid argument: check setup parameters.");
}

/**
* Read cycle parameters.
*/

const int nu1=config("nu1",1);
const int nu2=config("nu2",1);
const int mu=config("mu",1);

if(nu1<1 || nu2<1 || mu<1 || mu>2)
{
	throw invalid_argument("Received invalid argument: check cycle parameters.");
}

/**
* Read method parameters.
*/

const Real tol=config("tol",1e-8);
const int nmaxiter=config("nmaxiter",150);
const string fem=config("fem","CG");
const string multigrid=config("method","AMG");

if(tol<0 || nmaxiter<0)
{
	throw invalid_argument("Received invalid argument: check method parameters.");
}


parameter_setup ps(nlevel,theta);
parameter_cycle pc(nlevel,nu1,nu2,mu);
parameter_method pm(tol,nmaxiter);

/**
* Instantiate setup.
*/

setup *S;

if(fem=="CG")
{
	cout<<"CG setup"<<endl;
	S= new setup(A,ps);
}
else if(fem=="DG")
{
	cout<<"DG setup"<<endl;
	S= new setupDG(A,ps);
}
else
{
	throw invalid_argument("Received invalid argument: check method parameters.");
}

/**
* Instantiate cycle.
*/

cycle C(*S, f, pc);

/**
* Instantiate method.
*/

method M(C,pm);

/**
* Apply AMG.
*/

if(multigrid=="AMG")
{
	cout<<"Running AMG stand-alone."<<endl;
	M.AMGCycle();
}
else if(multigrid=="PCG")
{
	cout<<"Running AMG as preconditioner for conjugate gradient."<<endl;
	M.PCGCycle();
}
else
{
	throw invalid_argument("Received invalid argument: check method parameters.");
}

/**
* Print output.
*/

output O(testname,inputA,inputf,fem,multigrid,M.get_iter(),M.get_rho(),M.get_flag(),ps,pc,pm);

if(print=="file")
{
	cout<<"Printing on file."<<endl;
	O.print_on_file(outputdirectory);
}
else if(print=="screen")
	O.print_on_screen();
else
{
	throw invalid_argument("Received invalid argument: check output flags.");
}

/**
* Save solution.
*/

if(savesol=="Y")
{
	cout<<"Saving solution on file."<<endl;
	saveMarketVector(M.get_solution(),outputdirectory+"sol_"+testname);
}
else if(savesol=="N")
	cout<<""<<endl;
else
{
	throw invalid_argument("Received invalid argument: check output flags.");
}

return 0;
}
