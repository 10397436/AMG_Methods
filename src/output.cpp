/**
* @file   output.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include "output.h"

output::output(const string& testname,const string& inputA, const string& inputf, const string& fem, const string& method, const int& iter,const Real& rho, const bool& flag, const parameter_setup& ps, const parameter_cycle& pc, const parameter_method& pm)
{
_testname=testname;
_inputA=inputA;
_inputf=inputf;
_fem=fem;
_method=method;
_iter=iter;
_rho=rho;
_flag=flag;
_ps=ps;
_pc=pc;
_pm=pm;
}

void output::print_on_screen() const
{
cout<<endl;
cout<<"FILES"<<endl;
cout<<"Matrix A: "<<_inputA<<endl;
cout<<"Vector f: "<<_inputf<<endl;
cout<<endl;
cout<<"PARAMETERS"<<endl;
cout<<"nmatrix = "<<_ps.get_nmatrix()+1<<endl;
cout<<"theta = "<<_ps.get_theta()<<endl;
cout<<"nlevel = "<<_pc.get_nlevel()+1<<endl;
cout<<"nu1 = "<<_pc.get_nu1()<<endl;
cout<<"nu2 = "<<_pc.get_nu2()<<endl;
cout<<"mu = "<<_pc.get_mu()<<endl;
cout<<"tol = "<<_pm.get_tol()<<endl;
cout<<"maxiter = "<<_pm.get_maxiter()<<endl;
cout<<"fem = "<<_fem<<endl;
cout<<"method = "<<_method<<endl;
cout<<endl;
cout<<"RESULTS"<<endl;
if(_flag==1)
	cout<<"Method not convergent"<<endl;
else
{
	cout << "Method converges in "<<_iter<<" iterations"<<endl;
	cout << "Method has convergence factor rho = "<<_rho<<endl;
}
}

void output::print_on_file(const string& directory) const
{
ofstream myfile;
myfile.open(directory+_testname);
myfile<<"FILES"<<endl;
myfile<<"Matrix A: "<<_inputA<<endl;
myfile<<"Vector f: "<<_inputf<<endl;
myfile<<endl;
myfile<<"PARAMETERS"<<endl;
myfile<<"nmatrix = "<<_ps.get_nmatrix()+1<<endl;
myfile<<"theta = "<<_ps.get_theta()<<endl;
myfile<<"nlevel = "<<_pc.get_nlevel()+1<<endl;
myfile<<"nu1 = "<<_pc.get_nu1()<<endl;
myfile<<"nu2 = "<<_pc.get_nu2()<<endl;
myfile<<"mu = "<<_pc.get_mu()<<endl;
myfile<<"tol = "<<_pm.get_tol()<<endl;
myfile<<"maxiter = "<<_pm.get_maxiter()<<endl;
myfile<<"fem = "<<_fem<<endl;
myfile<<"method = "<<_method<<endl;
myfile<<endl;
myfile<<"RESULTS"<<endl;
if(_flag==1)
	myfile<<"Method not convergent"<<endl;
else
{
	myfile << "Method converges in "<<_iter<<" iterations"<<endl;
	myfile << "Method has convergence factor rho = "<<_rho<<endl;
}
myfile.close();
}
