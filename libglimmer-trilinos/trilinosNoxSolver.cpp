// Trilinos Objects
#include "Piro_Epetra_NOXSolver.hpp"
#include "trilinosModelEvaluator.hpp"

#include "Epetra_MpiComm.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;

// Objects that are global to the file
static RCP<Piro::Epetra::NOXSolver> Nsolver;
static RCP<trilinosModelEvaluator> model;
static RCP<Teuchos::ParameterList> paramList;
static RCP<Epetra_MpiComm> Comm_;

static EpetraExt::ModelEvaluator::InArgs inArgs;
static EpetraExt::ModelEvaluator::OutArgs outArgs;
static bool printProc;


extern "C" {
void noxinit_( int* nelems, double* statevector,
               int* mpi_comm_ignored, void* blackbox_res)
{

 try {

  // Build the epetra communicator
  Comm_=rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  Epetra_Comm& Comm=*Comm_;
  printProc = (Comm_->MyPID() == 0);
  
  if (printProc) cout << "NOXINIT CALLED    for nelem=" << *nelems << endl;

  paramList = rcp(new Teuchos::ParameterList);

  
  Teuchos::updateParametersFromXmlFile("input.xml", paramList.get());
  paramList->set("Lean Matrix Free",true); // Saves some GMRES steps
  if (printProc)
  cout << "NOXInit: param list from input.xml is:\n" << *paramList << endl;

  model = rcp(new trilinosModelEvaluator(*nelems, statevector, Comm, blackbox_res));
    
  Nsolver = rcp(new Piro::Epetra::NOXSolver(paramList, model));

  inArgs=Nsolver->createInArgs();
  outArgs=Nsolver->createOutArgs();

  // Ask the model for the converged solution from g(0)
  RCP<const Epetra_Map> xmap = Nsolver->get_g_map(0);
  RCP<Epetra_Vector> xout = rcp(new Epetra_Vector(*xmap));

  outArgs.set_g(0,xout);

 } //end try block
  catch (std::exception& e) { cout << e.what() << endl; exit(1); }
  catch (const char *s) { cout << s << endl; exit(1); }
  catch (...) { cout << "Caught unknown exception!" << endl; exit(1); }
}

/****************************************************/
void noxsolve_(int* nelems, double* statevector, void* blackbox_res)
{
  try {
    TEST_FOR_EXCEPTION(Nsolver==Teuchos::null, logic_error, 
                          "Exception: noxsolve called with solver=null: \n"
       << "You either did not call noxinit first, or called noxfinish already");
    if (printProc) cout << "NOXSolve called" << endl;

    // Solve    
    Nsolver->evalModel(inArgs,outArgs);

    // Copy out the solution
    RCP<Epetra_Vector> xout = outArgs.get_g(0); 
    if(xout == Teuchos::null) throw "evalModel is NOT returning a vector";

    for (int i=0; i<*nelems; i++) statevector[i] = (*xout)[i];
  }
  catch (std::exception& e) { cout << e.what() << endl; exit(1); }
  catch (const char *s) { cout << s << endl; exit(1); }
  catch (...) { cout << "Caught unknown exception!" << endl; exit(1); }

}

/****************************************************/ 
void noxfinish_()
{
 if (printProc) cout << "NOXFinish called" << endl;

 // Free memory
 Nsolver   = Teuchos::null;
 model     = Teuchos::null;
 paramList = Teuchos::null;
 Comm_     = Teuchos::null;
}

} //extern "C"
