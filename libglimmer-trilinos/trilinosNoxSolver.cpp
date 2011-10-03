// Trilinos Objects
#include "Piro_Epetra_NOXSolver.hpp"
#include "trilinosModelEvaluator.hpp"

#include "Epetra_MpiComm.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Teuchos_DefaultMpiComm.hpp"

#include "config.inc"

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
void FC_FUNC(noxinit,NOXINIT) ( int* nelems, double* statevector,
               int* mpi_comm_ignored, void* blackbox_res)
{

 try {

  // Build the epetra communicator
  Comm_=rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  Epetra_Comm& Comm=*Comm_;
  printProc = (Comm_->MyPID() == 0);
  Teuchos::MpiComm<int> tcomm(Teuchos::opaqueWrapper(MPI_COMM_WORLD));
  
  if (printProc) cout << "NOXINIT CALLED    for nelem=" << *nelems << endl;


    try { // Check that the parameter list is valid at the top
      RCP<Teuchos::ParameterList> pl =
        rcp(new Teuchos::ParameterList("Trilinos Options for NOX"));
      Teuchos::updateParametersFromXmlFileAndBroadcast(
                             "trilinosOptions.xml", pl.get(),tcomm);
 
      Teuchos::ParameterList validPL("Valid List");;
      validPL.sublist("Stratimikos"); validPL.sublist("Piro");
      pl->validateParameters(validPL, 0);
      paramList = Teuchos::sublist(pl,"Piro",true);
    }
    catch (std::exception& e) {
      cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n"
           << e.what() << "\nExiting: Invalid trilinosOptions.xml file."
           << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
      exit(1);
    }
  
  paramList->set("Lean Matrix Free",true); // Saves some GMRES steps
  if (printProc) cout << "NOXInit: param list is: (delete this debug line)\n" << *paramList << endl;

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
void FC_FUNC(noxsolve,NOXSOLVE) (int* nelems, double* statevector, void* blackbox_res)
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
void FC_FUNC(noxfinish,NOXFINISH) (void)
{
 if (printProc) cout << "NOXFinish called" << endl;

 // Free memory
 Nsolver   = Teuchos::null;
 model     = Teuchos::null;
 paramList = Teuchos::null;
 Comm_     = Teuchos::null;
}

} //extern "C"
