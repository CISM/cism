#include "trilinosModelEvaluator.hpp"

extern "C" {
  void calc_F(double* x, double* f, int N, void* bb);
  void apply_precond_nox(double* x, double* y, int n, void* bb);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

trilinosModelEvaluator::trilinosModelEvaluator (
                            int N_, double* statevector,
                            const Epetra_Comm& comm_, void* blackbox_res_)
			     : N(N_), comm(comm_), blackbox_res(blackbox_res_)
{
  
  cout << "In TrilinosModelEvaluator constructor:" << endl;
  xMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));

  cout << "  proc  = " << comm.MyPID() <<  ",  N = " << N 
       << ",  Ntot = " << xMap->NumGlobalElements() << endl;

  xVec = Teuchos::rcp(new Epetra_Vector(Copy, *xMap, statevector));

  precOp = Teuchos::rcp(new trilinosPreconditioner(N, xVec, xMap, blackbox_res));
}

/*******************************************************************************/
// NTS: Assuming that x, f and g all have the same map for all entries.
// This is silly, but sufficient for now.
// Return solution vector map
Teuchos::RCP<const Epetra_Map> trilinosModelEvaluator::get_x_map() const{
  return xMap;
}

// Return residual vector map
Teuchos::RCP<const Epetra_Map> trilinosModelEvaluator::get_f_map() const{
  return xMap;
}

// Return initial solution and x_dot init
Teuchos::RCP<const Epetra_Vector> trilinosModelEvaluator::get_x_init() const{
  return xVec;
}

Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
trilinosModelEvaluator::create_WPrec() const
{
  // bool is answer to: "Prec is already inverted?"
  return Teuchos::rcp(new EpetraExt::ModelEvaluator::Preconditioner(precOp,true));
}


/*******************************************************************************/
// Create InArgs
EpetraExt::ModelEvaluator::InArgs trilinosModelEvaluator::createInArgs() const{
  InArgsSetup inArgs;

  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(0);
  return inArgs;
}

/*******************************************************************************/
// Create OutArgs
EpetraExt::ModelEvaluator::OutArgs trilinosModelEvaluator::createOutArgs() const{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(0, 0);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_WPrec, true);

  return outArgs;
}

/*******************************************************************************/
// Evaluate model on InArgs
void trilinosModelEvaluator::evalModel(const InArgs& inArgs, const OutArgs& outArgs) const{

  //double nrm;
  // Get the solution vector x from inArgs and residual vector from outArgs
  RCP<const Epetra_Vector> x = inArgs.get_x();
  RCP<Epetra_Vector> f = outArgs.get_f();
  
  if (x == Teuchos::null) throw "trilinosModelEvaluator::evalModel: x was NOT specified!";

  // Save the current solution, which makes it initial guess for next nonlienar solve
  *xVec = *x;

  if (f != Teuchos::null) {
    f->PutScalar(0.0);
    calc_F(x->Values(), f->Values(), N, blackbox_res);

    //f->Norm2(&nrm);
    //cout << "AGS  Resid norm in eval_model total " << nrm << endl;
  }

  RCP<Epetra_Operator> WPrec = outArgs.get_WPrec();
  if (WPrec != Teuchos::null) {
     //cout << "evalModel called for WPrec -- doing nothing " <<  endl;
  }
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
trilinosPreconditioner::trilinosPreconditioner (
       int N_, RCP<Epetra_Vector> xVec_, RCP<Epetra_Map> xMap_, void* blackbox_res_)
       : N(N_), xVec(xVec_), xMap(xMap_), blackbox_res(blackbox_res_)
{
}

int trilinosPreconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  apply_precond_nox(Y(0)->Values(), X(0)->Values(), N, blackbox_res);

  return 0;
}


