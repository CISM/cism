#include <iostream>
#include "matrixInterface.hpp"

// Constructor
TrilinosMatrix_Interface::TrilinosMatrix_Interface
  (const Teuchos::RCP<const Epetra_Map>& rowMap,
   int bandwidth, const Epetra_Comm& comm)
  : rowMap_(rowMap), bandwidth_(bandwidth), matrixOrder_(-1), comm_(comm) {
  
  matrixOrder_ = rowMap->NumGlobalElements();

  operator_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *rowMap, bandwidth) );
  isFillCompleted_ = false;

  // create map of full vector
/*
  fullMap_ =  Teuchos::rcp(new Epetra_LocalMap(matrixOrder_, 1, comm));

  import_r2f = Teuchos::rcp(new Epetra_Import(*fullMap_, *rowMap_));
  import_f2r = Teuchos::rcp(new Epetra_Import(*rowMap_, *fullMap_));
*/
}

// Destructor
TrilinosMatrix_Interface::~TrilinosMatrix_Interface() {
}

// Accessor methods
bool TrilinosMatrix_Interface::isSparsitySet() const {return isFillCompleted_;}
int TrilinosMatrix_Interface::bandwidth() const {return bandwidth_;}
int TrilinosMatrix_Interface::matrixOrder() const {return matrixOrder_;}
const Epetra_Map& TrilinosMatrix_Interface::getRowMap() const {return *rowMap_;}
Teuchos::RCP<Epetra_CrsMatrix>& TrilinosMatrix_Interface::getOperator() {return operator_;}


// Fix the sparsity patter by calling FillComplete
void TrilinosMatrix_Interface::finalizeSparsity() {
  isFillCompleted_ = true;
  int ierr = operator_->FillComplete();
  assert (ierr==0);
}

// Update the operator and also the corresponding row map.
void TrilinosMatrix_Interface::updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator) {
  operator_ = newOperator;
  isFillCompleted_ = operator_->Filled();
}

Teuchos::RCP<Epetra_Vector> TrilinosMatrix_Interface::getPartitionedVec(double *fullRhs)
{
  if (fullMap_== Teuchos::null) fullMap_ =  Teuchos::rcp(new Epetra_LocalMap(matrixOrder_, 1, comm_));
  if (import_f2r == Teuchos::null) import_f2r = Teuchos::rcp(new Epetra_Import(*rowMap_, *fullMap_));

  Teuchos::RCP<Epetra_Vector> rhs_EV
    = Teuchos::rcp(new Epetra_Vector(*rowMap_));
  Epetra_Vector fullRhs_EV(View, *fullMap_, fullRhs); 
  rhs_EV->Import(fullRhs_EV, *import_f2r, Insert);
  return rhs_EV;
}

void TrilinosMatrix_Interface::spreadVector(const Epetra_Vector& vec, double* fullVec)
{
  if (fullMap_== Teuchos::null) fullMap_ =  Teuchos::rcp(new Epetra_LocalMap(matrixOrder_, 1, comm_));
  if (import_r2f == Teuchos::null) import_r2f = Teuchos::rcp(new Epetra_Import(*fullMap_, *rowMap_));

  Epetra_Vector fullVec_EV(View, *fullMap_, fullVec); 
  fullVec_EV.Import(vec, *import_r2f, Insert);
}

