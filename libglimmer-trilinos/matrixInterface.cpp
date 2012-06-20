#include <iostream>
#include "Teuchos_TestForException.hpp"
#include "matrixInterface.hpp"

// Constructor
TrilinosMatrix_Interface::TrilinosMatrix_Interface
  (const Teuchos::RCP<const Epetra_Map>& rowMap,
   int bandwidth, const Epetra_Comm& comm)
  : rowMap_(rowMap), bandwidth_(bandwidth), matrixOrder_(-1), comm_(comm) {
  
  matrixOrder_ = rowMap->NumGlobalElements();

  operator_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *rowMap, bandwidth) );
  isFillCompleted_ = false;
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
  TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::logic_error,
     "Error: Trilinos Fill Complete  returned nozero error code ( " << ierr << " )\n");

}

// Update the operator and also the corresponding row map.
void TrilinosMatrix_Interface::updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator) {
  operator_ = newOperator;
  isFillCompleted_ = operator_->Filled();
}
