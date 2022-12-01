#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolOps.h>

int main( int argc , char **argv ) {

  RDKit::ROMol *mol1 = RDKit::SmilesToMol( "CCC" );
  std::cout << "Number of atoms " << mol1->getNumAtoms() << std::endl;
  RDKit::ROMol *mol2 = RDKit::MolOps::addHs(*mol1);
  std::cout << "Number of atoms " << mol2->getNumAtoms() << std::endl;
  double *value = RDKit::MolOps::getAdjacencyMatrix(*mol1);
  // mol1->getProp("AdjacencyMatrix");
  for (int i=0; i < mol1->getPropList().size(); i++)
  {
    std::cout << mol1->getPropList()[i] << std::endl;
  };
}