#include <iostream>
#include <armadillo>

class Molecule
{
    private:
        std::string _smiles;

    public:
        Molecule(std::string smiles)
        {
            _smiles = smiles;

            int len = smiles.length();
            std::cout << len << std::endl;
        };

};

int main()
{
    Molecule("CCC");
}
