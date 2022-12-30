#pragma once
#include "input.h"
#include "utility.h"
#include "output.h"

class calcEnergy
{
    private:
        Molecule* _Mol;
        
        /**
        * Maps atom types to bond stretching parameters
        * First pair is the atom types involved in the bond, and the second pair are the params
        * First element of the second pair is the kb, while the second element of the second pair is the r0
        * Need 1,5 and 5,1 because C-H and H-C when considering the stretch bending contributions.
        **/
        std::map<std::pair<int, int>, std::pair<double, double>> bs_params = { {{1,1} , {4.258, 1.508}}, {{1,5} , {4.766, 1.093}}, {{5,1} , {4.766, 1.093}} };
        
        /**
        * Maps atom types to angle parameters
        * First tuple is the atom types involved in the angle, and the second pair are the params
        * First element of the second pair is the ka, while the second element of the second pair is the a0
        **/
        std::map<std::tuple<int, int, int>, std::pair<double, double>> ab_params = { {{1,1,5} , {0.636, 110.549}}, {{5, 1, 5}, {0.516, 108.836}} };
        
        /**
        * Maps atom types to stretch bond parameters
        * First tuple is the atom types involved in the interaction and second pair are the params
        * First element of the second pair is the kbaIJK, while the second element of the second pair is the kbaKJI
        **/
        std::map<std::tuple<int, int, int>, std::pair<double, double>> asb_params = { {{1,1,5} , {0.227, 0.070}}, {{5, 1, 5}, {0.115, 0.115}} };
        
        /**
        * Maps atom types to torsion parameters
        * First tuple is the atom types involved in the interaction and the second pair are the params
        * Elements of the second tuple are V1, V2, V3 respectively
        **/
        std::map<std::tuple<int, int, int, int>, std::tuple<double, double, double>> t_params = { {{5,1,1,5} , {0.284, -1.386, 0.314}} };
        
        /**
        * Maps atom types to vdw parameters
        * Int is the atom type involved in the interaction
        * Elements of the second tuple are alpha, N, A and G, respectively
        **/
        std::map<int, std::tuple<double, double, double, double>> vdw_params = { {1, {1.050, 2.490, 3.890, 1.282}}, {5, {0.250, 0.800, 4.200, 1.209}} };


    public:
        
        /**
         * @brief calcEnergy Constructor
         * @details Holds all the molecule information (atom types, bonds, angles, etc) and computes the total energy of the molecule based on the 
         * current coordinates. 
         * 
         * @param Mol (Molecule*): Pointer to an object of class molecule
         **/
        calcEnergy(Molecule* Mol);

        /**
         * @brief Calculates and returns the bond stretching energy between 2 atoms, i and j.
         * 
         * @param kb (double): Bond stretching force constant, defined per set of atoms i and j.
         * @param delta_R (double): Difference between actual and reference bond lengths.
         * @return double
         **/
        double bondStretching(double kb, double delta_R);

        /**
         * @brief Calculates the angle bending between atoms i, j and k, where j is the middle atom.
         * 
         * @param ka (double): Angle bending force constant, defined per set of atoms i, j, and k.
         * @param delta_Angle (double): Difference between actual and reference angle.
         * @return double
         **/
        double angleBending(double ka, double delta_Angle);

        /**
         * @brief Calculates the angle stretch bending between atoms i, j, k.
         * 
         * @param kba_ijk (double): Stretch bending force constant, defined per set of atoms i, j, and k.
         * @param kba_kji (double): Stretch bending force constant, defined per set of atoms k, j, and i.
         * @param delta_Rij (double): Difference between actual and reference bond lengths between atoms i and j.
         * @param delta_Rkj (double): Difference between actual and reference bond lengths between atoms j and k.
         * @param delta_Angle (double): Difference between actual and reference angle.
         * @return double
         **/
        double angleStretchBending(double kba_ijk, double kba_kji, double delta_Rij, double delta_Rkj, double delta_Angle);

        /**
         * NOT IMPLEMENTED 
         * Need helper function to calculate the wilson angle (X) between the bond j-l, and the plane i-j-k.
         * Calculates the out of plane bending at atom j in i,j,k.
         **/
        double outOfPlaneBending(double koop, double X);

        /**
         * @brief Calculates the torsional energy based on atoms i, j, k, l.
         * 
         * @param V1 (double): Constant defined per set of atoms i, j, k, l.
         * @param V2 (double): Constant defined per set of atoms i, j, k, l.
         * @param V3 (double): Constant defined per set of atoms i, j, k, l.
         * @param omega (double): Difference between actual and reference dihedral angles between atoms i, j, k, and l.
         * @return double
         **/
        double torsion(double V1, double V2, double V3, double omega);

        /**
         * @brief Calculates the vanderwaals interaction based on the buffered-14-7 form, based on atoms i and j.
         * 
         * @param dist_Rij (double): Distance between atoms i and j, with a minimum of 3 bonds away.
         * @param Ai (double): Constant defined per set of atoms i, j.
         * @param Aj (double): Constant defined per set of atoms i, j.
         * @param alphai (double): Constant defined per set of atoms i, j.
         * @param alphaj (double): Constant defined per set of atoms i, j.
         * @param Ni (double): Constant defined per set of atoms i, j.
         * @param Nj (double): Constant defined per set of atoms i, j.
         * @param Gi (double): Constant defined per set of atoms i, j.
         * @param Gj (double): Constant defined per set of atoms i, j.         * 
         * @return double
         **/
        double vdw(double dist_Rij, double Ai, double Aj, double alphai, double alphaj, double Ni, double Nj, double Gi, double Gj);

        /**
         * NOT IMPLEMENTED
         * Need helper function to calculate the partial atomic charges (qi, qj) based on bond charge increments (bci)
         * Calculates the electrostatic interaction between atoms i, j.
         */

        double electrostatic(double delta_R, double qi, double qj);

        /**
         * @brief Sums all the bonding contributions in the molecule.
         * 
         * @param molecule_coordinates (arma::mat): Coordinates of the molecule.
         * @param verbose (bool): Whether to output verbose information.
         * @return double 
         */
        double bondingContributions(arma::mat molecule_coordinates, bool verbose=false);

        /**
         * @brief Sums all the angle contributions in the molecule.
         * 
         * @param molecule_coordinates (arma::mat): Coordinates of the molecule.
         * @param verbose (bool): Whether to output verbose information.
         * @return double 
         */
        double angleContributions(arma::mat molecule_coordinates, bool verbose = false);

        /**
         * @brief Sums all the stretch bending contributions in the molecule.
         * 
         * @param molecule_coordinates (arma::mat): Coordinates of the molecule.
         * @return double 
         */
        double angleStretchBendingContributions(arma::mat molecule_coordinates);
        
        // Saturated aliphatics do not have oop contributions, so this is not implemented for now.
        void oopContributions();

        /**
         * @brief Sums all the torsional contributions in the molecule.
         * 
         * @param molecule_coordinates (arma::mat): Coordinates of the molecule.
         * @return double 
         */
        double torsionalContributions(arma::mat molecule_coordinates);

        /**
         * @brief Sums all the vanderwaals contributions in the molecule.
         * 
         * @param molecule_coordinates (arma::mat): Coordinates of the molecule.
         * @return double 
         */
        double vdwContributions(arma::mat molecule_coordinates);

        // Saturated aliphatics do not have electrostatic contributions, so this is not implemented for now.
        void electrostaticContributions();

        /**
         * @brief Sums all the energy contributions in the molecule.
         * 
         * @param molecule_coordinates (arma::mat): Coordinates of the molecule.
         * @param verbose (bool): Whether to output verbose information.
         * @return double 
         */
        double total(arma::mat molecule_coordinates, bool verbose=false);

        /**
         * @brief Sums all the energy contributions in the molecule.
         * 
         * @param coordinates (arma::mat): Coordinates of the molecule.
         * @param name (string): Name of the output file. The file type will be sdf.
         */
        void export_sdf(arma::mat coordinates, std::string name);
};