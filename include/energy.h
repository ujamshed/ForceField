#pragma once
#include "input.h"
#include "utility.h"
#include "output.h"

class calcEnergy
{
    private:
        Molecule* _Mol;

        // Maps atom types to bond stretching parameters
        // First pair is the atom types involved in the bond, and the second pair are the params
        // First element of the second pair is the kb, while the second element of the second pair is the r0
        // Need 1,5 and 5,1 because C-H and H-C when considering the stretch bending contributions.
        std::map<std::pair<int, int>, std::pair<double, double>> bs_params = { {{1,1} , {4.258, 1.508}}, {{1,5} , {4.766, 1.093}}, {{5,1} , {4.766, 1.093}} };
        
        // Maps atom types to angle parameters
        // First tuple is the atom types involved in the angle, and the second pair are the params
        // First element of the second pair is the ka, while the second element of the second pair is the a0
        std::map<std::tuple<int, int, int>, std::pair<double, double>> ab_params = { {{1,1,5} , {0.636, 110.549}}, {{5, 1, 5}, {0.516, 108.836}} };
        
        // Maps atom types to stretch bond parameters
        // First tuple is the atom types involved in the interaction and second pair are the params
        // First element of the second pair is the kbaIJK, while the second element of the second pair is the kbaKJI
        std::map<std::tuple<int, int, int>, std::pair<double, double>> asb_params = { {{1,1,5} , {0.227, 0.070}}, {{5, 1, 5}, {0.115, 0.115}} };

        // Maps atom types to torsion parameters
        // First tuple is the atom types involved in the interaction and the second pair are the params
        // Elements of the second tuple are V1, V2, V3 respectively
        std::map<std::tuple<int, int, int, int>, std::tuple<double, double, double>> t_params = { {{5,1,1,5} , {0.284, -1.386, 0.314}} };

        // Maps atom types to vdw parameters
        // Int is the atom type involved in the interaction
        // Elements of the second tuple are alpha, N, A and G, respectively
        std::map<int, std::tuple<double, double, double, double>> vdw_params = { {1, {1.050, 2.490, 3.890, 1.282}}, {5, {0.250, 0.800, 4.200, 1.209}} };


    public:
        
        calcEnergy(Molecule* Mol);

        // Need helper function to calculate the distance between the two atoms i and j.
        // Calculates the bond stretching energy between 2 atoms, i and j.
        double bondStretching(double kb, double delta_R);

        // delta_Angle is angle_ijk - reference angle.
        // Calculates the angle bending between atoms i, j and k, where j is the middle atom.
        double angleBending(double ka, double delta_Angle);

        // Calculates the angle stretch bending between atoms i, j, k. Distances between i, j and k, j need to be calculated, as well
        // as the angle between i, j, k.
        double angleStretchBending(double kba_ijk, double kba_kji, double delta_Rij, double delta_Rkj, double delta_Angle);

        // Need helper function to calculate the wilson angle (X) between the bond j-l, and the plane i-j-k.
        // Calculates the out of plane bending at atom j in i,j,k.
        double outOfPlaneBending(double koop, double X);

        // Calculates the torsional energy based on atoms i, j, k, l.
        double torsion(double V1, double V2, double V3, double omega);

        // Calculates the vanderwaals interaction based on the buffered-14-7 form, based on atoms i and j.
        double vdw(double dist_Rij, double Ai, double Aj, double alphai, double alphaj, double Ni, double Nj, double Gi, double Gj);

        // Need helper function to calculate the partial atomic charges (qi, qj) based on bond charge increments (bci)
        // Calculates the electrostatic interaction between atoms i, j.
        double electrostatic(double delta_R, double qi, double qj);

        double bondingContributions(arma::mat molecule_coordinates, bool verbose=false);

        double angleContributions(arma::mat molecule_coordinates, bool verbose = false);

        double angleStretchBendingContributions(arma::mat molecule_coordinates);
        
        // Saturated aliphatics do not have oop contributions, so this is not implemented for now.
        void oopContributions();

        double torsionalContributions(arma::mat molecule_coordinates);

        double vdwContributions(arma::mat molecule_coordinates);

        // Saturated aliphatics do not have electrostatic contributions, so this is not implemented for now.
        void electrostaticContributions();

        double total(arma::mat molecule_coordinates, bool verbose=false);

        // Function to export the molecule information in sdf format
        void export_sdf(arma::mat coordinates, std::string name);
};