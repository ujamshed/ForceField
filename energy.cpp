#include <cmath>

class calcEnergy
{
    private:

    public:
        
        // Empty constructor for now
        calcEnergy()=default;

        // Need helper function to calculate the distance between the 
        // Calculates the bond stretching energy between 2 atoms, i and j.
        double bondStretching(double kb, double delta_R)
        {
            // cs is the cubic strech constant at -2 A^-1
            double cs = -2;
            return 143.9325 * (kb / 2) * delta_R * delta_R * (1 + cs*delta_R + 7/12*cs*cs*delta_R*delta_R);
        }

        // Need helper function to calculate angle between the 3 atoms
        // delta_Angle is angle_ijk - reference angle.
        // Calculates the angle bending between atoms i, j and k, where j is the middle atom.
        double angleBending(double ka, double delta_Angle)
        {
            // cb is the cubic bend constant at -0.007 deg^-1 or -0.4 rad^-1
            double cb = -0.007;
            return 0.043844 * (ka/2) * delta_Angle * delta_Angle * (1 + cb* delta_Angle);
        }

        // Calculates the angle stretch bending between atoms i, j, k. Distances between i, j and k, j need to be calculated, as well
        // as the angle between i, j, k.
        double angleStretchBending(double kba_ijk, double kba_kji, double delta_Rij, double delta_Rkj, double delta_Angle)
        {
            return 2.51210*(kba_ijk*delta_Rij + kba_kji*delta_Rkj)*delta_Angle;
        }

        // Need helper function to calculate the wilson angle (X) between the bond j-l, and the plane i-j-k.
        // Calculates the out of plane bending at atom j in i,j,k.
        double outOfPlaneBending(double koop, double X)
        {
            return 0.043844 * (koop / 2) * X * X;
        }

        // Need helper function to calculate the torsional angle (omega) between atoms i, j, k, l.
        // Calculates the torsional energy based on atoms i, j, k, l.
        double torsion(double V1, double V2, double V3, double omega)
        {
            return 0.5* (V1 * (1 + cos(omega)) + V2 * (1 - cos(2*omega)) + V3 * (1 + cos(3*omega)));
        }

        double vdw()
        {
            
        }

        double electrostatic()
        {

        }
};

    