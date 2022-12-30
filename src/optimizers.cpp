#include "optimizers.h"

// Function to calculate central difference using input matrix
arma::mat grad(arma::mat Coordinates, double step_size, int num_atoms, calcEnergy& Molecule)
{
    arma::mat zero;
    zero.zeros(num_atoms,3);

    // Loop through all atoms.
    for (int i=0; i < num_atoms; i++)
    {
        // Loop through x y and z coordinates to get the gradient at each one keeping the rest constant.
        for (int k=0; k < 3; k++)
        {
            // Get current atom coordinates (x, y, z).
            arma::mat new_coordinates_plus = Coordinates;
            arma::mat new_coordinates_minus = Coordinates;

            // Step size adjustments for x, y or z.
            new_coordinates_plus(i, k) += step_size;
            new_coordinates_minus(i, k) -= step_size;

            // Calculate the plus and minus energies based on the movement of the atoms
            double plus_E = Molecule.total(new_coordinates_plus);
            double minus_E = Molecule.total(new_coordinates_minus);

            double force = (plus_E - minus_E) / (2*step_size);
            zero(i, k) = force;
        }
    }
    return zero;
}

// Steepest Descent
// Trivial, will not work for all starting arrangements
arma::mat steepest_descent(arma::mat Coordinates, double step_size, double tol, int num_atoms, calcEnergy& Molecule)
{
    arma::mat initial_guess = Coordinates;

    initial_guess.print("Initial Guess");

    arma::mat derivative = grad(initial_guess, step_size, num_atoms, Molecule);
    double initial_energy = Molecule.total(initial_guess);
    std::cout << "Initial Energy: " << initial_energy << std::endl;

    int count = 0;

    
    while (arma::norm(derivative, 2) > tol)
    {
        arma::mat new_positions = initial_guess - (derivative / arma::norm(derivative, 2)) * step_size;
    
        if (Molecule.total(new_positions) < Molecule.total(initial_guess))
        {
            initial_guess = new_positions;
            step_size *= 1.2;
            derivative = grad(initial_guess, step_size, num_atoms, Molecule);

            std::cout << std::endl;
            std::cout << "Step: " << count << std::endl;
            std::cout << "Energy: " << Molecule.total(initial_guess) << std::endl;
            std::cout << "Force Norm: " << arma::norm(derivative, 2) << std::endl;
            initial_guess.print("New positions");
            std::cout << "New step size: " << step_size << std::endl;
            std::cout << std::endl;

        }

        else
        {
            step_size *= 0.5;
        }

    count += 1;
    }
return initial_guess;
}

// Backtrack line-search with wolfe conditions
double line_search(int num_atoms, arma::mat coordinates, arma::mat search_direction, arma::mat gradient, calcEnergy& Molecule)
{
    // Coordinates and Gradient are in the size of (NUM_ATOMS X 3)
    double a = 1;
    double c1 = 1e-4;
    double c2 = 0.9;
    double initial_energy = Molecule.total(coordinates);

    arma::mat new_coordinates = coordinates + a * search_direction;
    arma::mat new_grad = grad(new_coordinates, 0.0001, num_atoms, Molecule);

    // Use the frobenius inner product between gradient and search direction
    while ( (Molecule.total(new_coordinates) >= (initial_energy + (c1*a*arma::dot(gradient,search_direction) )) || (arma::dot(new_grad, search_direction) <= (c2*arma::dot(gradient, search_direction)) )) )
    {
        a *= 0.5;
        new_coordinates = coordinates + a * search_direction; 
        new_grad = grad(new_coordinates, 0.0001, num_atoms, Molecule);
    }

    return a;
}

// BFGS
arma::mat BFGS(arma::mat Coordinates, double tol, int num_atoms, calcEnergy& Molecule)
{
    // Initial Guess
    arma::mat initial_guess = Coordinates;
    
    // Initial Gradient
    std::cout << "Inside BFGS" << std::endl;
    arma::mat initial_gradient = grad(initial_guess, 0.0001, num_atoms, Molecule);
    initial_gradient.print("Initial Gradient");

    // Initial Hessian
    // Should be 3N x 3N because the optimization degrees of freedom is 3N.
    arma::mat H = arma::eye(3*num_atoms, 3*num_atoms);

    int step = 0;

    while (arma::norm(initial_gradient, 2) > tol)
    {
        std::cout << "Step: " << step << std::endl;
        step++;
        // Starting here needs to go inside the while loop
        // Convert initial gradient to a column vector to multiply with the initial identity hessian to yield a (15x1) search direction vector.
        arma::vec colvec_gradient = initial_gradient.as_row().t();
        arma::vec search_direction = -H*colvec_gradient;

        // Reshaping the (15x1) search direction into 5x3 to use in our line search gives us the correct orientation for each atom and coordinate.
        arma::mat search_direction_mat = arma::reshape(search_direction, 3, num_atoms).t();
        
        double a = line_search(num_atoms, initial_guess, search_direction_mat, initial_gradient, Molecule);
        
        arma::vec s = search_direction*a;
        arma::mat new_coordinates = initial_guess + a*search_direction_mat;
        
        // Checking energy
        std::cout << "Energy: " << Molecule.total(new_coordinates) << std::endl;

        arma::mat new_gradient = grad(new_coordinates, 0.0001, num_atoms, Molecule);

        arma::mat y = new_gradient - initial_gradient;
        arma::vec y_vector = y.as_row().t();

        double r = 1/arma::dot(y_vector, s);
        arma::mat li = (arma::eye(3*num_atoms, 3*num_atoms) - (r*((s*y_vector.t()))));
        arma::mat ri = (arma::eye(3*num_atoms, 3*num_atoms) - (r*((y_vector*s.t()))));

        arma::mat hess_inter = li*H*ri;

        H = hess_inter + (r*(s*s.t()));

        // Update initial gradient to new gradient and initial guess to new coordinates.
        initial_gradient = new_gradient;
        initial_guess = new_coordinates;
    }
    initial_guess.print("Optimized coordinates");
    
    return initial_guess;
}
