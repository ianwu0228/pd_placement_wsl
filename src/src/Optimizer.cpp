#include "Optimizer.h"

#include <cmath>

SimpleConjugateGradient::SimpleConjugateGradient(BaseFunction &obj,
                                                 std::vector<Point2<double>> &var,
                                                 const double &alpha)
    : BaseOptimizer(obj, var),
      grad_prev_(var.size()),
      dir_prev_(var.size()),
      step_(0),
      alpha_(alpha) {}

void SimpleConjugateGradient::Initialize() {
    // Before the optimization starts, we need to initialize the optimizer.
    step_ = 0;
}

/**
 * @details Update the solution once using the conjugate gradient method.
 */
void SimpleConjugateGradient::Step() {
    const size_t &kNumModule = var_.size();

    // Compute the gradient direction
    obj_(var_);       // Forward, compute the function value and cache from the input
    obj_.Backward();  // Backward, compute the gradient according to the cache
    // cout << "obj value: " << obj_.value() << endl; 

    // Compute the Polak-Ribiere coefficient and conjugate directions
    double beta;                                  // Polak-Ribiere coefficient
    std::vector<Point2<double>> dir(kNumModule);  // conjugate directions
    if (step_ == 0) {
        // For the first step, we will set beta = 0 and d_0 = -g_0
        beta = 0.;
        for (size_t i = 0; i < kNumModule; ++i) {
            dir[i] = -obj_.grad().at(i);
        }
    } else {
        // For the remaining steps, we will calculate the Polak-Ribiere coefficient and
        // conjugate directions normally
        double t1 = 0.;  // Store the numerator of beta
        double t2 = 0.;  // Store the denominator of beta
        for (size_t i = 0; i < kNumModule; ++i) {
            Point2<double> t3 =
                obj_.grad().at(i) * (obj_.grad().at(i) - grad_prev_.at(i));
            t1 += t3.x + t3.y;
            t2 += std::abs(obj_.grad().at(i).x) + std::abs(obj_.grad().at(i).y);
        }
        // beta = t1 / (t2 * t2);
        const double epsilon = 1e-10;
        if (std::abs(t2) < epsilon) {
            beta = 0.0;
        } else {
            beta = t1 / (t2 * t2);
        }

        for (size_t i = 0; i < kNumModule; ++i) {
            dir[i] = -obj_.grad().at(i) + beta * dir_prev_.at(i);
            // cout << "onj_.grad().at(i) (x, y) = " << obj_.grad().at(i).x <<  obj_.grad().at(i).y << endl;
            // cout << "dir[" << i << "], dir.x =  " << dir[i].x << "dir.y = " << dir[i].y << endl; 
        }
    }

    // Assume the step size is constant
    // TODO(Optional): Change to dynamic step-size control

    // Update the solution
    // Please be aware of the updating directions, i.e., the sign for each term.
    for (size_t i = 0; i < kNumModule; ++i) {
        Point2<double> ori = var_[i]; // Store the current position
        var_[i] = var_[i] + alpha_ * dir[i]; //supports operator overloading, updates x and y at the same time
        Point2<double> new_pos = var_[i]; // Store the new position
        // cout << "dir[i]: " << dir[i].x << ", " << dir[i].y << endl;
        // cout << "Module " << i << ": (" << ori.x << ", " << ori.y << ") -> ("
        //      << new_pos.x << ", " << new_pos.y << ")" << endl;
    }

    // Update the cache data members
    grad_prev_ = obj_.grad();
    dir_prev_ = dir;
    step_++;
}
