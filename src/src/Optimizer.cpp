#include "Optimizer.h"
#include <iostream>
#include <cmath>
#include <fstream>   // for std::ofstream
#include <vector>    // for std::vector
#include <string>    // for std::string (if you use file paths with std::string)

SimpleConjugateGradient::SimpleConjugateGradient(BaseFunction &obj,
                                                 std::vector<Point2<double>> &var,
                                                 const double &alpha,
                                                 double boundary_left,
                                                 double boundary_right,
                                                 double boundary_top,
                                                 double boundary_bottom)
        : BaseOptimizer(obj, var),
        grad_prev_(var.size()),
        dir_prev_(var.size()),
        step_(0),
        alpha_(alpha) {
        boundary_left_ = boundary_left;
        boundary_right_ = boundary_right;
        boundary_top_ = boundary_top;
        boundary_bottom_ = boundary_bottom;
    }

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

    std::ofstream grad_file("plot_output/grad_vectors.txt");

    for (size_t i = 0; i < var_.size(); ++i) {
        const auto &pos = var_[i];
        const auto &g = obj_.grad()[i];
        double mag = std::sqrt(g.x * g.x + g.y * g.y);
        grad_file << pos.x << " " << pos.y << " " << -g.x << " " << -g.y << " " << mag << "\n";
    }


    grad_file.close();
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

        }
    }

    // Assume the step size is constant
    // TODO(Optional): Change to dynamic step-size control

    // === Dynamic Step-Size ===
    // chip width 60000 -> 120000 step size
    // chip width 2300 -> 4000 step size
    // ratio = 2
    // double s = 1000000.0; // target average move distance (tunable)
    // cout << "boundary_left_ = " << boundary_left_ << ", boundary_right_ = " << boundary_right_ << endl;
    double s = ((boundary_right_ - boundary_left_) * 2.3); // target average move distance (tunable)
    double norm = 0.0;
    for (size_t i = 0; i < kNumModule; ++i) {
        norm += dir[i].x * dir[i].x + dir[i].y * dir[i].y;
    }
    norm = std::sqrt(norm); // L2 norm

    double dynamic_alpha = (norm < 1e-12) ? 0.0 : s / norm;
    setAlpha(dynamic_alpha);  // or alpha_ = dynamic_alpha;

    // std::cout << "[Dynamic Step Size] alpha = " << dynamic_alpha << ", norm = " << norm << std::endl;

    // Update the solution
    // Please be aware of the updating directions, i.e., the sign for each term.
    for (size_t i = 0; i < kNumModule; ++i) {
        // Point2<double> ori = var_[i]; 
        var_[i] = var_[i] + alpha_ * dir[i]; 
        // Point2<double> new_pos = var_[i]; 
 
        // var_[i].x = max(boundary_left_ ,  min(var_[i].x, boundary_right_));
        // var_[i].y = max(boundary_bottom_, min(var_[i].y, boundary_top_));
        // if(var_[i].x < boundary_left_)
        //     var_[i].x = boundary_left_;
        // else if(var_[i].x > boundary_right_)
        //     var_[i].x = boundary_right_;
        // if(var_[i].y < boundary_bottom_)
        //     var_[i].y = boundary_bottom_;
        // else if(var_[i].y > boundary_top_)
        //     var_[i].y = boundary_top_;


    }

    // Update the cache data members
    grad_prev_ = obj_.grad();
    dir_prev_ = dir;
    step_++;
}
