#define _GLIBCXX_USE_CXX11_ABI 0  // Align the ABI version to avoid compatibility issues with `Placment.h`
#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include <vector>

#include "Placement.h"
#include "Point.h"

/**
 * @brief Base class for objective functions
 */
class BaseFunction {
   public:
    /////////////////////////////////
    // Conssutructors
    /////////////////////////////////

    BaseFunction(const size_t &input_size) : grad_(input_size) {}

    /////////////////////////////////
    // Accessors
    /////////////////////////////////

    const std::vector<Point2<double>> &grad() const { return grad_; }
    const double &value() const { return value_; }

    /////////////////////////////////
    // Methods
    /////////////////////////////////

    // Forward pass, compute the value of the function
    virtual const double &operator()(const std::vector<Point2<double>> &input) = 0;

    // Backward pass, compute the gradient of the function
    virtual const std::vector<Point2<double>> &Backward() = 0;

   protected:
    /////////////////////////////////
    // Data members
    /////////////////////////////////

    std::vector<Point2<double>> grad_;  // Gradient of the function
    double value_;                      // Value of the function
};

/**
 * @brief Example function for optimization
 *
 * This is a simple example function for optimization. The function is defined as:
 *      f(t) = 3*t.x^2 + 2*t.x*t.y + 2*t.y^2 + 7
 */
class ExampleFunction : public BaseFunction {
   public:
    /////////////////////////////////
    // Constructors
    /////////////////////////////////

    ExampleFunction(Placement &placement);

    /////////////////////////////////
    // Methods
    /////////////////////////////////

    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;

   private:
    /////////////////////////////////
    // Data members
    /////////////////////////////////

    std::vector<Point2<double>> input_;  // Cache the input for backward pass
    Placement &placement_;
};

/**
 * @brief Wirelength function
 */

class Wirelength : public BaseFunction {
    public:
        Wirelength(Placement &placement, double gamma);

        const double &operator()(const std::vector<Point2<double>> &input) override;
        const std::vector<Point2<double>> &Backward() override;

    private:
        Placement &placement_;
        double gamma_;
        std::vector<Point2<double>> input_;  // Cache input
};


/**
 * @brief Density function
 */


// class Density : public BaseFunction {
//     public:
//         Density(Placement &placement, int bin_rows = 50, int bin_cols = 50, double sigma_factor = 1.5, double target_density = 0.9);

//         const double &operator()(const std::vector<Point2<double>> &input) override;
//         const std::vector<Point2<double>> &Backward() override;
//         const double getBinCapacity() const { return bin_capacity_; }
//         const std::vector<std::vector<double>> &getBinDensity() const { return bin_density_; }
//     private:
//         Placement &placement_;

//         int bin_rows_, bin_cols_;
//         double chip_left_, chip_right_, chip_top_, chip_bottom_;
//         double bin_width_, bin_height_;
//         double sigma_, sigma_sq_; // Gaussian smoothing σ and σ²
//         double target_density_;
//         double bin_capacity_;

//         std::vector<std::vector<double>> bin_density_; // Smoothed density per bin
//         std::vector<Point2<double>> input_;            // Cached module positions
// };

class Density : public BaseFunction {
    public:
        Density(Placement &placement, int bin_rows = 50, int bin_cols = 50, double alpha = 5.0, double target_density = 0.9);


        
        const double &operator()(const std::vector<Point2<double>> &input) override;
        const std::vector<Point2<double>> &Backward() override;

        const double getBinCapacity() const { return bin_capacity_; }
        const std::vector<std::vector<double>> &getBinDensity() const { return bin_density_; }

        // Optional: expose smoothing trigger
        void smoothBinDensityLevel(int smoothing_pass = 1);

    private:
        Placement &placement_;

        int bin_rows_, bin_cols_;
        double chip_left_, chip_right_, chip_top_, chip_bottom_;
        double bin_width_, bin_height_;
        double alpha_;                  // Alpha for sigmoid steepness
        double target_density_;
        double bin_capacity_;

        std::vector<std::vector<double>> bin_density_; // Smoothed density per bin
        std::vector<Point2<double>> input_;            // Cached module positions

        // Sigmoid function used for smoothing density influence
        double sigmoid(double d, double lower, double upper) const;
        double sigmoid_derivative(double d, double lowwer, double upper) const;

};




/**
 * @brief Objective function for global placement
 */


class ObjectiveFunction : public BaseFunction {
    public:
        ObjectiveFunction(Placement &placement, double lambda, Wirelength wirelength, Density density);

        const double &operator()(const std::vector<Point2<double>> &input) override;
        const std::vector<Point2<double>> &Backward() override;

        void setLambda(double lambda);  // Optional: expose dynamic λ adjustment
        double getLambda() const;
        const Wirelength &getWirelength() const { return wirelength_; }
        const Density &getDensity() const { return density_; }

    private:
        Wirelength wirelength_;
        Density density_;
        double lambda_;
        // std::vector<Point2<double>> grad_;  // Combined gradient cache
        std::vector<Point2<double>> input_;            // Cached module positions
    
};




#endif  // OBJECTIVEFUNCTION_H
