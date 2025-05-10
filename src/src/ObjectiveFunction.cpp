#include "ObjectiveFunction.h"

#include "cstdio"

ExampleFunction::ExampleFunction(Placement &placement) : BaseFunction(1), placement_(placement)
{
    printf("Fetch the information you need from placement database.\n");
    printf("For example:\n");
    printf("    Placement boundary: (%.f,%.f)-(%.f,%.f)\n", placement_.boundryLeft(), placement_.boundryBottom(),
           placement_.boundryRight(), placement_.boundryTop());
}

const double &ExampleFunction::operator()(const std::vector<Point2<double>> &input)
{
    // Compute the value of the function
    value_ = 3. * input[0].x * input[0].x + 2. * input[0].x * input[0].y +
             2. * input[0].y * input[0].y + 7.;
    input_ = input;
    return value_;
}

const std::vector<Point2<double>> &ExampleFunction::Backward()
{
    // Compute the gradient of the function
    grad_[0].x = 6. * input_[0].x + 2. * input_[0].y;
    grad_[0].y = 2. * input_[0].x + 4. * input_[0].y;
    return grad_;
}


Wirelength::Wirelength(Placement &placement, double gamma)
    : BaseFunction(placement.numModules()), placement_(placement), gamma_(gamma) {}



const double &Wirelength::operator()(const std::vector<Point2<double>> &input) {
    value_ = 0.0; // the current value of the function, initialize to 0
    input_ = input; // cache the input for backward pass

    for (size_t netId = 0; netId < placement_.numNets(); ++netId) {
        Net &net = placement_.net(netId);
        const size_t pinCount = net.numPins();
        if (pinCount == 0) continue;

        std::vector<double> x, y;
        for (size_t k = 0; k < pinCount; ++k) {
            Pin &pin = net.pin(k);
            int moduleId = pin.moduleId();
            if (placement_.module(moduleId).isFixed()) {
                x.push_back(pin.x());
                y.push_back(pin.y());
            } else {
                x.push_back(input[moduleId].x + (pin.x() - placement_.module(moduleId).centerX()));
                y.push_back(input[moduleId].y + (pin.y() - placement_.module(moduleId).centerY()));
            }
        }

        auto wa = [&](const std::vector<double> &coord, double sign) {
            double max_coord = *std::max_element(coord.begin(), coord.end());
            double sum_exp = 0.0, sum_pos = 0.0;
            for (auto c : coord) {
                double e = std::exp(sign * (c - max_coord) / gamma_);
                sum_exp += e;
                sum_pos += c * e;
            }
            return sum_pos / sum_exp;
        };

        double wx = wa(x, 1.0) - wa(x, -1.0);
        double wy = wa(y, 1.0) - wa(y, -1.0);
        value_ += wx + wy;
    }

    return value_;
}


const std::vector<Point2<double>> &Wirelength::Backward() {
    // Reset gradient vector
    for (auto &g : grad_) {
        g = Point2<double>(0.0, 0.0);
    }

    for (size_t netId = 0; netId < placement_.numNets(); ++netId) {
        Net &net = placement_.net(netId);
        size_t pinCount = net.numPins();
        if (pinCount == 0) continue;

        std::vector<double> x(pinCount), y(pinCount);
        std::vector<int> moduleIds(pinCount);
        std::vector<bool> isFixed(pinCount);

        // Step 1: Collect pin positions (as in operator())
        for (size_t k = 0; k < pinCount; ++k) {
            Pin &pin = net.pin(k);
            int moduleId = pin.moduleId();
            moduleIds[k] = moduleId;

            if (placement_.module(moduleId).isFixed()) {
                x[k] = pin.x();
                y[k] = pin.y();
                isFixed[k] = true;
            } else {
                isFixed[k] = false;
                double dx = pin.x() - placement_.module(moduleId).centerX();
                double dy = pin.y() - placement_.module(moduleId).centerY();
                x[k] = input_[moduleId].x + dx;
                y[k] = input_[moduleId].y + dy;
            }
        }

        auto computeGrad = [&](const std::vector<double> &coord, bool isX) {
            double max_coord = *std::max_element(coord.begin(), coord.end());
            double min_coord = *std::min_element(coord.begin(), coord.end());

            // WA max
            std::vector<double> emax(pinCount), emax_sum_each(pinCount);
            double sum_emax = 0, sum_x_emax = 0;
            for (size_t i = 0; i < pinCount; ++i) {
                emax[i] = std::exp((coord[i] - max_coord) / gamma_);
                sum_emax += emax[i];
                sum_x_emax += coord[i] * emax[i];
            }
            double wa_max = sum_x_emax / sum_emax;

            // WA min
            std::vector<double> emin(pinCount), emin_sum_each(pinCount);
            double sum_emin = 0, sum_x_emin = 0;
            for (size_t i = 0; i < pinCount; ++i) {
                emin[i] = std::exp(-(coord[i] - min_coord) / gamma_);
                sum_emin += emin[i];
                sum_x_emin += coord[i] * emin[i];
            }
            double wa_min = sum_x_emin / sum_emin;

            // Now compute derivative contribution for each pin
            for (size_t i = 0; i < pinCount; ++i) {
                int moduleId = moduleIds[i];
                if (isFixed[i]) continue;

                // Derivative of WA max
                double d_wa_max = emax[i] / sum_emax * (1 + (coord[i] - wa_max) / gamma_);
                // Derivative of WA min
                double d_wa_min = -emin[i] / sum_emin * (1 + (wa_min - coord[i]) / gamma_);
                double grad_val = d_wa_max + d_wa_min;

                if (isX)
                    grad_[moduleId].x += grad_val;
                else
                    grad_[moduleId].y += grad_val;
            }
        };

        // Step 2: Compute gradients
        computeGrad(x, true);  // x-direction
        computeGrad(y, false); // y-direction
    }

    return grad_;
}


Density::Density(Placement &placement, int bin_rows, int bin_cols, double sigma_factor, double target_density)
    : BaseFunction(placement.numModules()), placement_(placement),
      bin_rows_(bin_rows), bin_cols_(bin_cols), target_density_(target_density)
{
    chip_left_ = placement.boundryLeft();
    chip_right_ = placement.boundryRight();
    chip_bottom_ = placement.boundryBottom();
    chip_top_ = placement.boundryTop();

    bin_width_ = (chip_right_ - chip_left_) / bin_cols_;
    bin_height_ = (chip_top_ - chip_bottom_) / bin_rows_;
    sigma_ = sigma_factor * std::min(bin_width_, bin_height_);
    sigma_sq_ = sigma_ * sigma_;

    bin_capacity_ = bin_width_ * bin_height_ * target_density_;

    // Initialize bin density grid
    bin_density_.resize(bin_rows_, std::vector<double>(bin_cols_, 0.0));
}


const double &Density::operator()(const std::vector<Point2<double>> &input) {
    value_ = 0.0;
    input_ = input;  // Cache for use in Backward()

    // Reset bin density grid
    for (int i = 0; i < bin_rows_; ++i)
        std::fill(bin_density_[i].begin(), bin_density_[i].end(), 0.0);

    const int num_modules = placement_.numModules();

    for (int i = 0; i < num_modules; ++i) {
        Module &mod = placement_.module(i);
        if (mod.isFixed()) continue;

        const double module_area = mod.area();
        const double cx = input[i].x;  // Center x from optimizer
        const double cy = input[i].y;  // Center y

        // Determine Gaussian influence range (3σ around center)
        double x_min = cx - 3 * sigma_;
        double x_max = cx + 3 * sigma_;
        double y_min = cy - 3 * sigma_;
        double y_max = cy + 3 * sigma_;

        // Convert to bin index range (clamped to grid)
        int bin_x_min = std::max(0, (int)((x_min - chip_left_) / bin_width_));
        int bin_x_max = std::min(bin_cols_ - 1, (int)((x_max - chip_left_) / bin_width_));
        int bin_y_min = std::max(0, (int)((y_min - chip_bottom_) / bin_height_));
        int bin_y_max = std::min(bin_rows_ - 1, (int)((y_max - chip_bottom_) / bin_height_));

        // Apply Gaussian to all bins in range
        for (int by = bin_y_min; by <= bin_y_max; ++by) {
            double bin_center_y = chip_bottom_ + (by + 0.5) * bin_height_;
            double dy = bin_center_y - cy;

            for (int bx = bin_x_min; bx <= bin_x_max; ++bx) {
                double bin_center_x = chip_left_ + (bx + 0.5) * bin_width_;
                double dx = bin_center_x - cx;

                double dist_sq = dx * dx + dy * dy;
                double g = std::exp(-dist_sq / (2 * sigma_sq_)) / (2 * M_PI * sigma_sq_);
                bin_density_[by][bx] += module_area * g;
            }
        }
    }

    // Compute penalty over all bins
    for (int by = 0; by < bin_rows_; ++by) {
        for (int bx = 0; bx < bin_cols_; ++bx) {
            double overflow = bin_density_[by][bx] - bin_capacity_;
            if (overflow > 0.0)
                value_ += overflow * overflow;  // Quadratic penalty
        }
    }

    return value_;
}



const std::vector<Point2<double>> &Density::Backward() {
    const size_t num_modules = placement_.numModules();

    // Reset gradients
    grad_.assign(num_modules, Point2<double>(0.0, 0.0));

    for (size_t i = 0; i < num_modules; ++i) {
        Module &mod = placement_.module(i);
        if (mod.isFixed()) continue;

        const double area = mod.area();
        const double cx = input_[i].x;
        const double cy = input_[i].y;

        // Influence range = 3σ
        double x_min = cx - 3 * sigma_;
        double x_max = cx + 3 * sigma_;
        double y_min = cy - 3 * sigma_;
        double y_max = cy + 3 * sigma_;

        int bin_x_min = std::max(0, (int)((x_min - chip_left_) / bin_width_));
        int bin_x_max = std::min(bin_cols_ - 1, (int)((x_max - chip_left_) / bin_width_));
        int bin_y_min = std::max(0, (int)((y_min - chip_bottom_) / bin_height_));
        int bin_y_max = std::min(bin_rows_ - 1, (int)((y_max - chip_bottom_) / bin_height_));

        for (int by = bin_y_min; by <= bin_y_max; ++by) {
            double bin_center_y = chip_bottom_ + (by + 0.5) * bin_height_;
            double dy = bin_center_y - cy;

            for (int bx = bin_x_min; bx <= bin_x_max; ++bx) {
                double bin_center_x = chip_left_ + (bx + 0.5) * bin_width_;
                double dx = bin_center_x - cx;

                double dist_sq = dx * dx + dy * dy;
                double g = std::exp(-dist_sq / (2 * sigma_sq_)) / (2 * M_PI * sigma_sq_);
                double dD_dx = area * g * (dx / sigma_sq_);
                double dD_dy = area * g * (dy / sigma_sq_);

                double overflow = bin_density_[by][bx] - bin_capacity_;
                if (overflow > 0.0) {
                    grad_[i].x += 2.0 * overflow * (-dD_dx);  // negative sign because bin moves opposite to module
                    grad_[i].y += 2.0 * overflow * (-dD_dy);
                }
            }
        }
    }

    return grad_;
}




// ObjectiveFunction::ObjectiveFunction(Placement &placement, double lambda)
//     : BaseFunction(placement.numModules()),
//       wirelength_(placement, /*gamma=*/5000.0),  // set γ as needed
//       density_(placement),                       // default: 50×50 grid
//       lambda_(lambda),
//       grad_(placement.numModules(), Point2<double>(0.0, 0.0)) {}

// const double &ObjectiveFunction::operator()(const std::vector<Point2<double>> &input) {
//     const double wl = wirelength_(input);
//     const double dp = density_(input);
//     value_ = wl + lambda_ * dp;
//     return value_;
// }


// const std::vector<Point2<double>> &ObjectiveFunction::Backward() {
//     const std::vector<Point2<double>> &grad_wl = wirelength_.Backward();
//     const std::vector<Point2<double>> &grad_dp = density_.Backward();

//     for (size_t i = 0; i < grad_.size(); ++i) {
//         grad_[i].x = grad_wl[i].x + lambda_ * grad_dp[i].x;
//         grad_[i].y = grad_wl[i].y + lambda_ * grad_dp[i].y;
//     }

//     return grad_;
// }


ObjectiveFunction::ObjectiveFunction(Placement &placement, double lambda)
    : BaseFunction(placement.numModules()),
      wirelength_(placement, /*gamma=*/5000.0),  // unused in dummy
      density_(placement),                       // unused in dummy
      lambda_(lambda),
      grad_(placement.numModules(), Point2<double>(0.0, 0.0)) {}

const double &ObjectiveFunction::operator()(const std::vector<Point2<double>> &input) {
    input_ = input;  // ← Add this at the start of operator()

    value_ = 0.0;
    for (const auto &p : input) {
        value_ += (p.x * p.x + p.y * p.y);            // dummy "wirelength"
        value_ += lambda_ * (p.x + p.y);              // dummy "density penalty"
    }
    value_ = 1.0;  // dummy constant offset
    return value_;
}

const std::vector<Point2<double>> &ObjectiveFunction::Backward() {
    for (size_t i = 0; i < grad_.size(); ++i) {
        grad_[i].x = 1.5;  // dummy gradient
        grad_[i].y = 1;  // dummy gradient
    }
    return grad_;
}



void ObjectiveFunction::setLambda(double lambda) {
    lambda_ = lambda;
}

double ObjectiveFunction::getLambda() const {
    return lambda_;
}
