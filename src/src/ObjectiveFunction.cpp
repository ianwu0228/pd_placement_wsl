#include "ObjectiveFunction.h"
#include "cstdio"
using namespace std;

// example function
/*
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

*/




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

        auto wa = [&](const vector<double> &coord, double sign) {
            double max_coord = *max_element(coord.begin(), coord.end());
            double min_coord = *min_element(coord.begin(), coord.end());
            double sum_exp = 0.0, sum_pos = 0.0;
            double e = 0.0;
            for (auto c : coord) {
                if(sign == 1) 
                {
                    e = exp(sign * (c - max_coord) / gamma_);
                }
                else
                {
                    e = exp(sign * (c - min_coord) / gamma_);

                }
                
                sum_exp += e;
                sum_pos += c * e;
            }
            return sum_pos / sum_exp;
        };

        double wx = wa(x, 1.0) - wa(x, -1.0);
        double wy = wa(y, 1.0) - wa(y, -1.0);
        value_ += wx + wy;
        // cout << "Wirelength for net " << netId << ": " << wx + wy << endl;
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

        auto computeGrad = [&](const vector<double> &coord, bool isX) {
            double max_coord = *max_element(coord.begin(), coord.end());
            double min_coord = *min_element(coord.begin(), coord.end());

            // WA max
            vector<double> emax(pinCount), emax_sum_each(pinCount);
            double sum_emax = 0, sum_x_emax = 0;
            for (size_t i = 0; i < pinCount; ++i) {
                emax[i] = exp((coord[i] - max_coord) / gamma_);
                sum_emax += emax[i];
                sum_x_emax += coord[i] * emax[i];
            }
            double wa_max = sum_x_emax / sum_emax;

            // WA min
            vector<double> emin(pinCount), emin_sum_each(pinCount);
            double sum_emin = 0, sum_x_emin = 0;
            for (size_t i = 0; i < pinCount; ++i) {
                emin[i] = exp(-(coord[i] - min_coord) / gamma_);
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
    // cout << "Wirelength grad for module 0: " << grad_[0].x << ", " << grad_[0].y << std::endl;

    return grad_;
}


Density::Density(Placement &placement, int bin_rows, int bin_cols, double alpha, double target_density)
    : BaseFunction(placement.numModules()), placement_(placement),
      bin_rows_(bin_rows), bin_cols_(bin_cols), alpha_(alpha), target_density_(target_density)
{

    chip_left_ = placement.boundryLeft();
    chip_right_ = placement.boundryRight();
    chip_bottom_ = placement.boundryBottom();
    chip_top_ = placement.boundryTop();

    bin_width_ = (chip_right_ - chip_left_) / bin_cols_;
    bin_height_ = (chip_top_ - chip_bottom_) / bin_rows_;

    bin_capacity_ = bin_width_ * bin_height_ * target_density_;

    // Initialize bin density grid
    bin_density_.resize(bin_rows_, std::vector<double>(bin_cols_, 0.0));
}


double Density::sigmoid(double d, double lower, double upper) const {
    double f = 1/(1 + exp(-alpha_ * (d - lower)));
    f *= 1/(1 + exp(-alpha_ * (upper - d)));
    return f;
}


double Density::sigmoid_derivative(double d, double lower, double upper) const {
    // Debug prints
    if (std::abs(d) > 1e6 || std::abs(lower) > 1e6 || std::abs(upper) > 1e6) {
        cout << "Large input values: d=" << d << ", lower=" << lower << ", upper=" << upper << endl;
    }

    // Limit exponential arguments
    double exp1 = std::min(std::max(-alpha_ * (d - lower), -500.0), 500.0);
    double exp2 = std::min(std::max(-alpha_ * (upper - d), -500.0), 500.0);

    double term1_num = alpha_ * exp(exp1);
    double term1_den1 = 1 + exp(exp2);
    double term1_den2 = pow(1 + exp(exp1), 2);
    
    // Check for division by zero
    if (term1_den1 < 1e-300 || term1_den2 < 1e-300) {
        return 0.0;  // Return safe value
    }
    
    double term1 = term1_num / (term1_den1 * term1_den2);

    double term2_num = alpha_ * exp(exp2);
    double term2_den1 = 1 + exp(exp1);
    double term2_den2 = pow(1 + exp(exp2), 2);
    
    if (term2_den1 < 1e-300 || term2_den2 < 1e-300) {
        return 0.0;  // Return safe value
    }
    
    double term2 = term2_num / (term2_den1 * term2_den2);

    return term1 - term2;
}



double Density::bell_shaped_potential(double dx, double wv, double wb) {
    dx = std::abs(dx);
    double a = 4.0 / ((wv + 2 * wb) * (wv + 4 * wb));
    double b = 2.0 / (wb * (wv + 4 * wb));
    
    if (dx <= (wv / 2.0 + wb))
        return 1.0 - a * dx * dx;
    else if (dx <= wv / 2.0 + 2 * wb)
        return b * pow(dx - wv / 2.0 - 2 * wb, 2.0);
    else
        return 0.0;
}

vector<vector<double>> Density::generateGaussianKernel(int size, double sigma) {
    vector<vector<double>> kernel(size, vector<double>(size));
    double sum = 0.0;
    int center = size / 2;

    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j) {
            double x = i - center;
            double y = j - center;
            kernel[i][j] = exp(-(x*x + y*y) / (2 * sigma * sigma));
            sum += kernel[i][j];
        }

    for (auto& row : kernel)
        for (auto& val : row)
            val /= sum;

    return kernel;
}

void Density::applyGaussianSmoothing(std::vector<std::vector<double>> &density, int size, double sigma) {
    auto kernel = generateGaussianKernel(size, sigma);
    std::vector<std::vector<double>> smoothed = density;

    int offset = size / 2;
    for (int y = 0; y < bin_rows_; ++y) {
        for (int x = 0; x < bin_cols_; ++x) {
            double val = 0.0;
            for (int i = 0; i < size; ++i)
                for (int j = 0; j < size; ++j) {
                    int ny = y + i - offset;
                    int nx = x + j - offset;
                    if (ny >= 0 && ny < bin_rows_ && nx >= 0 && nx < bin_cols_)
                        val += density[ny][nx] * kernel[i][j];
                }
            smoothed[y][x] = val;
        }
    }
    density = smoothed;
}


double Density::bell_shaped_potential_derivative(double dx, double wv, double wb) {
    double abs_dx = abs(dx);
    double sign = (dx >= 0) ? 1.0 : -1.0;
    double a = 4.0 / ((wv + 2 * wb) * (wv + 4 * wb));
    double b = 2.0 / (wb * (wv + 4 * wb));

    if (abs_dx <= wv / 2.0 + wb)
        return -2 * a * dx;
    else if (abs_dx <= wv / 2.0 + 2 * wb)
        return 2 * b * (abs_dx - wv / 2.0 - 2 * wb) * sign;
    else
        return 0.0;
}



const double& Density::operator()(const std::vector<Point2<double>> &input) {

    value_ = 0.0;
    input_ = input;  
    // Reset bin density grid
    for (int i = 0; i < bin_rows_; ++i)
        std::fill(bin_density_[i].begin(), bin_density_[i].end(), 0.0);

    const int num_modules = placement_.numModules();
    for(int i = 0; i < num_modules; ++i)
    {
        Module &mod = placement_.module(i);
        if(mod.isFixed()) continue;

        // some constants
        const double mod_area = mod.area();
        const double mod_center_x = input[i].x;
        const double mod_center_y = input[i].y;
        const double mod_h = mod.height();
        const double mod_w = mod.width();

        const double influence_coefficient = 2;
        double influence_range_x = mod_w * influence_coefficient;
        double influence_range_y = mod_h * influence_coefficient;
        
        // set up the affecting region for the cell
        double x_min = mod_center_x - influence_range_x;
        double x_max = mod_center_x + influence_range_x;
        double y_min = mod_center_y - influence_range_y;
        double y_max = mod_center_y + influence_range_y;

        // map to specific grid(s)
        int bin_x_min = max(0, (int)((x_min - chip_left_) / bin_width_));
        int bin_x_max = min(bin_cols_ - 1, (int)((x_max - chip_left_) / (bin_width_)));
        int bin_y_min = max(0, (int)((y_min - chip_bottom_) / (bin_height_)));
        int bin_y_max = min(bin_rows_ - 1, (int)((y_max - chip_bottom_) / (bin_height_)));

        // accumulate density by applying sigmoid function
        // x axis
        

            for(int by = bin_y_min; by <= bin_y_max; ++by)
            {
                double bin_center_y = chip_bottom_ + (by + 0.5) * bin_height_;
                double dy = bin_center_y - mod_center_y;
                double sy = bell_shaped_potential(dy, mod_h, bin_height_);
                for(int bx = bin_x_min; bx <= bin_x_max; ++bx)
                {
                    double bin_center_x = chip_left_ + (bx + 0.5) * bin_width_;
                    double dx = bin_center_x - mod_center_x;
                    double sx = bell_shaped_potential(dx, mod_w, bin_width_);

                    double density_influence = sx * sy;
                    bin_density_[by][bx] += density_influence;
                }
                
            }
        
    }
    applyGaussianSmoothing(bin_density_, 5, 2);


    // calculate the value for the term density
    // density term sum(D_b(x, y) - M_b)^2
    // D_b, density for bin b
    // M_b, max density for bin b

    
    

    // Step 1: Normalize bin_density_ to [0, 1] range
    double p_min = numeric_limits<double>::max();
    double p_max = numeric_limits<double>::lowest();

    for (int by = 0; by < bin_rows_; ++by) {
        for (int bx = 0; bx < bin_cols_; ++bx) {
            double val = bin_density_[by][bx];
            p_min = min(p_min, val);
            p_max = max(p_max, val);
        }
    }

    // Avoid divide-by-zero
    double epsilon = 1e-8;
    double range = max(p_max - p_min, epsilon);
    norm_density.resize(bin_rows_, std::vector<double>(bin_cols_, 0.0));
    p_prime_prime.resize(bin_rows_, std::vector<double>(bin_cols_, 0.0));
    overflow_array.resize(bin_rows_, std::vector<double>(bin_cols_, 0.0));
    density_term_array.resize(bin_rows_, std::vector<double>(bin_cols_, 0.0));

    // Normalize to [0, 1]
    // norm_density(bin_rows_, vector<double>(bin_cols_));
    for (int by = 0; by < bin_rows_; ++by)
        for (int bx = 0; bx < bin_cols_; ++bx)
            norm_density[by][bx] = (bin_density_[by][bx] - p_min) / range;

    // Step 2: Compute normalized average
    double p_avg = 0.0;
    for (int by = 0; by < bin_rows_; ++by)
        for (int bx = 0; bx < bin_cols_; ++bx)
            p_avg += norm_density[by][bx];
    p_avg /= (bin_rows_ * bin_cols_);

    // Step 3: Apply level smoothing to normalized density
    delta_for_smoothing_ = 5.0;
    value_ = 0.0;

    for (int by = 0; by < bin_rows_; ++by) {
        for (int bx = 0; bx < bin_cols_; ++bx) {
            double p = norm_density[by][bx];
            double p_prime = 0.0;

            if (p >= p_avg)
                p_prime = p_avg + pow(p - p_avg, delta_for_smoothing_);
            else
                p_prime = p_avg - pow(p_avg - p, delta_for_smoothing_);

            p_prime_prime[by][bx] = p_prime;
            double overflow = (p_prime - target_density_);  // compare to target
            overflow_array[by][bx] = overflow;
            value_ += overflow * overflow;
            density_term_array[by][bx] = overflow * overflow;
            // cout << "norm_density[" << by << "][" << bx << "] = " << norm_density[by][bx] << ", p_prime = " << p_prime << ", overflow = " << overflow << endl;
        }
    }

    return value_ ;
}



const std::vector<Point2<double>> &Density::Backward() {
    const size_t num_modules = placement_.numModules();

    // Reset gradients
    for (auto &g : grad_)
        g = Point2<double>(0.0, 0.0);


    vector<vector<double>> normalized(bin_rows_, vector<double>(bin_cols_));
    for (int y = 0; y < bin_rows_; ++y) {
        for (int x = 0; x < bin_cols_; ++x) {
            normalized[y][x] = norm_density[y][x] / target_density_;  // or target_density_
        }
    }

    std::vector<std::vector<std::pair<double, double>>> grad_map(bin_rows_, std::vector<std::pair<double, double>>(bin_cols_));

    for (int y = 1; y < bin_rows_ - 1; ++y) {
        for (int x = 1; x < bin_cols_ - 1; ++x) {
            double dx = (norm_density[y][x + 1] - norm_density[y][x - 1]) / (2.0 * bin_width_);
            double dy = (norm_density[y + 1][x] - norm_density[y - 1][x]) / (2.0 * bin_height_);
            grad_map[y][x] = {dx, dy};
        }
    }
    
    for (size_t i = 0; i < num_modules; ++i) {
        Module &mod = placement_.module(i);
        if (mod.isFixed()) continue;

        double cx = input_[i].x;
        double cy = input_[i].y;
        double w = mod.width();
        double h = mod.height();
        double area = mod.area();

        double influence_range_x = w * 4.0;
        double influence_range_y = h * 4.0;

        double x_min = cx - influence_range_x / 2.0;
        double x_max = cx + influence_range_x / 2.0;
        double y_min = cy - influence_range_y / 2.0;
        double y_max = cy + influence_range_y / 2.0;

        int bin_x_min = std::max(1, (int)((x_min - chip_left_) / bin_width_));
        int bin_x_max = std::min(bin_cols_ - 2, (int)((x_max - chip_left_) / bin_width_));
        int bin_y_min = std::max(1, (int)((y_min - chip_bottom_) / bin_height_));
        int bin_y_max = std::min(bin_rows_ - 2, (int)((y_max - chip_bottom_) / bin_height_));

        for (int by = bin_y_min; by <= bin_y_max; ++by) {
            double bin_center_y = chip_bottom_ + (by + 0.5) * bin_height_;
            double dy = bin_center_y - cy;
            double wy = bell_shaped_potential(dy, h, bin_height_);

            for (int bx = bin_x_min; bx <= bin_x_max; ++bx) {
                double bin_center_x = chip_left_ + (bx + 0.5) * bin_width_;
                double dx = bin_center_x - cx;
                double wx = bell_shaped_potential(dx, w, bin_width_);

                double weight = wx * wy;

                grad_[i].x += weight * grad_map[by][bx].first * area;
                grad_[i].y += weight * grad_map[by][bx].second * area;
            }
        }
    }

    return grad_;
}




// const std::vector<Point2<double>> &Density::Backward() {
//     const size_t num_modules = placement_.numModules();

//     // Reset gradients
//     for (auto &g : grad_)
//         g = Point2<double>(0.0, 0.0);


//     vector<vector<double>> normalized(bin_rows_, vector<double>(bin_cols_));
//     for (int y = 0; y < bin_rows_; ++y) {
//         for (int x = 0; x < bin_cols_; ++x) {
//             normalized[y][x] = norm_density[y][x] / target_density_;  // or target_density_
//         }
//     }

//     std::vector<std::vector<std::pair<double, double>>> grad_map(bin_rows_, std::vector<std::pair<double, double>>(bin_cols_));


//     for (int by = 0; by < bin_rows_; ++by) {
//         double bin_center_y = chip_bottom_ + (by + 0.5) * bin_height_;
//         for (int bx = 0; bx < bin_cols_; ++bx) {
//             double bin_center_x = chip_left_ + (bx + 0.5) * bin_width_;

//             double dL_dD = 2.0 * (norm_density[by][bx] - target_density_);  // assuming L2 loss

//             for (size_t v = 0; v < placement_.numModules(); ++v) {
//                 Module &mod = placement_.module(v);
//                 if (mod.isFixed()) continue;

//                 double cx = input_[v].x;
//                 double cy = input_[v].y;
//                 double w = mod.width();
//                 double h = mod.height();

//                 double dx = bin_center_x - cx;
//                 double dy = bin_center_y - cy;

//                 double Px = bell_shaped_potential(dx, w, bin_width_);
//                 double dPx = bell_shaped_potential_derivative(dx, w, bin_width_);

//                 double Py = bell_shaped_potential(dy, h, bin_height_);
//                 double dPy = bell_shaped_potential_derivative(dy, h, bin_height_);

//                 // Skip if the influence is 0
//                 if (Px == 0 && Py == 0) continue;

//                 // Chain rule: ∂D/∂x = dPx * Py, ∂D/∂y = Px * dPy
//                 grad_map[by][bx].first += dL_dD * dPx * Py;
//                 grad_map[by][bx].second += dL_dD * Px * dPy;
//             }
//         }
//     }


//     vector<vector<pair<double, double>>> grad_map_norm(bin_rows_, vector<pair<double, double>>(bin_cols_));

//     double px_min = numeric_limits<double>::max();
//     double px_max = numeric_limits<double>::lowest();
//     double py_min = numeric_limits<double>::max();
//     double py_max = numeric_limits<double>::lowest();
//     double epsilon = 1e-8;

//     for (int by = 0; by < bin_rows_; ++by) {
//         for (int bx = 0; bx < bin_cols_; ++bx) {
//             double val1 = grad_map[by][bx].first;
//             double val2 = grad_map[by][bx].second;
//             px_min = min(px_min, val1);
//             px_max = max(px_max, val1);
//             py_min = min(py_min, val1);
//             py_max = max(py_max, val1);
//         }
//     }
//     double rangex = max(px_max - px_min, epsilon);
//     double rangey = max(py_max - py_min, epsilon);
    

//     for (int by = 0; by < bin_rows_; ++by)
//         for (int bx = 0; bx < bin_cols_; ++bx)
//         {
//             grad_map_norm[by][bx].first = (grad_map[by][bx].first - px_min) / rangex;
//             grad_map_norm[by][bx].second = (grad_map[by][bx].second - py_min) / rangey;

//         }
    
//     for (size_t i = 0; i < num_modules; ++i) {
//         Module &mod = placement_.module(i);
//         if (mod.isFixed()) continue;

//         double cx = input_[i].x;
//         double cy = input_[i].y;
//         double w = mod.width();
//         double h = mod.height();
//         double area = mod.area();

//         double influence_range_x = w * 4.0;
//         double influence_range_y = h * 4.0;

//         double x_min = cx - influence_range_x / 2.0;
//         double x_max = cx + influence_range_x / 2.0;
//         double y_min = cy - influence_range_y / 2.0;
//         double y_max = cy + influence_range_y / 2.0;

//         int bin_x_min = std::max(1, (int)((x_min - chip_left_) / bin_width_));
//         int bin_x_max = std::min(bin_cols_ - 2, (int)((x_max - chip_left_) / bin_width_));
//         int bin_y_min = std::max(1, (int)((y_min - chip_bottom_) / bin_height_));
//         int bin_y_max = std::min(bin_rows_ - 2, (int)((y_max - chip_bottom_) / bin_height_));

//         for (int by = bin_y_min; by <= bin_y_max; ++by) {
//             double bin_center_y = chip_bottom_ + (by + 0.5) * bin_height_;
//             double dy = bin_center_y - cy;
//             double wy = bell_shaped_potential(dy, h, bin_height_);

//             for (int bx = bin_x_min; bx <= bin_x_max; ++bx) {
//                 double bin_center_x = chip_left_ + (bx + 0.5) * bin_width_;
//                 double dx = bin_center_x - cx;
//                 double wx = bell_shaped_potential(dx, w, bin_width_);

//                 double weight = wx * wy;

//                 grad_[i].x += weight * grad_map_norm[by][bx].first * area;
//                 grad_[i].y += weight * grad_map_norm[by][bx].second * area;
//             }
//         }
//     }

//     return grad_;
// }



ObjectiveFunction::ObjectiveFunction(Placement &placement, double lambda, Wirelength wirelength, Density density)
    : BaseFunction(placement.numModules()),
        wirelength_(wirelength),  // set γ as needed
        density_(density),                       // default: 50×50 grid
        lambda_(lambda)/*,
        grad_(placement.numModules(), Point2<double>(0.0, 0.0))*/ {}

const double &ObjectiveFunction::operator()(const std::vector<Point2<double>> &input) {
        input_ = input;                     // cache for Backward()
        wirelength_(input);                // ensure internal input_ is set
        density_(input);                   // ensure internal input_ is set
        const double wl = wirelength_.value();
        const double dp = density_.value();
        // cout << "wl: " << wl << ", dp: " << dp << std::endl;
        value_ = wl + lambda_ * dp;
        return value_;
}



const std::vector<Point2<double>> &ObjectiveFunction::Backward() {
    // Store gradients in temporary vectors to prevent overwrites
    auto grad_wl = wirelength_.Backward();  // Make a copy
    auto grad_dp = density_.Backward();     // Make a copy

    for (size_t i = 0; i < grad_.size(); ++i) {
        grad_[i].x = grad_wl[i].x + lambda_ * grad_dp[i].x;
        grad_[i].y = grad_wl[i].y + lambda_ * grad_dp[i].y;
        // cout << "grad_wl[" << i << "]: " << grad_wl[i].x << ", " << grad_wl[i].y << std::endl;
        // cout << "lambda * grad_dp[" << i << "]: " << lambda_ * grad_dp[i].x << ", " << lambda_ * grad_dp[i].y << std::endl;
    }

    return grad_;
}

void ObjectiveFunction::setLambda(double lambda) {
    lambda_ = lambda;
}

double ObjectiveFunction::getLambda() const {
    return lambda_;
}
