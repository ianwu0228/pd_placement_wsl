#include "GlobalPlacer.h"

#include <cstdio>
#include <vector>
#include <set>
#include <algorithm>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"
#include <random>

GlobalPlacer::GlobalPlacer(Placement &placement)
    : _placement(placement) {
}

void GlobalPlacer::place() {

    // pre-process "shift the bottom-left cell to (0,0)"


    ////////////////////////////////////////////////////////////////////
    // This section is an example for analytical methods.
    // The objective is to minimize the following function:
    //      f(x,y) = 3*x^2 + 2*x*y + 2*y^2 + 7
    //
    // If you use other methods, you can skip and delete it directly.
    ////////////////////////////////////////////////////////////////////
    // std::vector<Point2<double>> t(1);                   // Optimization variables (in this example, there is only one t)
    // const size_t num_modules = _placement.numModules();

    const size_t num_modules = _placement.numModules();
    std::vector<Point2<double>> t(num_modules);

    // Initialize random number generator once outside the loop
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create small random offsets (±5% of chip size around center)
    double offset_x = (_placement.boundryRight() - _placement.boundryLeft()) * 0.05;  // 5% of chip width
    double offset_y = (_placement.boundryTop() - _placement.boundryBottom())* 0.05;    // 5% of chip height
    std::uniform_real_distribution<> dis_x(-offset_x, offset_x);
    std::uniform_real_distribution<> dis_y(-offset_y, offset_y);

    // Center coordinates
    double center_x = (_placement.boundryRight() + _placement.boundryLeft()) / 2;
    double center_y = (_placement.boundryTop() + _placement.boundryBottom()) / 2;
    cout << "center_x: " << center_x << ", center_y: " << center_y << endl;
    // std::vector<Point2<double>> t(num_modules);
    for (size_t i = 0; i < num_modules; ++i) {
        if (_placement.module(i).isFixed()) continue;
        
        // Add small random offset to center position
        t[i] = Point2<double>(
            center_x + dis_x(gen),
            center_y + dis_y(gen)
        );
        
        // print out the original position of all cells
        cout << _placement.module(i).name() << " (" << t[i].x << ", " << t[i].y << ")" << endl;
    }
   

    // int bin_rows, bin_cols = (int)((_placement.boundryRight() - _placement.boundryLeft())/3); // 800 for ibm05
    // int bin_rows = (int)((_placement.boundryRight() - _placement.boundryLeft())/3);
    // int bin_cols = (int)((_placement.boundryRight() - _placement.boundryLeft())/3); 
    int bin_rows = 200;
    int bin_cols = 200;

    Wirelength wirelength_(_placement, /*gamma=*/700.0);  // Wirelength function
    Density density_(_placement, /*bin_rows=*/bin_rows, /*bin_cols=*/bin_cols, /*sigma_factor=*/1.5, /*target_density=*/0.9);  // Density function
    ObjectiveFunction obj(_placement, /*lambda=*/0.0000000001, wirelength_, density_);

    const double kAlpha = 5;                         // Constant step size
    SimpleConjugateGradient optimizer(obj, t, kAlpha, _placement.boundryLeft(), _placement.boundryRight(), _placement.boundryTop(), _placement.boundryBottom());  // Optimizer

    // Initialize the optimizer
    optimizer.Initialize();


    // Perform optimization, the termination condition is that the number of iterations reaches 100
    // TODO: You may need to change the termination condition, which is determined by the overflow ratio.
    // for (size_t i = 0; i < 100; ++i) {
    //     optimizer.Step();
    //     // printf("iter = %3lu, f = %9.4f, x = %9.4f, y = %9.4f\n", i, obj(t), t[0].x, t[0].y);
    //     // printf("iter = %3lu, objective function value = %9.4f\n", i, obj(t));
    //     // for (size_t j = 0; j < t.size(); ++j) {
    //     //     if (_placement.module(j).isFixed()) continue;
    //     //     // printf("iter %2lu, module %2zu: pos = (%.2f, %.2f)\n", i, j, t[i].x, t[i].y);
    //     // }

    //     // printf("iter %3lu: value = %.2f, WL = %.2f, DP = %.2f\n", i, obj(t), wirelength_(t), density_(t));
    //     //plot the density during the optimization

    // }
    double init_lambda = 0.0000000001;
    double lambda = init_lambda;
    int i = -1;
    do{
        i++;
        lambda *= 2;
        obj.setLambda(max(lambda, init_lambda * 4000));
        double min_density = std::numeric_limits<double>::max();
        double max_density = std::numeric_limits<double>::lowest();
        vector<vector<double>> bin_density = density_.getBinDensity();
        if (i % 1 == 0) {
            // Create output directory

            int bin_rows = bin_density.size();
            int bin_cols = bin_density[0].size();

            // Histogram bins
            int count_0_05 = 0, count_05_1 = 0, count_1_15 = 0, count_15_2 = 0;
            int count_2_10 = 0, count_10_20 = 0, count_20_50 = 0, count_50_100 = 0, count_over_100 = 0;

            for (int y = 0; y < bin_rows; ++y) {
                for (int x = 0; x < bin_cols; ++x) {
                    double d = bin_density[y][x];

                    // Update min/max
                    min_density = std::min(min_density, d);
                    max_density = std::max(max_density, d);

                    // Histogram classification
                    if (d < 0.5) count_0_05++;
                    else if (d < 1.0) count_05_1++;
                    else if (d < 1.5) count_1_15++;
                    else if (d < 2.0) count_15_2++;
                    else if (d < 10.0) count_2_10++;
                    else if (d < 20.0) count_10_20++;
                    else if (d < 50.0) count_20_50++;
                    else if (d < 100.0) count_50_100++;
                    else count_over_100++;
                }
            }
            cout << "Density Histogram:" << endl;
            // cout << "  [0.0 ~ 0.5)    : " << count_0_05 << " bins" << endl;
            // cout << "  [0.5 ~ 1.0)    : " << count_05_1 << " bins" << endl;
            // cout << "  [1.0 ~ 1.5)    : " << count_1_15 << " bins" << endl;
            // cout << "  [1.5 ~ 2.0)    : " << count_15_2 << " bins" << endl;
            // cout << "  [2.0 ~ 10.0)   : " << count_2_10 << " bins" << endl;
            // cout << "  [10.0 ~ 20.0)  : " << count_10_20 << " bins" << endl;
            // cout << "  [20.0 ~ 50.0)  : " << count_20_50 << " bins" << endl;
            // cout << "  [50.0 ~ 100.0) : " << count_50_100 << " bins" << endl;
            // cout << "  [>100.0]       : " << count_over_100 << " bins" << endl;
            // cout << "Min density      : " << min_density << endl;
            cout << "Max density      : " << max_density << endl;
            cout << "wirelength = " << wirelength_(t) << endl;



            ////////////////////////////////////////////////////////////////////////////////////////
            // Check if we have fixed and movable modules
            
            // system("mkdir -p plot_output");
            
            // ////////////////////////////////// Density Map Plot //////////////////////////////////
            // string densityname = "plot_output/density_" + std::to_string(i) + ".plt";
            // string densitypng = "plot_output/density_" + std::to_string(i) + ".png";
            
            // ofstream densityfile(densityname.c_str(), ios::out);
            // densityfile << "set terminal png size 800,800 enhanced font 'Arial,12'" << endl;
            // densityfile << "set output '" << densitypng << "'" << endl;
            // densityfile << "set title \"Density Map - Iteration " << i << "\"" << endl;
            // densityfile << "set view map" << endl;
            // densityfile << "set size ratio 1" << endl;
            // densityfile << "unset key" << endl;
            // densityfile << "set palette defined (0 'white', 0.5 'yellow', 1 'red', 2 'dark-red')" << endl;
            // densityfile << "set cbrange [0:2]" << endl;
            // densityfile << "set cblabel 'Density'" << endl;
            // densityfile << "set xrange [0:" << bin_density[0].size()-1 << "]" << endl;
            // densityfile << "set yrange [0:" << bin_density.size()-1 << "]" << endl;

            
        
            // densityfile << "set pm3d map" << endl;
            // densityfile << "splot '-' using 1:2:3 notitle" << endl;

            // for (size_t y = 0; y < bin_density.size(); ++y) {
            //     for (size_t x = 0; x < bin_density[y].size(); ++x) {
            //         densityfile << x << " " << y << " " << bin_density[y][x] << endl;
            //     }
            //     densityfile << endl;
            // }
            // densityfile << "e" << endl;

            // densityfile.close();
            
            // bool hasFixed = false;
            // bool hasMovable = false;
            // for (size_t j = 0; j < t.size(); ++j) {
            //     if (_placement.module(j).isFixed()) {
            //         hasFixed = true;
            //     } else {
            //         hasMovable = true;
            //     }
            //     if (hasFixed && hasMovable) break;
            // }

            // // Cell Distribution Plot
            // string cellname = "plot_output/cells_" + std::to_string(i) + ".plt";
            // string cellpng = "plot_output/cells_" + std::to_string(i) + ".png";

            // ofstream cellfile(cellname.c_str(), ios::out);
            // cellfile << "set terminal png size 800,800 enhanced font 'Arial,12'" << endl;
            // cellfile << "set output '" << cellpng << "'" << endl;
            // cellfile << "set title \"Cell Distribution - Iteration " << i 
            //         << "\\nWL = " << wirelength_(t) << ", DP = " << density_(t) << "\"" << endl;
            // cellfile << "set size ratio 1" << endl;
            // cellfile << "set xrange [" << _placement.boundryLeft() << ":" 
            //         << _placement.boundryRight() << "]" << endl;
            // cellfile << "set yrange [" << _placement.boundryBottom() << ":" 
            //         << _placement.boundryTop() << "]" << endl;

            // // Set point styles
            // cellfile << "set style line 1 lc rgb 'red' pt 7 ps 0.3" << endl;
            // cellfile << "set style line 2 lc rgb 'blue' pt 7 ps 0.3" << endl;
            // cellfile << "set style line 3 lc rgb 'black' lt 1 lw 2" << endl;

            // // Construct plot command based on what types of modules exist
            // string plotCmd = "plot ";
            // if (hasFixed) {
            //     plotCmd += "'-' w p ls 1 title 'Fixed'";
            //     if (hasMovable) plotCmd += ", ";
            // }
            // if (hasMovable) {
            //     plotCmd += "'-' w p ls 2 title 'Movable'";
            // }
            // plotCmd += ", '-' w l ls 3 title 'Boundary'";
            // cellfile << plotCmd << endl;

            // // Plot fixed modules if they exist
            // if (hasFixed) {
            //     for (size_t j = 0; j < t.size(); ++j) {
            //         if (_placement.module(j).isFixed()) {
            //             cellfile << t[j].x << " " << t[j].y << endl;
            //         }
            //     }
            //     cellfile << "e" << endl;
            // }

            // // Plot movable modules if they exist
            // if (hasMovable) {
            //     for (size_t j = 0; j < t.size(); ++j) {
            //         if (!_placement.module(j).isFixed()) {
            //             cellfile << t[j].x << " " << t[j].y << endl;
            //         }
            //     }
            //     cellfile << "e" << endl;
            // }

            // // Always plot boundary
            // plotBoxPLT(cellfile, _placement.boundryLeft(), _placement.boundryBottom(), 
            //         _placement.boundryRight(), _placement.boundryTop());
            // cellfile << "e" << endl;
            // cellfile.close();

            // // Execute gnuplot
            // char cmd[256];
            // sprintf(cmd, "gnuplot %s %s", densityname.c_str(), cellname.c_str());
            // system(cmd);

            // // Combine the two plots side by side
            // sprintf(cmd, "convert +append %s %s plot_output/combined_%zu.png", 
            //         densitypng.c_str(), cellpng.c_str(), i);
            // system(cmd);

            // printf("Generated plots for iteration %zu\n", i);
            /////////////////////////////////////////////////////////////////////////////////////////////
        }

        optimizer.Step();
        density_(t);
        density_.setSmoothingDelta(max(1.0,  density_.getSmoothingDelta() * 0.9));
        if(i != 0 && max_density < 12) break;
        
    }while(true);

    ////////////////////////////////////////////////////////////////////
    // Global placement algorithm
    ////////////////////////////////////////////////////////////////////

    // TODO: Implement your global placement algorithm here.
    // const size_t num_modules = _placement.numModules();  // You may modify this line.
    std::vector<Point2<double>> positions(num_modules);  // Optimization variables (positions of modules). You may modify this line.

    ////////////////////////////////////////////////////////////////////
    // Write the placement result into the database. (You may modify this part.)
    ////////////////////////////////////////////////////////////////////
    size_t fixed_cnt = 0;
    int num_in_chip_boundary = 0;
    for (size_t i = 0; i < num_modules; i++) {
        // If the module is fixed, its position should not be changed.
        // In this programing assignment, a fixed module may be a terminal or a pre-placed module.
        if (_placement.module(i).isFixed()) {
            fixed_cnt++;
            continue;
        }

        if(t[i].x < _placement.boundryLeft() || t[i].x > _placement.boundryRight() ||
           t[i].y < _placement.boundryBottom() || t[i].y > _placement.boundryTop()) {
            num_in_chip_boundary++;
        }


        // clamp the position of the module to the chip boundary
        t[i].x = max(_placement.boundryLeft() ,  min(t[i].x, _placement.boundryRight()));
        t[i].y = max(_placement.boundryBottom(), min(t[i].y, _placement.boundryTop()));


        _placement.module(i).setPosition(t[i].x, t[i].y);
        
        

        // /////////////////////////////////// random placement start ///////////////////////
        // std::random_device rd;
        // std::mt19937 gen(rd());
        // std::uniform_real_distribution<> dis_x(_placement.boundryLeft(), _placement.boundryRight());
        // std::uniform_real_distribution<> dis_y(_placement.boundryBottom(), _placement.boundryTop());
        // // Get module dimensions to ensure it stays within boundaries
        // double width = _placement.module(i).width();
        // double height = _placement.module(i).height();
        
        // // Generate random position while keeping the module within chip boundaries
        // double x = dis_x(gen);
        // double y = dis_y(gen);
        
        // // Adjust if the module would go outside boundaries
        // x = std::min(x, _placement.boundryRight() - width);
        // y = std::min(y, _placement.boundryTop() - height);
        // x = std::max(x, _placement.boundryLeft());
        // y = std::max(y, _placement.boundryBottom());
        
        // _placement.module(i).setPosition(x, y);
        // /////////////////////////////////// random placement end ///////////////////////


        // print out the final position of all cells
        // cout << _placement.module(i).name() << " (" << _placement.module(i).centerX() << ", " << _placement.module(i).centerY() << ")" << endl;
    }
    printf("INFO: %lu / %lu modules are fixed.\n", fixed_cnt, num_modules);
    cout << "INFO: " << num_in_chip_boundary << " modules are out of chip boundary." << endl;


}

void GlobalPlacer::plotPlacementResult(const string outfilename, bool isPrompt) {
    ofstream outfile(outfilename.c_str(), ios::out);
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl
            << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl
            << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT(outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop());
    outfile << "EOF" << endl;
    outfile << "# modules" << endl
            << "0.00, 0.00" << endl
            << endl;
    for (size_t i = 0; i < _placement.numModules(); ++i) {
        Module &module = _placement.module(i);
        plotBoxPLT(outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height());
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if (isPrompt) {
        char cmd[200];
        sprintf(cmd, "gnuplot %s", outfilename.c_str());
        if (!system(cmd)) {
            cout << "Fail to execute: \"" << cmd << "\"." << endl;
        }
    }
}

void GlobalPlacer::plotBoxPLT(ofstream &stream, double x1, double y1, double x2, double y2) {
    stream << x1 << ", " << y1 << endl
           << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl
           << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl
           << endl;
}


