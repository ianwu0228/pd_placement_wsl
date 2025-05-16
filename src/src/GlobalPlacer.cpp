#include "GlobalPlacer.h"

#include <cstdio>
#include <vector>
#include <set>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"
#include <random>

GlobalPlacer::GlobalPlacer(Placement &placement)
    : _placement(placement) {
}

void GlobalPlacer::place() {
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
    double offset_x = _placement.boundryRight() * 0.05;  // 5% of chip width
    double offset_y = _placement.boundryTop() * 0.05;    // 5% of chip height
    std::uniform_real_distribution<> dis_x(-offset_x, offset_x);
    std::uniform_real_distribution<> dis_y(-offset_y, offset_y);

    // Center coordinates
    double center_x = _placement.boundryRight() / 2;
    double center_y = _placement.boundryTop() / 2;

    // std::vector<Point2<double>> t(num_modules);
    for (size_t i = 0; i < num_modules; ++i) {
        if (_placement.module(i).isFixed()) continue;
        
        // Add small random offset to center position
        t[i] = Point2<double>(
            center_x + dis_x(gen),
            center_y + dis_y(gen)
        );
        
        cout << _placement.module(i).name() << " (" << t[i].x << ", " << t[i].y << ")" << endl;
    }

   
    // ExampleFunction foo(_placement);                    // Objective function
    ObjectiveFunction obj(_placement, /*lambda=*/100000);

    const double kAlpha = 10;                         // Constant step size
    SimpleConjugateGradient optimizer(obj, t, kAlpha);  // Optimizer

    // Set initial point
    // t[0] = 4.;  // This set both t[0].x and t[0].y to 4.

    // Initialize the optimizer
    optimizer.Initialize();

    Wirelength wirelength_(_placement, /*gamma=*/5000.0);  // Wirelength function
    Density density_(_placement, /*bin_rows=*/100, /*bin_cols=*/100, /*sigma_factor=*/1.5, /*target_density=*/0.9);  // Density function
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
    for (size_t i = 0; i < 100; ++i) {
        optimizer.Step();
        
        // Get current density map
        vector<vector<double>> bin_density = density_.getBinDensity();
        // Plot every 10 iterations
        if (i % 10 == 0) {
            // Create output directory
            system("mkdir -p plot_output");
            
            ////////////////////////////////// Density Map Plot //////////////////////////////////
            string densityname = "plot_output/density_" + std::to_string(i) + ".plt";
            string densitypng = "plot_output/density_" + std::to_string(i) + ".png";
            
            ofstream densityfile(densityname.c_str(), ios::out);
            densityfile << "set terminal png size 800,800 enhanced font 'Arial,12'" << endl;
            densityfile << "set output '" << densitypng << "'" << endl;
            densityfile << "set title \"Density Map - Iteration " << i << "\"" << endl;
            densityfile << "set view map" << endl;
            densityfile << "set size ratio 1" << endl;
            densityfile << "unset key" << endl;
            densityfile << "set palette defined (0 'white', 0.5 'yellow', 1 'red', 2 'dark-red')" << endl;
            densityfile << "set cbrange [0:2]" << endl;
            densityfile << "set cblabel 'Density'" << endl;
            densityfile << "set xrange [0:" << bin_density[0].size()-1 << "]" << endl;
            densityfile << "set yrange [0:" << bin_density.size()-1 << "]" << endl;
            
            // Write density data
            densityfile << "$data << EOD" << endl;
            for (size_t y = 0; y < bin_density.size(); ++y) {
                for (size_t x = 0; x < bin_density[y].size(); ++x) {
                    densityfile << x << " " << y << " " << bin_density[y][x] << endl;
                }
                densityfile << endl;
            }
            densityfile << "EOD" << endl;
            
            densityfile << "plot '$data' using 1:2:3 with image" << endl;
            densityfile.close();
            
            ////////////////////////////////// Cell Distribution Plot //////////////////////////////////
            string cellname = "plot_output/cells_" + std::to_string(i) + ".plt";
            string cellpng = "plot_output/cells_" + std::to_string(i) + ".png";
            
            ofstream cellfile(cellname.c_str(), ios::out);
            cellfile << "set terminal png size 800,800 enhanced font 'Arial,12'" << endl;
            cellfile << "set output '" << cellpng << "'" << endl;
            cellfile << "set title \"Cell Distribution - Iteration " << i 
                    << "\\nWL = " << wirelength_(t) << ", DP = " << density_(t) << "\"" << endl;
            cellfile << "set size ratio 1" << endl;
            cellfile << "set xrange [" << _placement.boundryLeft() << ":" 
                    << _placement.boundryRight() << "]" << endl;
            cellfile << "set yrange [" << _placement.boundryBottom() << ":" 
                    << _placement.boundryTop() << "]" << endl;
            
            // Set point styles
            cellfile << "set style line 1 lc rgb 'red' pt 7 ps 0.3" << endl;
            cellfile << "set style line 2 lc rgb 'blue' pt 7 ps 0.3" << endl;
            cellfile << "set style line 3 lc rgb 'black' lt 1 lw 2" << endl;
            
            // Plot command
            cellfile << "plot '-' w p ls 1 title 'Fixed', "
                    << "     '-' w p ls 2 title 'Movable', "
                    << "     '-' w l ls 3 title 'Boundary'" << endl;
            
            // Plot fixed modules
            for (size_t j = 0; j < t.size(); ++j) {
                if (_placement.module(j).isFixed()) {
                    cellfile << t[j].x << " " << t[j].y << endl;
                }
            }
            cellfile << "e" << endl;
            
            // Plot movable modules
            for (size_t j = 0; j < t.size(); ++j) {
                if (!_placement.module(j).isFixed()) {
                    cellfile << t[j].x << " " << t[j].y << endl;
                }
            }
            cellfile << "e" << endl;
            
            // Plot boundary
            plotBoxPLT(cellfile, _placement.boundryLeft(), _placement.boundryBottom(), 
                    _placement.boundryRight(), _placement.boundryTop());
            cellfile << "e" << endl;
            cellfile.close();

            // Execute gnuplot
            char cmd[256];
            sprintf(cmd, "gnuplot %s %s", densityname.c_str(), cellname.c_str());
            system(cmd);

            // Combine the two plots side by side
            sprintf(cmd, "convert +append %s %s plot_output/combined_%zu.png", 
                    densitypng.c_str(), cellpng.c_str(), i);
            system(cmd);

            printf("Generated plots for iteration %zu\n", i);
        }
    }

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
    for (size_t i = 0; i < num_modules; i++) {
        // If the module is fixed, its position should not be changed.
        // In this programing assignment, a fixed module may be a terminal or a pre-placed module.
        if (_placement.module(i).isFixed()) {
            fixed_cnt++;
            continue;
        }
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

        cout << _placement.module(i).name() << " (" << _placement.module(i).centerX() << ", " << _placement.module(i).centerY() << ")" << endl;
    }
    printf("INFO: %lu / %lu modules are fixed.\n", fixed_cnt, num_modules);
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


