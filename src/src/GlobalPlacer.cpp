#include "GlobalPlacer.h"

#include <cstdio>
#include <vector>
#include <set>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"

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
    const size_t num_modules = _placement.numModules();
    std::vector<Point2<double>> t(num_modules);
    for (size_t i = 0; i < num_modules; ++i) {
        if (_placement.module(i).isFixed()) continue;
        t[i] = Point2<double>(_placement.module(i).centerX(), _placement.module(i).centerY());
        // cout << _placement.module(i).name() << " (" << t[i].x << ", " << t[i].y << ")" << endl;
    }

   
    // ExampleFunction foo(_placement);                    // Objective function
    ObjectiveFunction obj(_placement, /*lambda=*/0.0001);

    const double kAlpha = 10;                         // Constant step size
    SimpleConjugateGradient optimizer(obj, t, kAlpha);  // Optimizer

    // Set initial point
    // t[0] = 4.;  // This set both t[0].x and t[0].y to 4.

    // Initialize the optimizer
    optimizer.Initialize();

    Wirelength wirelength_(_placement, /*gamma=*/5000.0);  // Wirelength function
    Density density_(_placement, /*bin_rows=*/50, /*bin_cols=*/50, /*sigma_factor=*/1.5, /*target_density=*/0.9);  // Density function
    // Perform optimization, the termination condition is that the number of iterations reaches 100
    // TODO: You may need to change the termination condition, which is determined by the overflow ratio.
    for (size_t i = 0; i < 500; ++i) {
        optimizer.Step();
        // printf("iter = %3lu, f = %9.4f, x = %9.4f, y = %9.4f\n", i, obj(t), t[0].x, t[0].y);
        printf("iter = %3lu, objective function value = %9.4f\n", i, obj(t));
        for (size_t j = 0; j < t.size(); ++j) {
            if (_placement.module(j).isFixed()) continue;
            // printf("iter %2lu, module %2zu: pos = (%.2f, %.2f)\n", i, j, t[i].x, t[i].y);
        }

        // printf("iter %3lu: value = %.2f, WL = %.2f, DP = %.2f\n", i, obj(t), wirelength_(t), density_(t));

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
        // _placement.module(i).setPosition(10, 10);

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
