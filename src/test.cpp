#include<stdio.h>
using namespace std;
#include <iostream>
#include<stdlib.h>
#include <vector>

#include <iostream>


class MyFunc {
private:
    double value_;
public:
    const double& operator()(double x) {
        value_ = x * x;  // modifies internal state
        return value_;   // gives a read-only reference to caller
    }
};


int main()
{
    MyFunc f;




}