#include <iostream>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <iomanip>

using namespace std;


int main() {

	ofstream fout;
    fout.open("C:\\Users\\vovan\\source\\repos\\Hello_Parallel\\RunForestRun\\a.txt");

    if (fout.is_open())
    {
        double temp = 123456789.;
        fout << fixed << setprecision(0) << temp;
    }

    #pragma omp parallel section
    {
        #pragma omp section
        {

        }
        #pragma omp section
        {

        }
    }

	return 0;
}