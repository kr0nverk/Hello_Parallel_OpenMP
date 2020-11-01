#include <iostream>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <iomanip>

using namespace std;


int main() {
/*
    ofstream fout("C:\\Users\\vovan\\source\\repos\\Hello_Parallel\\RunForestRun\\a.txt");
    if (fout.is_open()) {
        double temp = 123456789.;
        fout << fixed << setprecision(0) << temp;
    }
    fout.close();
    */

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            char temp;
            int count = 0;
            ifstream inf("C:\\Users\\vovan\\source\\repos\\Hello_Parallel\\RunForestRun\\a.txt");
            if (inf.is_open()) {
                while (inf.get(temp) && count < 2) {
                    cout << temp;
                    count++;
                }
                inf.close();
            }
        }
    }
	return 0;
}