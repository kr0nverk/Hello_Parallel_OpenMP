#include <iostream>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <iomanip>

using namespace std;


int main() {
    /*
    ofstream fout("D:\\Documents\\Source\\repos\\Hello_Parallel_OpenMP\\RunForestRun\\a.txt");
    if (fout.is_open()) {
        double temp = 8883456789.;
        fout << fixed << setprecision(0) << temp;
    }
    fout.close();

    char temp;
            int count = 0;
            ifstream inf("D:\\Documents\\Source\\repos\\Hello_Parallel_OpenMP\\RunForestRun\\a.txt");
            if (inf.is_open()) {
                while (inf.get(temp) && count < 2) {
                    cout << temp;
                    count++;
                }
                inf.close();
            }
 */

    char temp2;
    string file_name = "D:\\Documents\\Source\\repos\\Hello_Parallel_OpenMP\\RunForestRun\\a.txt";
    ifstream in(file_name);

    if (in.is_open()) {
#pragma omp parallel sections
        {

#pragma omp section
            {
                while (in.get(temp2)) {
#pragma omp ordered
                }
            }

#pragma omp section
            {

            }
        }
        in.close();
    }
	return 0;
}