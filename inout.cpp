#include <fstream>
using namespace std;

void readdata (char *f, int *x, int n) {
   ifstream in (f, ios::in|ios::binary);
   in.read ((char *)x, n * sizeof(int));
   in.close();
}

void readdata (char *f, float *x, int n) {
   ifstream in (f, ios::in|ios::binary);
   in.read ((char *)x, n * sizeof(float));
   in.close();
}

void readdata (char *f, double *x, int n) {
   ifstream in (f, ios::in|ios::binary);
   in.read ((char *)x, n * sizeof(double));
   in.close();
}

void writedata (char *f, int *x, int n) {
   ofstream out (f, ios::out|ios::binary);
   out.write ((char *)x, n * sizeof(int));
   out.close();
}

void writedata (char *f, float *x, int n) {
   ofstream out (f, ios::out|ios::binary);
   out.write ((char *)x, n * sizeof(float));
   out.close();
}

void writedata (char *f, double *x, int n) {
   ofstream out (f, ios::out|ios::binary);
   out.write ((char *)x, n * sizeof(double));
   out.close();
}
