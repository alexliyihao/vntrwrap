/*
compile:
g++ \
-O3 \
-Wall \
-std=c++03 \
find_neighbors.cpp \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_software/Sources/Eagle_v2.4.1/src \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/include \
-Wl,-rpath,/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/lib \
-o find_neighbors \
-L /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/lib \
-l boost_iostreams \
-lz

running:
./find_neighbors \
int batch number\
int total batches \
float Max z range \
path txt.gz output from normalize_mosdepth \
output path, which is also a .txt.gz file > keep stdout if you need

This code is originally designed for very large batches so pleaes stick with 0 and 1, use
e.g.
./find_neighbors 0 1 2.0 \
/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/mosdepth/GRCh37_ID_scale_zdepths.txt.gz \
/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/find_neighbor/GRCh37_find_neighbor > running_log_find_neighbor_GRCh37_zmax=2.0.txt
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>

#include "FileUtils.cpp"
#include "StringUtils.cpp"
#include "Timer.cpp"

using namespace std;

inline float sq(float x) { return x*x; }

float crop(float z, float zMax) {
  return min(zMax, max(-zMax, z));
}

int main(int argc, char *argv[]) {

  if (argc != 6) {
    cerr << "ERROR: 5 arguments required" << endl;
    cerr << "- arg1: batch number" << endl;
    cerr << "- arg2: total batches" << endl;
    cerr << "- arg3: Max z range" << endl;
    cerr << "- arg4: ID_scale_zdepths file" << endl;
    cerr << "- arg5: output prefix" << endl;
    /*
    cerr << "- arg3: fraction of regions to use" << endl;
    cerr << "- arg4: zMax: crop value for z-scores" << endl;
    */
    return 1;
  }
  int b; sscanf(argv[1], "%d", &b);
  int B; sscanf(argv[2], "%d", &B);
  const float fracR = 1; //sscanf(argv[3], "%f", &fracR);
  float zMax; sscanf(argv[3], "%f", &zMax);
  const char *dataFile = argv[4]; //"ID_scale_zdepths.txt.gz";
  const char *outPrefix = argv[5]; //"ID_scale_neighbors";

  cout << "Computing nearest neighbors for batch " << b << " mod " << B << endl;
  cout << "Cropping 'z-score' values to zMax = " << zMax << endl;

  Timer timer;

  // read numbers of individuals and regions
  FileUtils::AutoGzIfstream fin; fin.openOrExit(dataFile);
  int N, R; fin >> N >> R; // N = 3915, R = 939
  int Nbatch = 0;
  // divide the total N by total number of batches, use the modulo as the seperator
  for (int n = 0; n < N; n++)
    if (n % B == b)
      Nbatch++;
  cout << "the actual N_batch is " << Nbatch << endl;

  string line; getline(fin, line); // throw away first header line: mu

  fin >> N >> R;
  vector <float> sigma2ratios(R);
  for (int r = 0; r < R; r++)
    fin >> sigma2ratios[r];

  vector <float> sigma2sort = sigma2ratios;
  sort(sigma2sort.begin(), sigma2sort.end());
  float sigma2min = sigma2sort[(int) (R*(1-fracR))];
  float sigma2max = 1000; //sigma2sort[0] * 100;
  int Ruse = 0, Rextreme = 0;
  for (int r = 0; r < R; r++) {
    if (sigma2ratios[r] >= sigma2min && sigma2ratios[r] <= sigma2max)
      Ruse++;
    if (sigma2ratios[r] > sigma2max)
      Rextreme++;
  }
  cout << "Removed " << Rextreme << " of " << R << " regions with sigma2ratio > " << sigma2max
       << endl;
  cout << "Keeping " << Ruse << " of " << R-Rextreme << " remaining regions with sigma2ratio >= "
       << sigma2min << endl;

  cout << "Reading data for " << Nbatch << " / " << N << " indivs in batch at " << R << " regions"
       << endl;

  // allocate memory for individuals in batch
  vector <string> IDs(N);
  vector <float> scales(N);
  float *zs = new float[R*(long long) Nbatch];

  // store z-scores for individuals in batch
  string ID; float scale, z;
  for (int n = 0; n < N; n++) {
    fin >> ID >> scale;
    IDs[n] = ID;
    scales[n] = scale;
    if (n % B == b) {
      int i = n/B;
      for (int r = 0; r < R; r++) {
      	fin >> z;
      	zs[r*Nbatch+i] = crop(z, zMax);
        // cout << "cropped zs: "<<  zs[r*Nbatch+i] << endl;
      }
    }
    else
      getline(fin, line);
    if (n%100==0)
      cout << "." << flush;
  }
  fin.close();
  cout << endl << "Read data for " << Nbatch << " / " << N << " indivs in batch ("
       << timer.update_time() << " sec)" << endl;

  // stream z-scores for all individuals
  float *dists = new float[N*(long long) Nbatch];
  memset(dists, 0, N*(long long) Nbatch*sizeof(dists[0]));
  fin.openOrExit(dataFile);
  getline(fin, line); getline(fin, line); // throw away header lines
  for (int n = 0; n < N; n++) {
    fin >> ID >> scale; // already stored
    for (int r = 0; r < R; r++) {
      fin >> z;
      if (sigma2ratios[r] < sigma2min || sigma2ratios[r] > sigma2max) continue;
      z = crop(z, zMax);
      for (int i = 0; i < Nbatch; i++){
        dists[n*Nbatch+i] += sq(z - zs[r*Nbatch+i]);
        // cout << "computed a distance" <<  dists[n*Nbatch+i] <<endl;
      }
    }
    if (n%100==0)
      cout << "." << flush;
  }
  fin.close();
  cout << endl << "Computed distances for " << Nbatch << " / " << N << " indivs in batch ("
       << timer.update_time() << " sec)" << endl;

  char buf[200];
  //sprintf(buf, "ID_scale_neighbors.top%g.zMax%g.batch%d.txt.gz", fracR, zMax, b);
  sprintf(buf, "%s.zMax%g.txt.gz", outPrefix, zMax);
  FileUtils::AutoGzOfstream fout; fout.openOrExit(buf);
  fout << std::setprecision(2) << std::fixed;
  vector < pair <float, int> > distIDs(N); // this ID is INDEX of IDS vector
  for (int i = 0; i < Nbatch; i++) {
    int n_i = i*B+b;
    fout << IDs[n_i] << "\t" << scales[n_i];
    // sort and output best matches
    for (int n = 0; n < N; n++)
      distIDs[n] = make_pair(dists[n*Nbatch+i], n);
    distIDs[n_i].first = 1e9;
    sort(distIDs.begin(), distIDs.end());

    const int N_OUTPUT = 500;
    for (int j = 0; j < N_OUTPUT; j++) {
      int n = distIDs[j].second;
      fout << "\t" << IDs[n] << "\t" << scales[n] << "\t" << dists[n*Nbatch+i]/(2*Ruse);
    }
    fout << endl;
  }
  fout.close();
  cout << "Found neighbors and wrote output (" << timer.update_time() << " sec)" << endl;

  delete[] dists;
  delete[] zs;
}
