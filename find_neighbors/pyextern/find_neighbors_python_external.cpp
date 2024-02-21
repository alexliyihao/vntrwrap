/*
g++ \
-O3 \
-Wall \
-std=c++03 \
find_neighbors_python_external.cpp \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_software/Sources/Eagle_v2.4.1/src \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/include \
-Wl,-rpath,/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/lib \
-shared -fPIC -o find_neighbors_python_external.so \
-L /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/lib \
-l boost_iostreams \
-lz
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

float crop(float z, float z_max) {
  return min(z_max, max(-z_max, z));
}

void run_find_neighbors(
    int batch_index,
    int n_batch,
    float z_max,
    const char * input_file,
    const char * output) {

  const float fracR = 1;

  cout << "Computing nearest neighbors for batch " <<batch_index << " mod " << n_batch << endl;
  cout << "Cropping 'z-score' values to z_max = " << z_max << endl;

  Timer timer;

  // read numbers of individuals and regions
  FileUtils::AutoGzIfstream fin; fin.openOrExit(input_file);
  size_t N, R; fin >> N >> R; // N = 3915, R = 939
  int Nbatch = 0;
  // divide the total N by total number of batches, use the modulo as the seperator
  for (size_t n = 0; n < N; n++)
    if (n % n_batch == batch_index)
      Nbatch++;
  cout << "the actual N_batch is " << Nbatch << endl;

  string line; getline(fin, line); // throw away first header line: mu

  fin >> N >> R;
  vector <float> sigma2ratios(R);
  for (size_t r = 0; r < R; r++)
    fin >> sigma2ratios[r];

  vector <float> sigma2sort = sigma2ratios;
  sort(sigma2sort.begin(), sigma2sort.end());
  float sigma2min = sigma2sort[(int) (R*(1-fracR))];
  float sigma2max = 1000; //sigma2sort[0] * 100;
  size_t Ruse = 0, Rextreme = 0;
  for (size_t r = 0; r < R; r++) {
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
  for (size_t n = 0; n < N; n++) {
    fin >> ID >> scale;
    IDs[n] = ID;
    scales[n] = scale;
    if (n % n_batch ==batch_index) {
      int i = n/n_batch;
      for (size_t r = 0; r < R; r++) {
	fin >> z;
	zs[r*Nbatch+i] = crop(z, z_max);
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
  fin.openOrExit(input_file);
  getline(fin, line); getline(fin, line); // throw away header lines
  for (size_t n = 0; n < N; n++) {
    fin >> ID >> scale; // already stored
    for (size_t r = 0; r < R; r++) {
      fin >> z;
      if (sigma2ratios[r] < sigma2min || sigma2ratios[r] > sigma2max) continue;
      z = crop(z, z_max);
      for (size_t i = 0; i < Nbatch; i++){
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
  sprintf(buf, "%s_zMax_%g.txt.gz", output, z_max);
  FileUtils::AutoGzOfstream fout; fout.openOrExit(buf);
  fout << std::setprecision(2) << std::fixed;
  vector < pair <float, int> > distIDs(N); // this ID is INDEX of IDS vector
  for (size_t i = 0; i < Nbatch; i++) {
    int n_i = i*n_batch + batch_index;
    fout << IDs[n_i] << "\t" << scales[n_i];
    // sort and output best matches
    for (size_t n = 0; n < N; n++)
      distIDs[n] = make_pair(dists[n*Nbatch+i], n);
    distIDs[n_i].first = 1e9;
    sort(distIDs.begin(), distIDs.end());

    const size_t N_OUTPUT = 500;
    for (size_t j = 0; j < N_OUTPUT; j++) {
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

extern "C" {
    void find_neighbors(
            int batch_index,
            int n_batch,
            float z_max,
            const char * input_file,
            const char * output){
        return run_find_neighbors(
            batch_index,
            n_batch,
            z_max,
            input_file,
            output);
    }
}
