/*
g++ \
-O3 \
-Wall \
-std=c++03 \
normalize_mosdepth_python_external.cpp \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_software/Sources/Eagle_v2.4.1/src \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/include \
-Wl,-rpath,/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/lib \
-shared -fPIC -o normalize_mosdepth_python_external.so \
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

void run_normalize_mosdepth(
    const char * mosdepth_prefix,
    const char * bed_source,
    const char * example_output,
    const char * output_path,
    int MAX_N) {

  Timer timer;

  FileUtils::AutoGzIfstream fin;
  /***** DECIDE WHICH REGIONS TO EXTRACT (based on first 10 batches = ~250 indivs) *****/
  unsigned long long R = 0; // number of regions;
  vector <double> mean_depths;
  size_t N = 0;
  for (size_t batch = 0; batch < 10; batch++) {
    char buf[200];
    sprintf(buf, "%s_batch_%zu.txt.gz", mosdepth_prefix, batch+1);
    fin.openOrExit(buf);
    string ID; // throw away for now
    if (R == 0) {
      N++;
      string line;
      getline(fin, line);
      istringstream iss(line);
      iss >> ID;
      int depth;
      while (iss >> depth) {
	mean_depths.push_back(depth);
	R++;
      }
    }
    while (fin >> ID) {
      N++;
      int depth;
      for (size_t r = 0; r < R; r++) {
	fin >> depth;
	mean_depths[r] += depth*0.01;
      }
    }
    fin.close();
    cout << "Read batch " << batch << endl;
  }

  vector <bool> extract(R);
  unsigned long long Rextract = 0;
  for (size_t r = 0; r < R; r++) {
    mean_depths[r] /= N;
    if (20 <= mean_depths[r] && mean_depths[r] <= 100) {
      extract[r] = true;
      Rextract++;
    }
  }

  cout << "Read " << N << " indivs (" << timer.update_time() << " sec)" << endl;
  cout << "Extracting " << Rextract << " / " << R << " regions" << endl;


  /***** EXCLUDE REGIONS OVERLAPPING VNTRs *****/

  // read VNTRs; mark 1kb regions overlapping list of repeat regions to mask
  vector < vector <bool> > overlapsVNTR(23, vector <bool> (300000));
  fin.openOrExit(bed_source);
  int numVNTRs = 0, totLen = 0;
  string chrStr;
  while (fin >> chrStr) {
    if (chrStr == "chrX") continue;
    numVNTRs++;
    size_t chr, bpStart, bpEnd, bpLen; string vntrStr;
    sscanf(chrStr.c_str(), "chr%zu", &chr);
    fin >> bpStart >> bpEnd >> vntrStr >> bpLen;
    for (size_t kb = bpStart/1000; kb <= bpEnd/1000; kb++)
      overlapsVNTR[chr][kb] = true;
    totLen += bpEnd - bpStart;
  }
  fin.close();
  cout << "Read " << numVNTRs << " autosomal VNTRs spanning " << totLen*1e-6 << " Mb" << endl;

  // get region spans; update extract list
  int Roverlap = 0, RextractOverlap = 0;
  fin.openOrExit(example_output);
  for (size_t r = 0; r < R; r++) {
    int chr, bpStart, bpEnd; double dep;
    fin >> chr >> bpStart >> bpEnd >> dep;
    if (overlapsVNTR[chr][bpStart/1000]) {
      Roverlap++;
      if (extract[r]) {
	RextractOverlap++;
	extract[r] = false;
	Rextract--;
      }
    }
    if (r==R-1) cout << "Last region: " << chr << ":" << bpStart << "-" << bpEnd << endl;
  }
  fin.close();
  cout << "Excluding " << Roverlap << " / " << R << " regions overlapping VNTRs" << endl;
  cout << "Excluded " << RextractOverlap << " in extract set; " << Rextract << " left" << endl;


  /***** extract data for all individuals; normalize each individual (across regions) *****/

  N = 0;
  const double batch_size = 25.0; /* for ceiling*/
  const size_t batch_num = ceil(MAX_N/batch_size);
  float *depths = new float[MAX_N*Rextract];
  vector <string> IDs(MAX_N);
  vector <float> scales(MAX_N);

  for (size_t batch = 0; batch < batch_num; batch++) {
    //if (376<=batch && batch<=423 || 1000<=batch && batch<1700) continue;
    char buf[200];
    sprintf(buf, "%s_batch_%zu.txt.gz", mosdepth_prefix, batch+1);
    fin.openOrExit(buf);
    string ID;
    while (fin >> ID) {
      float *depthsRow = depths + N*Rextract;
      /*cout << "Current ID " << ID << endl;*/
      IDs[N] = ID;
      int rSub = 0;
      float mean_depthExtract = 0;
      int depth;
      for (size_t r = 0; r < R; r++) {
      	fin >> depth;
        if (isnan(depth)){
          cout << "na: depth: " << r << batch << endl;
        }
      	if (extract[r]) {
      	  depthsRow[rSub++] = depth;
      	  mean_depthExtract += depth;
      	}
      }
      mean_depthExtract /= Rextract;
      if (isnan(mean_depthExtract)){
        cout << "na: mean_depthExtract: " << batch << endl;
      }
      scales[N] = mean_depthExtract;
      float invScale = 1 / mean_depthExtract;
      if (isnan(invScale)){
        cout << "na: invScale: "<< batch <<endl;
      }
      for (size_t s = 0; s < Rextract; s++){
      	depthsRow[s] *= invScale;
        if (isnan(depthsRow[s])){
          cout << "na: depthsRow[s]: " << s << batch <<endl;
        }
      }
      N++;
      /*cout << "Current N " << N << endl;*/
    }
    fin.close();
    cout << "Read batch " << batch << " (" << timer.update_time() << " sec)" << endl;
  }
  cout << "Read " << N << " indivs; normalizing by region" << endl;

  /***** normalize each region (across individuals) *****/

  float x[N];
  const float ratioMult = 100; // just for more readable output
  vector <float> mus(Rextract), sigma2s(Rextract), sigma2ratios(Rextract);
  for (size_t s = 0; s < Rextract; s++) {
    for (size_t n = 0; n < N; n++) {
      x[n] = depths[n*Rextract+s];
      if(isnan(x[n])){
        cout << "current info x: NA: n=" << n  << " s="<< s << endl;
      }
    }
    float mu = 0, s2 = 0;
    for (size_t n = 0; n < N; n++) mu += x[n];
    mu /= N;
    for (size_t n = 0; n < N; n++) s2 += (x[n]-mu)*(x[n]-mu);
    s2 /= (N-1);
    mus[s] = mu; sigma2s[s] = s2; sigma2ratios[s] = ratioMult*s2/mu;
    //float invStd = 1 / sqrtf(s2);
    float invRootMean = 1 / sqrtf(mu);
    for (size_t n = 0; n < N; n++) depths[n*Rextract+s] = (x[n]-mu) * invRootMean/*invStd*/;
  }
  cout << "Normalized by region (" << timer.update_time() << " sec)" << endl;

  sort(sigma2ratios.begin(), sigma2ratios.end());
  float sigma2ratioMin = sigma2ratios[(int) (0.9*Rextract)];
  vector <bool> want(Rextract); int Rwant = 0;
  for (size_t s = 0; s < Rextract; s++)
    if (ratioMult*sigma2s[s]/mus[s] > sigma2ratioMin) {
      want[s] = true;
      Rwant++;
    }
  cout << "Restricting to " << Rwant << " regions with sigma2ratio > " << sigma2ratioMin
       << endl;
  float sigma2ratioMedian = sigma2ratios[Rextract/2];
  cout << "Rescaling to approximate z-scores based on median sigma2ratio = " << sigma2ratioMedian
       << endl;
  float scale = 1/sqrtf(sigma2ratioMedian/ratioMult);

  /***** write output *****/

  FileUtils::AutoGzOfstream fout; fout.openOrExit(output_path);
  fout << std::setprecision(3) << std::fixed;

  fout << N << "\t" << Rwant;
  for (size_t s = 0; s < Rextract; s++)
    if (want[s])
      fout << "\t" << mus[s];
  fout << endl;

  fout << N << "\t" << Rwant;
  for (size_t s = 0; s < Rextract; s++)
    if (want[s])
      fout << "\t" << ratioMult*sigma2s[s]/mus[s];
  fout << endl;

  fout << std::setprecision(2) << std::fixed;
  for (size_t n = 0; n < N; n++) {
    fout << IDs[n] << "\t" << 0.01f*scales[n];
    for (size_t s = 0; s < Rextract; s++)
      if (want[s])
	fout << "\t" << scale*depths[n*Rextract+s];
    fout << endl;
  }
  fout.close();
  cout << "Wrote output (" << timer.update_time() << " sec)" << endl;

  delete[] depths;
}

extern "C" {
    void normalize_mosdepth(
        const char * mosdepth_prefix,
        const char * bed_source,
        const char * example_output,
        const char * output_path,
        int MAX_N){
        return run_normalize_mosdepth(
            mosdepth_prefix,
            bed_source,
            example_output,
            output_path,
            MAX_N);
    }
}
