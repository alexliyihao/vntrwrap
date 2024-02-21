/*
compile:
g++ \
-O3 \
-Wall \
-std=c++03 \
normalize_mosdepth_inflow_rewritten.cpp \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_software/Sources/Eagle_v2.4.1/src \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/include \
-Wl,-rpath,/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/lib \
-o normalize_mosdepth_inflow_rewritten \
-L /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/lib \
-l boost_iostreams \
-lz

running:
./normalize_mosdepth_inflow_rewritten \
the <prefix> of mosdepth batch output <prefix>_batch_<batchnumber>.txt.gz\
the repeat mask bed you are using \
choose one example mosdepth per-sample output the path of .regions.bed.gz files (C++ need that to set up the pointer size...)\
total sample size \
output path, which is a .txt.gz file > keep stdout if you need

e.g.
./normalize_mosdepth_inflow_rewritten \
/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/mosdepth/mosdepth_whicap_GRCh37 \
/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/pipeline_columbia_local/sources/WES_read_depth/repeat_mask_list.hg37.ucsc_bed \
/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/mosdepth/washei25748.BQSR.recaled.bam_GRCh37/washei25748.BQSR.recaled.bam.regions.bed.gz \
3915 \
/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/mosdepth/GRCh37_ID_scale_zdepths.txt.gz > running_log_normalize_mosdepth_GRCh37.txt

*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <set>
#include "FileUtils.cpp"
#include "StringUtils.cpp"
#include "Timer.cpp"


using namespace std;

int main(int argc, char *argv[]) {

  if (argc != 6) {
    cerr << "ERROR: 5 arguments required" << endl;
    cerr << "- arg1: prefix of mosdepth input (no more than 170 characters), used in <prefix>_batch_<batchnumber>.txt.gz" << endl;
    cerr << "- arg2: bed file path e.g. /path/to/repeat_mask_list.hg38.ucsc_bed" << endl;
    cerr << "- arg3: example input e.g. /path/to/name_regions.bed.gz" << endl;
    cerr << "- arg4: N_sample(int)" << endl;
    cerr << "- arg5: output path e.g. /path/to/ID_scale_zdepths.txt.gz" << endl;
    return 1;
  }
  const char *mosdepth_prefix = argv[1];
  const char *bed_source = argv[2];
  const char *example_output = argv[3];
  int MAX_N; sscanf(argv[4], "%d", &MAX_N);
  const char *output_path = argv[5];

  Timer timer;

  FileUtils::AutoGzIfstream fin;
  /***** DECIDE WHICH REGIONS TO EXTRACT (based on first 10 batches = ~250 indivs) *****/
  long long R = 0; // number of regions;
  vector <double> mean_depths;
  int N = 0;
  for (int batch = 0; batch < 10; batch++) {
    char buf[200];
    sprintf(buf, "%s_batch_%d.txt.gz", mosdepth_prefix, batch+1);
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
      	mean_depths.push_back(depth*0.01);
      	R++;
      }
    }
    while (fin >> ID) {
      N++;
      int depth;
      for (int r = 0; r < R; r++) {
      	fin >> depth;
      	mean_depths[r] += depth*0.01;
      }
    }
    fin.close();
    cout << "Read batch " << batch << endl;
  }

  vector <bool> extract(R);
  long long Rextract = 0;
  for (int r = 0; r < R; r++) {
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
  // we only works on Chr 6, if need more, modify here.
  string chrStr;
  std::set<int> valid_chr;
  std::set<int>::iterator it = valid_chr.begin();
  it = valid_chr.insert(it, 6);
//  for (int i = 1; i < 23; i++) {
//      it = valid_chr.insert(it, i);
//  }
  while (fin >> chrStr) {
    if (chrStr == "chrX") continue;
    int chr, bpStart, bpEnd, bpLen; string vntrStr;
    sscanf(chrStr.c_str(), "chr%d", &chr);
    //cout << "chr read" << chr << endl;
    if (valid_chr.count(chr) == 0) continue;
    //cout << "chr recorded" << chr << endl;
    numVNTRs++;
    fin >> bpStart >> bpEnd >> vntrStr >> bpLen;
    for (int kb = bpStart/1000; kb <= bpEnd/1000; kb++)
      overlapsVNTR[chr][kb] = true;
    totLen += bpEnd - bpStart;
  }
  fin.close();
  cout << "Read " << numVNTRs << " autosomal VNTRs spanning " << totLen*1e-6 << " Mb" << endl;

  // get region spans; update extract list
  int Roverlap = 0, RextractOverlap = 0;
  fin.openOrExit(example_output);
  for (int r = 0; r < R; r++) {
    int chr, bpStart, bpEnd; double dep;
    fin >> chr >> bpStart >> bpEnd >> dep;
    // this part is added for we are only working on one chr,
    // if read a non-related chromosome, skip this line
    // and reset loop to the previous counter
    if (valid_chr.count(chr) == 0) {
      r--;
      continue;
    }
    //cout << "example_output included chr " << chr << endl;
    if (overlapsVNTR[chr][bpStart/1000]) {
      Roverlap++;
      if (extract[r]) {
      	RextractOverlap++;
      	extract[r] = false;
      	Rextract--;
        //cout << "chr" << chr << " starting from " << bpStart << " dropped"<< endl;
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
  const int batch_num = ceil(MAX_N/batch_size);
  float *depths = new float[MAX_N*Rextract];
  vector <string> IDs(MAX_N);
  vector <float> scales(MAX_N);

  for (int batch = 0; batch < batch_num; batch++) {
    //if (376<=batch && batch<=423 || 1000<=batch && batch<1700) continue;
    char buf[200];
    sprintf(buf, "%s_batch_%d.txt.gz", mosdepth_prefix, batch+1);
    fin.openOrExit(buf);
    string ID;
    while (fin >> ID) {
      float *depthsRow = depths + N*Rextract;
      /*cout << "Current ID " << ID << endl;*/
      IDs[N] = ID;
      int rSub = 0;
      float mean_depthExtract = 0;
      int depth;
      for (int r = 0; r < R; r++) {
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
      for (int s = 0; s < Rextract; s++){
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
  for (int s = 0; s < Rextract; s++) {
    for (int n = 0; n < N; n++) {
      x[n] = depths[n*Rextract+s];
      if(isnan(x[n])){
        cout << "current info x: NA: n=" << n  << " s="<< s << endl;
      }
    }
    float mu = 0, s2 = 0;
    for (int n = 0; n < N; n++) mu += x[n];
    mu /= N;
    for (int n = 0; n < N; n++) s2 += (x[n]-mu)*(x[n]-mu);
    s2 /= (N-1);
    mus[s] = mu; sigma2s[s] = s2; sigma2ratios[s] = ratioMult*s2/mu;
    //float invStd = 1 / sqrtf(s2);
    float invRootMean = 1 / sqrtf(mu);
    for (int n = 0; n < N; n++) depths[n*Rextract+s] = (x[n]-mu) * invRootMean/*invStd*/;
  }
  cout << "Normalized by region (" << timer.update_time() << " sec)" << endl;

  sort(sigma2ratios.begin(), sigma2ratios.end());
  float sigma2ratioMin = sigma2ratios[(int) (0.9*Rextract)];
  vector <bool> want(Rextract); int Rwant = 0;
  for (int s = 0; s < Rextract; s++)
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
  for (int s = 0; s < Rextract; s++)
    if (want[s])
      fout << "\t" << mus[s];
  fout << endl;

  fout << N << "\t" << Rwant;
  for (int s = 0; s < Rextract; s++)
    if (want[s])
      fout << "\t" << ratioMult*sigma2s[s]/mus[s];
  fout << endl;

  fout << std::setprecision(2) << std::fixed;
  for (int n = 0; n < N; n++) {
    fout << IDs[n] << "\t" << 0.01f*scales[n];
    for (int s = 0; s < Rextract; s++)
      if (want[s])
	fout << "\t" << scale*depths[n*Rextract+s];
    fout << endl;
  }
  fout.close();
  cout << "Wrote output (" << timer.update_time() << " sec)" << endl;

  delete[] depths;
}
