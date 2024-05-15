/*
g++ \
-O3 \
-Wall \
-std=c++03 \
realign_GRCh37.cpp \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_software/Sources/Eagle_v2.4.1/src \
-I /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/include \
-Wl,-rpath,/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/lib \
-o realign_GRCh37 \
-L /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/external_C++_library/boost_1_58_0/install/lib \
-l boost_iostreams \
-lz
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cctype>
#include <cassert>

using namespace std;

const char MIN_QUAL = ':';
const int BAD_SCORE = 9999;

const int starts[7] = {161032593, 161038148,161043694,161049238, 161054784, 161060331, 161065878};

void compute_scores(int scores[7], int pos, const string &cigar, const string &seq,
		    const string &qual, const string &ref) {
  for (int rep_bwa = 0; rep_bwa < 7; rep_bwa++)
    if (pos > starts[rep_bwa]-500 && pos < starts[rep_bwa]+500) {
      // initialize scores
      for (int rep = 0; rep < 7; rep++)
	scores[rep] = 0;
      // parse cigar string
      int seq_at = 0;
      int ref_at = pos-161032032;
      istringstream iss(cigar);
      int len; char op;
      while (iss >> len >> op) {
	if (op == 'S' || op == 'I')
	  seq_at += len;
	else if (op == 'D')
	  ref_at += len;
	else if (op == 'M') {
	  for (int i = 0; i < len; i++) {
	    if (qual[seq_at] >= MIN_QUAL) {
	      int conf = (qual[seq_at]-'!') / ((qual[seq_at]==':'||qual[seq_at]=='F') ? 1 : 2);
	      for (int rep = 0; rep < 7; rep++)
		if (seq[seq_at] != ref[ref_at + starts[rep] - starts[rep_bwa]])
		  scores[rep] += conf;
	    }
	    seq_at++;
	    ref_at++;
	  }
	}
	else break; // bad cigar
      }
      return;
    }
  // bad cigar or bad pos
  for (int rep = 0; rep < 7; rep++)
    scores[rep] = BAD_SCORE;
}

int main(void) {

  ifstream fin("/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/pipeline_columbia_local/extract_KIV_2/ref_hg19.fasta");
  string ref;
  fin >> ref; fin >> ref;
  fin.close();

  int ctr_1B_KIV3 = 0, ctr_1B_KIV2 = 0, ctr_1B_tied = 0, ctr_1A = 0;

  string qname[2]; int pos[2]; string cigar[2], seq[2], qual[2];
  int r = 0;
  int scores[2][7];
  while (cin >> qname[r&1] >> pos[r&1] >> cigar[r&1] >> seq[r&1] >> qual[r&1]) {
    if (qname[0] == qname[1]) {
      for (int k = 0; k < 2; k++)
	compute_scores(scores[k], pos[k], cigar[k], seq[k], qual[k], ref);
      //cout << qname[0];
      for (int rep = 0; rep < 7; rep++) {
	scores[0][rep] = min(BAD_SCORE, scores[0][rep] + scores[1][rep]);
	//cout << "\t" << scores[0][rep];
      }
      //cout << endl;
      int *s = scores[0];
      if (s[0]<s[1] && s[0]<s[2] && s[0]<s[3] && s[0]<s[4] && s[0]<s[5] && s[0]<s[6])
	ctr_1B_KIV3++;
      else if (s[4]<s[0] && s[4]<s[1] && s[4]<s[2] && s[4]<s[3] && s[4]<s[5] && s[4]<s[6])
	ctr_1B_KIV2++;
      else if (s[4]==s[0] && s[4]<s[1] && s[4]<s[2] && s[4]<s[3] && s[4]<s[5] && s[4]<s[6])
	ctr_1B_tied++;
      else {
	int min_1A = BAD_SCORE;
	for (int rep = 1; rep <= 6; rep++)
	  if (rep != 4)
	    min_1A = min(min_1A, s[rep]);
	if (min_1A < s[0] && min_1A < s[4])
	  ctr_1A++;
      }
      qname[0] = qname[1] = "";
    }
    r++;
  }
  cout << "\t" << ctr_1B_KIV3 << "\t" << ctr_1B_KIV2 << "\t" << ctr_1B_tied << "\t" << ctr_1A
       << endl;

  return 0;
}
