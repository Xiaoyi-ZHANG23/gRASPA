// parser_fixed.cpp  —  REFERENCE SOLUTION (answer key; do not give to the agent under test)
//
// Identical to challenge/parser_under_test.cpp except for the one-line-per-marker fix:
// the section-marker scan now keeps the FIRST occurrence instead of the LAST
// (`&& <marker>_startline == 0`). That is the real fix applied to
// src_clean/read_data.cpp's OverWrite_Mixing_Rule() (commit ea5a825).
//
// Build:  g++ -O2 -std=c++17 parser_fixed.cpp -o parser_fixed
// Run:    ./parser_fixed <path/to/force_field.def>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
using namespace std;

static vector<string> split_ws(const string& s) {
  string r;
  for (char c : s) r += (c == '\t') ? string("    ") : string(1, c);
  vector<string> out; string cur;
  for (char c : r) {
    if (c == ' ') { if (!cur.empty()) { out.push_back(cur); cur.clear(); } }
    else cur += c;
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}

static bool is_comment_or_blank(const string& s) {
  for (char c : s) { if (c == '#') return true; if (c != ' ' && c != '\t') return false; }
  return true;
}

int main(int argc, char** argv) {
  if (argc < 2) { fprintf(stderr, "usage: %s <force_field.def>\n", argv[0]); return 2; }
  const string path = argv[1];

  size_t defint_startline = 0, Ndefint = 0;
  size_t mixrule_startline = 0, Nmixrule = 0, counter = 0;
  {
    ifstream f(path);
    if (!f) { fprintf(stderr, "cannot open %s\n", path.c_str()); return 2; }
    string line; vector<string> t;
    while (getline(f, line)) {
      // FIX: first occurrence wins (do not let a later/trailing marker overwrite it).
      if (line.find("number of defined interactions", 0) != string::npos && defint_startline == 0)
        defint_startline = counter;
      if (defint_startline > 0 && counter == defint_startline + 1) {
        t = split_ws(line); if (!t.empty()) sscanf(t[0].c_str(), "%zu", &Ndefint);
      }
      // FIX: first occurrence wins.
      if (line.find("mixing rules to overwrite", 0) != string::npos && mixrule_startline == 0)
        mixrule_startline = counter;
      if (mixrule_startline > 0 && counter == mixrule_startline + 1) {
        t = split_ws(line); if (!t.empty()) sscanf(t[0].c_str(), "%zu", &Nmixrule);
      }
      counter++;
    }
  }

  const bool UseLJ1264 = false;
  if (!UseLJ1264) Ndefint = 0;

  vector<string> overrides;
  if (Nmixrule > 0) {
    ifstream f(path);
    string line; size_t c = 0, taken = 0;
    while (getline(f, line) && taken < Nmixrule) {
      if (c > mixrule_startline + 1 && !is_comment_or_blank(line)) {
        overrides.push_back(line);
        taken++;
      }
      c++;
    }
  }

  printf("Ndefint=%zu\n", Ndefint);
  printf("Nmixrule=%zu\n", Nmixrule);
  printf("overrides_applied=%zu\n", overrides.size());
  for (const string& o : overrides) {
    vector<string> tk = split_ws(o);
    for (size_t i = 0; i < tk.size(); i++)
      printf("%s%s", tk[i].c_str(), i + 1 < tk.size() ? " " : "\n");
  }
  if (Ndefint == 0 && Nmixrule == 0)
    printf("RESULT: early-return, no overrides applied\n");
  return 0;
}
