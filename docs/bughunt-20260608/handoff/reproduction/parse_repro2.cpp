#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
using namespace std;
static vector<string> splitws(string s){
  string r; for(char c: s){ if(c=='\t') r+="    "; else r+=c; }
  vector<string> o; string cur;
  for(char c: r){ if(c==' '){ if(!cur.empty()){o.push_back(cur);cur.clear();} } else cur+=c; }
  if(!cur.empty()) o.push_back(cur); return o;
}
static void run(const string&path,bool fixed){
  ifstream f(path); string str; vector<string> t;
  size_t defint_startline=0,Ndefint=0,mixrule_startline=0,Nmixrule=0,counter=0;
  while(getline(f,str)){
    if(str.find("number of defined interactions",0)!=string::npos && (!fixed||defint_startline==0)) defint_startline=counter;
    if(defint_startline>0&&counter==defint_startline+1){t=splitws(str); if(!t.empty())sscanf(t[0].c_str(),"%zu",&Ndefint);}
    if(str.find("mixing rules to overwrite",0)!=string::npos && (!fixed||mixrule_startline==0)) mixrule_startline=counter;
    if(mixrule_startline>0&&counter==mixrule_startline+1){t=splitws(str); if(!t.empty())sscanf(t[0].c_str(),"%zu",&Nmixrule);}
    counter++;
  }
  bool Use1264=false; if(!Use1264)Ndefint=0;
  bool early=(Ndefint==0&&Nmixrule==0);
  printf("  %-18s mixrule_startline=%2zu  Nmixrule=%2zu  -> %s\n",
    fixed?"NEW+FIX:":"NEW(buggy):", mixrule_startline,Nmixrule,
    early?"EARLY RETURN -> 0 overrides  [WRONG]":"process Nmixrule overrides  [CORRECT]");
}
int main(int c,char**v){ run(v[1],false); run(v[1],true); return 0; }
