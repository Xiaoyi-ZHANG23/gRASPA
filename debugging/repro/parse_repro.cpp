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
int main(int argc,char**argv){
  string path=argv[1];
  { // NEW first-pass (186e4d3 .. HEAD) -- exact copy of control flow
    ifstream f(path); string str; vector<string> t;
    size_t defint_startline=0,Ndefint=0,mixrule_startline=0,Nmixrule=0,counter=0;
    while(getline(f,str)){
      if(str.find("number of defined interactions",0)!=string::npos) defint_startline=counter;
      if(defint_startline>0&&counter==defint_startline+1){t=splitws(str); if(!t.empty())sscanf(t[0].c_str(),"%zu",&Ndefint);}
      if(str.find("mixing rules to overwrite",0)!=string::npos) mixrule_startline=counter;
      if(mixrule_startline>0&&counter==mixrule_startline+1){t=splitws(str); if(!t.empty())sscanf(t[0].c_str(),"%zu",&Nmixrule);}
      counter++;
    }
    bool Use1264=false; if(!Use1264)Ndefint=0;
    bool early=(Ndefint==0&&Nmixrule==0);
    printf("  NEW (186e4d3..HEAD): mixrule_startline=%zu  Nmixrule=%zu  -> %s\n",
      mixrule_startline,Nmixrule, early?"EARLY RETURN -> 0 overrides applied  [BUG]":"process N overrides");
  }
  { // OLD first-pass (pre-12-6-4) -- has the break
    ifstream f(path); string str; vector<string> t;
    size_t startline=0,Noverwrite=0,counter=0;
    while(getline(f,str)){
      if(str.find("mixing rules to overwrite",0)!=string::npos) startline=counter;
      if(startline>0&&counter==startline+1){t=splitws(str); if(!t.empty())sscanf(t[0].c_str(),"%zu",&Noverwrite); break;}
      counter++;
    }
    printf("  OLD (pre-12-6-4)   : startline=%zu  Noverwrite=%zu  -> %s\n",
      startline,Noverwrite, (Noverwrite==0)?"EARLY RETURN -> 0 overrides applied":"process Noverwrite overrides  [CORRECT]");
  }
  return 0;
}
