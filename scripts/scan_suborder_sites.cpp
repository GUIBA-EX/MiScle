#include <thread>
#include <mutex>
#include <atomic>
#include <map>
// Independent exhaustive all-window scanner: no tail prefilter, no Python scanner imports.
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
using namespace std;
int mask(char c){switch(c){case'A':return 1;case'C':return 2;case'G':return 4;case'T':return 8;case'R':return 5;case'Y':return 10;case'S':return 6;case'W':return 9;case'K':return 12;case'M':return 3;case'B':return 14;case'D':return 13;case'H':return 11;case'V':return 7;default:return 15;}}
int comp(int x){return ((x&1)?8:0)|((x&2)?4:0)|((x&4)?2:0)|((x&8)?1:0);}
struct Rec{string id;vector<int>s;};
vector<Rec> fasta(string f){ifstream in(f);string l;vector<Rec> r;while(getline(in,l)){if(l.empty())continue;if(l[0]=='>'){istringstream ss(l.substr(1));string id;ss>>id;r.push_back({id,{}});}else for(char c:l)if(!isspace(c))r.back().s.push_back(mask(toupper(c)));}return r;}
struct Hit{int pos,mm,last3,amb;};
vector<int> primer(string q,bool reverse,bool gt){vector<int> m;for(char c:q)m.push_back(mask(c));if(reverse){std::reverse(m.begin(),m.end());for(auto &x:m)x=comp(x);}if(gt)for(auto&x:m){int old=x;if(!reverse){if(old&4)x|=1;if(old&8)x|=2;}else{if(old&2)x|=8;if(old&1)x|=4;}}return m;}
vector<Hit> scan(const vector<int>&s,const vector<int>&q,bool reverse,int limit){vector<Hit> h;for(int i=0;i+(int)q.size()<=(int)s.size();i++){int seed=0;for(int k=0;k<7;k++){int j=reverse?k:(int)q.size()-1-k;seed+=(s[i+j]&q[j])==0;}if(seed>limit)continue;int mm=0,m3=0,amb=0;for(int j=0;j<(int)q.size();j++){int x=s[i+j];amb+=x!=1&&x!=2&&x!=4&&x!=8;int mis=(x&q[j])==0;mm+=mis;if(reverse?j<3:j>=(int)q.size()-3)m3+=mis;if(mm>limit)break;}if(mm<=limit)h.push_back({i,mm,m3,amb});}return h;}
int main(int argc,char**argv){string root=argv[1];ifstream in(root+"/oligos.tsv");ofstream out(root+"/background_sites.tsv");auto records=fasta(root+"/../background.fasta");auto rev=records;for(auto &r:rev){reverse(r.s.begin(),r.s.end());for(auto &x:r.s)x=comp(x);}map<pair<string,string>,vector<string>> grouped;string name,side,q;while(in>>name>>side>>q)grouped[{side,q}].push_back(name);vector<pair<pair<string,string>,vector<string>>> jobs(grouped.begin(),grouped.end());atomic<int> next(0),done(0);mutex lock;vector<thread> workers;for(int t=0;t<4;t++)workers.emplace_back([&](){while(true){int i=next++;if(i>=(int)jobs.size())break;auto job=jobs[i];auto side=job.first.first,q=job.first.second;ostringstream block;for(int gt=0;gt<2;gt++){auto pr=primer(q,side=="R",gt);for(int ori=0;ori<2;ori++){auto &db=ori?rev:records;for(auto &rec:db)for(auto &h:scan(rec.s,pr,side=="R",gt?2:4))if(!h.amb)for(auto &name:job.second)block<<name<<'\t'<<rec.id<<'\t'<<ori<<'\t'<<gt<<'\t'<<h.pos<<'\t'<<h.mm<<'\t'<<h.last3<<'\n';}}{lock_guard<mutex> g(lock);out<<block.str();int n=++done;if(n%100==0)cerr<<n<<" / "<<jobs.size()<<" unique oligos done\n";}}});for(auto &t:workers)t.join();cerr<<"COMPLETE "<<jobs.size()<<" unique sequences\n";}