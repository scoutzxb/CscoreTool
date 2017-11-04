#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;

void parseStr(string& s, char c, vector<string>& vs)
{
	vs.clear();
	for (int i=0;;++i)
	{
		string::size_type pos=s.find(c);
		if (pos==string::npos)	
		{
			vs.push_back(s);
			break;
		}
		else
		{
			string s1=s.substr(0,pos);
			vs.push_back(s1);
			string s2=s.substr(pos+1,s.size()-pos-1);
			s=s2;
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc<4)
	{
		cout<<"usage: generateEqualLengthBed <chrSizes.txt> <out.bed> <windowsize>"<<endl;
		return 0;	
	}
	ifstream fin(argv[1]);
	ofstream fout(argv[2]);
	int windowsize=atoi(argv[3]);
	
	while (!fin.eof())
	{
		string s;
		getline(fin,s);
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<2) continue;
		string chrName=vs[0];
		int chrSize=atoi(vs[1].c_str());
		int n=(chrSize+windowsize-1)/windowsize;
		for (int i=0;i<n;++i)
		{
			int l=i*windowsize+1;
			int u=(i+1)*windowsize;
			if (u>chrSize) u=chrSize;
			fout<<chrName<<"\t"<<l<<"\t"<<u<<endl;	
		}	
	}	
}