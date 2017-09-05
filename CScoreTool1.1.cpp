//This program is to perform HiC compartmentalization in raw matrix file
//the process is to fit nij=Bi*Bj*hij*(1+Ci*Cj)
//inter-chr reads are nor used. i.e. results from each chrs are independent
//P denotes parallel computation

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <math.h>
#include <omp.h>
#include "twister.h"
using namespace std;

const double maxc=0.9999;
const double steplength=0.04;

typedef vector<vector<double> > Matrix ;

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

double gammaln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

void LUdecomp(const Matrix& mat, Matrix& md)
{
	int n=mat.size();
	md=mat;
	for (int j=0;j<n;++j)
	{
		for (int i=0;i<=j;++i)
		{
			double s=0;
			for (int k=0;k<i;++k)
			{
				s+=md[i][k]*md[k][j];
			}
			md[i][j]=mat[i][j]-s;
		}
		for (int i=j+1;i<n;++i)
		{
			double s=0;
			for (int k=0;k<j;++k)
			{
				s+=md[i][k]*md[k][j];
			}
			md[i][j]=(mat[i][j]-s)/md[j][j];
		}
	}
}


//|m|
double det(const Matrix& mat)
{
	if (mat.size()==0)
	{
		return 1;
		cerr<<"Not a Matrix!"<<endl;
	}
	else if (mat.size()==1)
	{
		return mat[0][0];
	} 
	else if (mat.size()==2)
	{
		return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
	}
	else if (mat.size()==3)
	{
		return mat[0][0]*mat[1][1]*mat[2][2]+mat[0][1]*mat[1][2]*mat[2][0]+mat[0][2]*mat[1][0]*mat[2][1]
			-(mat[0][2]*mat[1][1]*mat[2][0]+mat[0][1]*mat[1][0]*mat[2][2]+mat[0][0]*mat[1][2]*mat[2][1]);
	}
	else
	{
		Matrix md(mat);
		LUdecomp(mat,md);
		double s=1.0;
		for (int i=0;i<mat.size();++i)
		{
			s*=md[i][i];
		}
		return s;
	}
}

double innerprod(const vector<double>& x)
{
	int n=x.size();
	double d=0;
	for(int i=0;i<n;++i)
	{
		d+=x[i]*x[i];	
	}
	return d;
}

double innerprod(const vector<double>& x1, const vector<double>& x2)
{
	int n=x1.size();
	int n2=x2.size();
	if (n!=n2)
	{
		cerr<<"unmatched dimension in inner production!"<<endl;
		return 0;	
	}
	double d=0;
	for(int i=0;i<n;++i)
	{
		d+=x1[i]*x2[i];	
	}
	return d;
}

//m^-1
/*
Matrix inverse(const Matrix& mat)
{
	Matrix m1;
	m1=mat;
	if (mat.size()<=0) 
	{
		cerr<<"Not a Matrix!"<<endl;
	}
	Matrix md(mat);
	LUdecomp(mat,md);
	int n=mat.size();
	for (int j=0;j<n;++j)
	{
		double s=0;
		for (int k=0;k<j;++k)
		{
			s+=md[j][k]
		}
	}
}*/

//(m^-1)x

vector<double> invp(const Matrix& mat, const vector<double>& x)
{
	if (mat.size()!=x.size() || mat.size()==0)
	{
		cerr<<"Unmatched dimensions for scalar production!"<<endl;
		return x;
	}
	else if (mat.size()==1)
	{
		vector<double> d(1,x[0]/mat[0][0]);
		return d;
	}
	else
	{
		Matrix md(mat);
		LUdecomp(mat,md);
		double s=0;
		int n=mat.size();
		vector<double> c(x);
		for (int j=0;j<n;++j)
		{
			double s=0;
			for (int k=0;k<j;++k)
			{
				s+=md[j][k]*c[k];
			}
			c[j]=x[j]-s;
		}
		double sum=0;
		vector<double> d(x);
		for (int j=n-1;j>=0;--j)
		{
			double s=0;
			for (int k=j+1;k<n;++k)
			{
				s+=md[j][k]*d[k];
			}
			d[j]=(c[j]-s)/md[j][j];
		}
		return d;
	}
}

double scpro(const Matrix& mat, const vector<double>& x)
{
	if (mat.size()!=x.size() || mat.size()==0)
	{
		cerr<<"Unmatched dimensions for scalar production!"<<endl;
		return 0;
	}
	else if (mat.size()==1)
	{
		return x[0]*x[0]/mat[0][0];
	}
	else if (mat.size()==2)
	{
		return (mat[1][1]*x[0]*x[0]-(mat[1][0]+mat[0][1])*x[0]*x[1]+mat[0][0]*x[1]*x[1])/det(mat);
	}
	else if (mat.size()==3)
	{
		return ((mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])*x[0]*x[0]
			-(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0]+mat[0][1]*mat[2][2]-mat[0][2]*mat[2][1])*x[0]*x[1]
			+(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]+mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1])*x[0]*x[2]
			+(mat[0][0]*mat[2][2]-mat[0][2]*mat[2][0])*x[1]*x[1]
			-(mat[0][0]*mat[1][2]-mat[0][2]*mat[1][0]+mat[0][0]*mat[2][1]-mat[0][1]*mat[2][0])*x[1]*x[2]
			+(mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0])*x[2]*x[2])
			/det(mat);
	}
	else
	{
		Matrix md(mat);
		LUdecomp(mat,md);
		double s=0;
		int n=mat.size();
		vector<double> c(x);
		for (int j=0;j<n;++j)
		{
			double s=0;
			for (int k=0;k<j;++k)
			{
				s+=md[j][k]*c[k];
			}
			c[j]=x[j]-s;
		}
		double sum=0;
		vector<double> d(x);
		for (int j=n-1;j>=0;--j)
		{
			double s=0;
			for (int k=j+1;k<n;++k)
			{
				s+=md[j][k]*d[j];
			}
			d[j]=(c[j]-s)/md[j][j];
			sum+=x[j]*d[j];
		}
		return sum;
	}
}

void raninit(vector<double>& x)
{
	int n=x.size();
	vector<double> v(n);
	double s=0;
	for (int i=0;i<n;++i)
	{
		v[i]=-log(1-rndu());
		s+=v[i];
	}	
	for (int i=0;i<n;++i)
	{
		v[i]*=maxc/s;
		x[i]=sqrt(v[i])*(2*rndu()-1);	
	}
}


//variable-length window
class VLW
{
public:
	string chrName;
	int chrStart;
	int chrEnd;
	int vlwsnumber;
	int color; //0:R,1:O,2:Y,3:G,4:B,5:P
	bool null;
};

bool leftNonOverlap(const VLW& loc1, const VLW& loc2)
{
	if (loc1.chrName<loc2.chrName) return true;
	else if (loc1.chrName==loc2.chrName && loc1.chrEnd<=loc2.chrStart) return true;
	return false;
}

void readVLWs(const string& ifname, vector<VLW>& vlws)
{
	vlws.clear();
	ifstream fin(ifname.c_str());
	while (!fin.eof())
	{
		string s;
		getline(fin,s);
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<4) continue;
		VLW vlw1;
		vlw1.chrName=vs[0];
		vlw1.chrStart=atoi(vs[1].c_str());
		vlw1.chrEnd=atoi(vs[2].c_str());
		if (vs[3]=="RED") vlw1.color=0;
		if (vs[3]=="ORA") vlw1.color=1;
		if (vs[3]=="YEL") vlw1.color=2;
		if (vs[3]=="GRE") vlw1.color=3;
		if (vs[3]=="BLU") vlw1.color=4;
		if (vs[3]=="PUR") vlw1.color=5;
		vlws.push_back(vlw1);
	}
}

int readVLWindow(const string& ifname, vector<VLW>& vlws)
{
	vlws.clear();
	ifstream fin(ifname.c_str());
	if (!fin.is_open())
	{
		cout<<"Bad window File Name"<<endl;
		return -1;	
	}
	string s;
	while (!fin.eof())
	{
		getline(fin,s);
		if (s==""||s.substr(0,3)!="chr") continue;
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<3) continue;
		VLW vlw1;
		vlw1.chrName=vs[0];
		vlw1.chrStart=atoi(vs[1].c_str())-1;
		vlw1.chrEnd=atoi(vs[2].c_str());
		vlw1.null=true;
		vlws.push_back(vlw1);
	}
}

int findInVLWs(const VLW& vlw1, const vector<VLW>& vlws, int start, int end)
{
	if (start>=end) return -1;
	int j=(start+end)/2;
	if (vlw1.chrName==vlws[j].chrName && vlw1.chrStart>=vlws[j].chrStart && vlw1.chrEnd<=vlws[j].chrEnd) return j;
	else if (leftNonOverlap(vlw1,vlws[j])) return findInVLWs(vlw1,vlws,start,j);
	else if (leftNonOverlap(vlws[j],vlw1)) return findInVLWs(vlw1,vlws,j+1,end);
	else return -1;
}

int vlwdis(const vector<VLW>& vlws, int i, int j)
{
	if (vlws[i].chrName!=vlws[j].chrName) return -1;
	if (i>j) swap(i,j);
	return (vlws[j].chrStart+vlws[j].chrEnd-vlws[i].chrStart-vlws[i].chrEnd)/2;
}

void AssignVLWs(vector<VLW>& vlws, const vector<VLW>& vlws1) //assign color and VLW number for each bin in vlws
{
	for (int i=0;i<vlws.size();++i)
	{
		int j=findInVLWs(vlws[i],vlws1,0,vlws1.size());
		vlws[i].vlwsnumber=j;
		if (j>=0) vlws[i].color=vlws1[j].color;
		else vlws[i].color=-1;
	}
}
int findInvlws(int j,const vector<int>& idxvlws,int start,int end)
{
	if (start>=end) return -1;
	int k=(start+end)/2;
	if (j==idxvlws[k]) return k;
	else if (j>idxvlws[k]) return findInvlws(j,idxvlws,k+1,end);
	else return findInvlws(j,idxvlws,start,k);
}

class Analysis
{
public:
	double d0;
	int dmin;
	string chrana;
//	double steplength;
	int anaTp;
	int nsession;
	vector<VLW> regions;
	map<string,pair<int,int> > chrIdxRange;
	vector<map<int,double> > interactions;
	vector<double> bias;
	vector<double> cscore;
	double maxdevc;
	vector<double> hh;
//	vector<double> h;
	vector<double> Nd;
	vector<double> Nr;
	double Nt;
	double delta;
	vector<double> sumbh;
	vector<double> sumbch;
	map<string,double> sumbchr;
	double sumb;
	double likelihood;
	int readRegions(const string& fileName);
	int readInteractions(const string& fileName);
	int readMatrix(const string& fileName);
	int initBias();
	int initCscore();
	int updatehh();
	int updateBias();
	int updateCscore();
	int updateLikelihood();
};

int Analysis::readRegions(const string& fileName)
{
	readVLWindow(fileName,regions);
	cout<<regions.size()<<endl;
	cout<<"read VLWs File"<<endl;
	sort(regions.begin(),regions.end(),leftNonOverlap);
	cout<<"VLWs sorted"<<endl;
	d0=regions[0].chrEnd-regions[0].chrStart;
	chrIdxRange.clear();
	for (int i=0;i<regions.size();++i)
	{
		string cn=regions[i].chrName;
		map<string,pair<int,int> >::iterator itr1=chrIdxRange.find(cn);
		if (itr1==chrIdxRange.end()) 
		{
			chrIdxRange.insert(make_pair(cn,make_pair(i,i+1)));
		}
		else itr1->second.second=i+1;
	}
	cout<<"chr ranges assigned"<<endl;
	return 0;
}

int Analysis::readInteractions(const string& fileName)
{
	ifstream fin(fileName.c_str());
	if (!fin.is_open())
	{
		cout<<"Bad interaction File Name"<<endl;
		return -1;	
	}
	int count=0;
	interactions.clear();
	Nd.clear();
	Nd.resize(132,0);
	hh.clear();
	hh.resize(132,0);
	Nr.clear();
	Nr.resize(regions.size(),0);
	interactions.clear();
	interactions.resize(regions.size());
	while (!fin.eof())
	{
		string s;
		getline(fin,s);
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<6) continue;
		++count;
		if (count%1000000==0) cout<<count<<endl;
		VLW vlw1;
		VLW vlw2;
		vlw1.chrName=vs[1];
		vlw1.chrStart=atoi(vs[2].c_str());
		vlw1.chrEnd=vlw1.chrStart+1;
		vlw2.chrName=vs[4];
		vlw2.chrStart=atoi(vs[5].c_str());
		vlw2.chrEnd=vlw2.chrStart+1;	
		int i1=findInVLWs(vlw1,regions,0,regions.size());
		int i2=findInVLWs(vlw2,regions,0,regions.size());
		if (i1<0 ||i2<0) continue;
		if (regions[i1].chrName!=regions[i2].chrName) continue;
		if (chrana!="" && (regions[i1].chrName!=chrana ||regions[i2].chrName!=chrana)) continue;
		double d=d0*abs(i2-i1);
		int dd=(log10(d+0.01)-3)/steplength;
		if (dd<dmin) continue;
		if (dd>=hh.size()) dd=hh.size()-1;
		map<int,double>::iterator itr1=interactions[i1].find(i2);
		if (itr1==interactions[i1].end()) interactions[i1].insert(make_pair(i2,1));
		else itr1->second++;
		map<int,double>::iterator itr2=interactions[i2].find(i1);
		if (itr2==interactions[i2].end()) interactions[i2].insert(make_pair(i1,1));
		else itr2->second++;
		regions[i1].null=false;
		regions[i2].null=false;
		Nd[dd]++;
		Nr[i1]++;
		Nr[i2]++;
	}
	return 0;
}

int Analysis::readMatrix(const string& fileName)
{
	ifstream fin(fileName.c_str());
	int count=0;
	interactions.clear();
	Nd.clear();
	Nd.resize(132,0);
	hh.clear();
	hh.resize(132,0);
	Nr.clear();
	Nr.resize(regions.size(),0);
	interactions.clear();
	interactions.resize(regions.size());
	while (!fin.eof())
	{
		string s;
		getline(fin,s);
		vector<string> vs;
		parseStr(s,'\t',vs);
		if (vs.size()<4) continue;
		++count;
		if (count%1000000==0) cout<<count<<endl;
		VLW vlw1;
		VLW vlw2;
		vlw1.chrName=vs[0];
		vlw1.chrStart=atoi(vs[1].c_str());
		vlw1.chrEnd=vlw1.chrStart+1;
		vlw2.chrName=vs[0];
		vlw2.chrStart=atoi(vs[2].c_str());
		vlw2.chrEnd=vlw2.chrStart+1;
		double x=atof(vs[3].c_str());
		int i1=findInVLWs(vlw1,regions,0,regions.size());
		int i2=findInVLWs(vlw2,regions,0,regions.size());
		if (i1<0 ||i2<0) continue;
		if (regions[i1].chrName!=regions[i2].chrName) continue;
		if (chrana!="" && (regions[i1].chrName!=chrana ||regions[i2].chrName!=chrana)) continue;
		double d=d0*abs(i2-i1);
		int dd=(log10(d+0.01)-3)/steplength;
		if (dd<dmin) continue;
		if (dd>=hh.size()) dd=hh.size()-1;
		map<int,double>::iterator itr1=interactions[i1].find(i2);
		if (itr1==interactions[i1].end()) interactions[i1].insert(make_pair(i2,1));
		else itr1->second+=x;
		map<int,double>::iterator itr2=interactions[i2].find(i1);
		if (itr2==interactions[i2].end()) interactions[i2].insert(make_pair(i1,1));
		else itr2->second+=x;
		regions[i1].null=false;
		regions[i2].null=false;
		Nd[dd]++;
		Nr[i1]++;
		Nr[i2]++;
	}
	return 0;
}

int Analysis::initBias()
{
	bias.clear();
	bias.resize(regions.size(),0);
	for (int i=0;i<bias.size();++i)
	{
		if (!regions[i].null) bias[i]=1;	
	}
	return 0;	
}

int Analysis::initCscore()
{
	cscore.clear();
	cscore.resize(regions.size(),0);
	for (int i=0;i<cscore.size();++i)
	{
		if (!regions[i].null) cscore[i]=(2.0*rndu()-1)*maxc;	
	}
	maxdevc=0;
	return 0;		
}

int Analysis::updatehh()
{
	vector<double> sumcc(hh.size(),0);
	vector<vector<double> > sumccpl(nsession,vector<double>(hh.size(),0));
	for (int k=0;k<nsession;++k)
	{
		for (int i=k;i<regions.size();i+=nsession)
		{
			if (regions[i].null) continue;
			for(int j=i+1;j<regions.size();++j)
			{
				if (regions[j].chrName!=regions[i].chrName) break;	
				if (regions[j].null) continue;
				double d=d0*abs(j-i);
				int dd=(log10(d+0.01)-3)/steplength;
				if (dd<dmin) continue;
				if (dd>=hh.size()) dd=hh.size()-1;
				sumccpl[k][dd]+=bias[i]*bias[j]*(1+cscore[i]*cscore[j]);
			}	
		}
	}
	for (int k=0;k<nsession;++k)
	{
		for (int dd=dmin;dd<sumcc.size();++dd)
		{
			sumcc[dd]+=sumccpl[k][dd];	
		}	
	}
//	hh[dmin]=1;
	for (int i=dmin;i<hh.size();++i)
	{
		if (Nd[i]>0) hh[i]=Nd[i]/sumcc[i];
	}
	sumbh.clear();
	sumbh.resize(regions.size(),0);
	sumbch.clear();
	sumbch.resize(regions.size(),0);
	#pragma omp parallel for
	for (int i=0;i<regions.size();++i)
	{
		if (regions[i].null) continue;
		string cn=regions[i].chrName;
		int r1=chrIdxRange[cn].first;
		int r2=chrIdxRange[cn].second;
		double F=0,G=0;
		for (int j=r1;j<r2;++j)
		{
			double d=abs(j-i)*d0;
			int dd=(log10((double)d+0.01)-3)/steplength;
			if (dd<dmin) continue;
			if (dd>=hh.size()) dd=hh.size()-1;
			double u=bias[j]*hh[dd];
			F+=u;
			G+=u*cscore[j];
		}
		sumbh[i]=F;
		sumbch[i]=G;
	}
	double L=0;
	for (int i=0;i<regions.size();++i)
	{
		if (regions[i].null) continue;
		L+=Nr[i]*log(bias[i]);
		for (map<int,double>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
		{
			int j=itr->first;
			if (j>i) break;
			double nij=itr->second;
			L+=nij*log(1+cscore[i]*cscore[j]);
		}
	}
	for (int i=dmin;i<hh.size();++i)
	{
		if (hh[i]>0) L+=Nd[i]*log(hh[i]);
	}
//	L-=sumcc[dmin];
	cout<<"L="<<L<<"\tmaxdevc="<<maxdevc<<endl;
	if (L<likelihood+1) return 1;
	likelihood=L;
	return 0;	
}

int Analysis::updateLikelihood()
{
	
}

int Analysis::updateBias()
{
	for (int i=0;i<regions.size();++i)
	{
		if (regions[i].null) continue;
		string cn=regions[i].chrName;
		int r1=chrIdxRange[cn].first;
		int r2=chrIdxRange[cn].second;
		double x=Nr[i]/(sumbh[i]+cscore[i]*sumbch[i]);
		double bb=x-bias[i];
		#pragma omp parallel for
		for (int j=r1;j<r2;++j)
		{
			double d=d0*abs(j-i);
			int dd=(log10((double)d+0.01)-3)/steplength;
			if (dd<dmin) continue;
			if (dd>=hh.size()) dd=hh.size()-1;
			double u=bb*hh[dd];
			sumbh[j]+=u;
			sumbch[j]+=u*cscore[i];
		}
		bias[i]=x;	
	}
}

int Analysis::updateCscore()
{	
	maxdevc=0;
	for (int i=0;i<regions.size();++i)
	{
		if (regions[i].null) continue;
		string cn=regions[i].chrName;
		int r1=chrIdxRange[cn].first;
		int r2=chrIdxRange[cn].second;
		double F=bias[i]*sumbch[i];
		double x=cscore[i];
		for (int nlp1=0;nlp1<100;++nlp1)
		{
			double fl=0;
			double gl=0;
			int j;
			double u,y,z,nij;
			for (map<int,double>::iterator itr=interactions[i].begin();itr!=interactions[i].end();++itr)
			{
				j=itr->first;
				nij=itr->second;
				u=cscore[j]/(1+x*cscore[j]);
				y=nij*u;
				z=y*u;
				fl+=y;
				gl+=z;
			}
			double v=fl-F;
			if (x==maxc && v>0) break;
			if (x==-maxc && v<0) break;
			if (fabs(v)<1.0e-8) break;
			x+=v/gl;
			if (x>maxc) x=maxc;
			if (x<-maxc) x=-maxc;
		}
		double cc=x-cscore[i];//bias[i]*(x-cscore[i]);
		if (fabs(cc)>maxdevc) maxdevc=fabs(bias[i]*cc);
//		sumbc+=cc;
//		bc[cn]+=cc;
		#pragma omp parallel for
		for (int j=r1;j<r2;++j)
		{
			double d=d0*abs(j-i);
			int dd=(log10((double)d+0.01)-3)/steplength;
			if (dd<dmin) continue;
			if (dd>=hh.size()) dd=hh.size()-1;
			sumbch[j]+=cc*bias[i]*hh[dd];
		}
		cscore[i]=x;
	}
}

int main(int argc, char* argv[])
{
	if (argc<6)
	{
		cout<<"Usage: CscoreTool1.1 <windows.bed> <input.summary> <OutputPrefix> <session> <minDis> [chrName]"<<endl;
		return 0;	
	}
	cout<<setprecision(17);
	SetSeed(CreateSeed());
	string chrs="";
	if (argc>6)
	chrs=argv[6];
	
	Analysis analysis;
	analysis.chrana=chrs;
	analysis.nsession=atoi(argv[4]);
	int minDis=atoi(argv[5]);
	analysis.dmin=(log10((double)minDis+0.01)-3)/steplength;
	cout<<"dmin="<<analysis.dmin<<endl;

	string outputprefix=argv[3];
	string windowFileName=argv[1];
	int f1=analysis.readRegions(windowFileName);
	if (f1<0) return 0;
	string inputFileName=argv[2];
	int f2=analysis.readInteractions(inputFileName);
	if (f2<0) return 0;
	int count=0;
	analysis.initBias();
	analysis.initCscore();
	analysis.likelihood=-1.0e308;
	omp_set_num_threads(analysis.nsession);
	for (int i=0;i<1000;++i)
	{
		int flg=analysis.updatehh();
		if (flg>0) break;		
		analysis.updateBias();
	//	analysis.updateLikelihood();
		
		analysis.updateCscore();	
	}
	
	string ofname0=outputprefix+"bias.txt";
	ofstream fbias(ofname0.c_str());
	for (int i=0;i<analysis.bias.size();++i)
	{
		fbias<<i;//idxvlws[i];
		fbias<<"\t"<<analysis.bias[i];
		fbias<<endl;
	}
	
	string ofname1=outputprefix+"_cscore.txt";
	ofstream fegc(ofname1.c_str());
	for (int i=0;i<analysis.cscore.size();++i)
	{
		fegc<<i;//idxvlws[i];
		fegc<<"\t"<<analysis.cscore[i];
		fegc<<endl;
	}
	string ofname2=outputprefix+"_hh.txt";
	ofstream fhh(ofname2.c_str());
	for (int i=0;i<analysis.hh.size();++i)
	{
		if (analysis.hh[i]<=0) continue;
		fhh<<i;//idxvlws[i];
		fhh<<"\t"<<analysis.hh[i];
		fhh<<endl;
	}
	
	string ofname3=outputprefix+"_cscore.bedgraph";
	ofstream fbdg(ofname3.c_str());
	fbdg<<"track type=\"bedGraph\" name=\""<<ofname3<<"\""<<endl;
	for (int i=0;i<analysis.cscore.size();++i)
	{
		if (analysis.bias[i]>0.2 && analysis.bias[i]<5)
		{
			fbdg<<analysis.regions[i].chrName<<"\t"<<analysis.regions[i].chrStart<<"\t"<<analysis.regions[i].chrEnd;
			fbdg<<"\t"<<analysis.cscore[i];
			fbdg<<endl;
		}
	}
}