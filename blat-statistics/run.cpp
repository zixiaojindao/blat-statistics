#include<iostream>
#include<bitset>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include "psl.h"
#include<algorithm>
#include<numeric>
using namespace std;

void PrintHelp()
{
	cout<<"./blat-statistics path2output.psl path2ref.fa path2query.fa"<<endl;
}

inline bool isACGT(char c)
{
	return c == 'A' || c == 'C' || c == 'G' || c == 'T'
		|| c == 'a' || c == 'c' || c == 'g' || c == 't';
}

int NX0(vector<int>& lens, double px0)
{
	int nx0Acc = 0;
	int nx0 = 0;
	for(vector<int>::iterator itc = lens.begin(); itc != lens.end(); ++itc)
	{ 
		nx0Acc += *itc;
		if(nx0Acc > px0)
		{
			nx0 = *itc;
			break;
		}
	}	
	return nx0;
}

int main(int argc, char** argv)
{
	if(argc < 4)
	{
		PrintHelp();
		exit(EXIT_FAILURE);
	}
	//open file
	string blatFilePath = argv[1];
	string refFilePath = argv[2];
	string queryFilePath = argv[3];
	ifstream fpBlat(blatFilePath);
	if(!fpBlat)
	{
		cout<<"error: open "<<blatFilePath<<endl;
		exit(EXIT_FAILURE);
	}
	ifstream fpRef(refFilePath);
	if(!fpRef)
	{
		cout<<"error: open "<<refFilePath<<endl;
		exit(EXIT_FAILURE);
	}
	ifstream fpQuery(queryFilePath);
	if(!fpQuery)
	{
		cout<<"error: open "<<queryFilePath<<endl;
		exit(EXIT_FAILURE);
	}

	string staOutPut = "stac.txt";
	ofstream fpsta(staOutPut);
	//read ref length
	string line;
	long long refLen = 0; 
	while(getline(fpRef, line))
		if(line.length() >0 && line[0]!='>')
			refLen += line.length();
	fpRef.close();

	//read psl and calculate pid and score
	string nonUseLine;
	int header = 5;
	vector<bool> ref(refLen, 0);
	string oldQName = "";
	int percent95 = 0;
	int percent50 = 0;
	int totalBlatQuery = 0;
	for(int i = 0; i < header; ++i)
	{
		getline(fpBlat, nonUseLine);
		fpsta<<nonUseLine;
		if(i == 2)
			fpsta<<"\tpid\tsocre";
		fpsta<<endl;
	}
	vector<int> contigsLen;
	double bestPid = 0;
	vector<int> blockSizes;
	vector<int> tStarts;
	string strand = "+";
	while(true)
	{
		psl pslNode;
		if(!(fpBlat>>pslNode))
			break;
		pslNode.PercentIdentity();
		pslNode.pslScore();
		fpsta<<pslNode;
		if(oldQName != pslNode.qName)
		{
			oldQName = pslNode.qName;
			totalBlatQuery += 1;
			if(bestPid > 0.95)
				percent95 += 1;
			if(bestPid > 0.5)
				percent50 += 1;
			for(int i = 0; i < tStarts.size(); ++i)
				for(int j = 0; j < blockSizes[i]; ++j)
					//if(strand == "+")
						ref[tStarts[i] - 1 + j] = true;
			bestPid = 0;
			blockSizes.clear();
			tStarts.clear();
			bestPid = pslNode.percentIdentity;
			blockSizes = pslNode.blockSizes;
			tStarts = pslNode.tStarts;
			strand = pslNode.strand;
		}
		else
		{
			if(bestPid < pslNode.percentIdentity)
			{
				bestPid = pslNode.percentIdentity;
				blockSizes = pslNode.blockSizes;
				tStarts = pslNode.tStarts;
				strand = pslNode.strand;
			}
		}
	}
	if(bestPid > 0.95)
		percent95 += 1;
	if(bestPid > 0.5)
		percent50 += 1;
	for(int i = 0; i < tStarts.size(); ++i)
		for(int j = 0; j < blockSizes[i]; ++j)
			//if(strand == "+")
				ref[tStarts[i] - 1 + j] = true;
	int sum = 0;
	for(vector<bool>::iterator it = ref.begin(); it != ref.end(); ++it)
		if(*it)
			sum += 1;

	int tl = -1;
	int totalQuery = 0;
	while(getline(fpQuery, line))
	{
		if(line.length() > 0 && line[0] == '>')
		{
			if(tl != -1)
				contigsLen.push_back(tl);
			tl = 0;
			totalQuery += 1;
		}
		else
		{
			tl += count_if(line.begin(), line.end(), isACGT);
		}
	}
	if(tl != -1)
		contigsLen.push_back(tl);

	sort(contigsLen.rbegin(), contigsLen.rend());
	int contigSum = accumulate(contigsLen.begin(), contigsLen.end(), 0);
	int n90 = NX0(contigsLen, contigSum * 0.9); 
	int n50 = NX0(contigsLen, contigSum * 0.5);
	fpsta<<"------------------------------------------------"<<endl;
	fpsta<<"ref length="<<refLen<<endl;
	fpsta<<"blat query # ="<<totalBlatQuery<<endl;
	fpsta<<"total query # ="<<totalQuery<<endl;
	fpsta<<">95% map ratio :"<<percent95<<" / "<<totalQuery<<"="<<percent95 / (double)totalQuery * 100<<"%"<<endl;
	fpsta<<">50% map ratio :"<<percent50<<" / "<<totalQuery<<"="<<percent50 / (double)totalQuery * 100<<"%"<<endl;
	fpsta<<"ref coverage ="<<sum / (double)refLen * 100<<"%"<<endl;
	fpsta<<"max query ="<<(contigsLen.size() == 0? 0:contigsLen[0])<<endl;
	fpsta<<"n50 ="<<n50<<endl;
	fpsta<<"n90 ="<<n90<<endl;
	fpBlat.close();
	fpsta.close();
	fpQuery.close();
}