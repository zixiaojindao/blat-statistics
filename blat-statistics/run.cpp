#include<iostream>
#include<bitset>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include "psl.h"
using namespace std;

void PrintHelp()
{
	cout<<"./blat-statistics path2output.psl path2ref.fa"<<endl;
}

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		PrintHelp();
		exit(EXIT_FAILURE);
	}
	string blatFilePath = argv[1];
	string refFilePath = argv[2];
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
	string staOutPut = "stac.txt";
	ofstream fpsta(staOutPut);

	string line;
	long long refLen = 0; 
	while(getline(fpRef, line))
		if(line.length() >0 && line[0]!='>')
			refLen += line.length();
	//cout<<refLen<<endl;
	fpRef.close();
	string nonUseLine;
	int header = 5;
	vector<bool> ref(refLen, 0);
	string oldQName = "";
	int percent95 = 0;
	int percent50 = 0;
	int totalQuery = 0;
	for(int i = 0; i < header; ++i)
	{
		getline(fpBlat, nonUseLine);
		fpsta<<nonUseLine;
		if(i == 2)
			fpsta<<"\tpid\tsocre";
		fpsta<<endl;
	}
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
			totalQuery += 1;
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
	fpsta<<"------------------------------------------------"<<endl;
	fpsta<<"ref length="<<refLen<<endl;
	fpsta<<"query seq # ="<<totalQuery<<endl;
	fpsta<<">95% map ratio ="<<percent95 / (double)totalQuery * 100<<"%"<<endl;
	fpsta<<">50% map ration ="<<percent50 / (double)totalQuery * 100<<"%"<<endl;
	fpsta<<"ref coverage ="<<sum / (double)refLen * 100<<"%"<<endl;
	fpBlat.close();
	fpsta.close();
}