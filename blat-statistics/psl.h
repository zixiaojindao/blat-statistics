#pragma once
#include<string>
#include<vector>
#include<iostream>
#include<cmath>
class psl 
{
public:
	int match;
	int misMatch;
	int repMatch;
	int Ns;
	int qGapCount;
	int qGapBases;
	int tGapCount;
	int tGapBases;
	std::string strand;
	std::string qName;
	int  qSize;
	int  qStart;
	int	 qEnd;
	std::string tName;	
	int tSize;
	int tStart;
	int tEnd;
	int blockCount;
	double percentIdentity;
	int score;
	std::vector<int> blockSizes;
	std::vector<int> qStarts;
	std::vector<int> tStarts;
public:
	psl(void)
	{
		percentIdentity = -1;
		score = -1;
	}
	~psl(void)
	{
	}
	friend std::istream& operator >> (std::istream &in, psl& pslNode)
	{
		in>>pslNode.match;
		in>>pslNode.misMatch;
		in>>pslNode.repMatch;
		in>>pslNode.Ns;
		in>>pslNode.qGapCount;
		in>>pslNode.qGapBases;
		in>>pslNode.tGapCount;
		in>>pslNode.tGapBases;
		in>>pslNode.strand;
		in>>pslNode.qName;
		in>>pslNode.qSize;
		in>>pslNode.qStart;
		in>>pslNode.qEnd;
		in>>pslNode.tName;
		in>>pslNode.tSize;
		in>>pslNode.tStart;
		in>>pslNode.tEnd;
		in>>pslNode.blockCount;
		int ele = 0;
		char comma;
		for(int i = 0; i < pslNode.blockCount; ++i)
		{
			in>>ele>>comma;
			pslNode.blockSizes.push_back(ele);
		}
		for(int i = 0; i < pslNode.blockCount; ++i)
		{
			in>>ele>>comma;
			pslNode.qStarts.push_back(ele);
		}
		for(int i = 0; i < pslNode.blockCount; ++i)
		{
			in>>ele>>comma;
			pslNode.tStarts.push_back(ele);
		}
		return in;
	}
	friend std::ostream& operator << (std::ostream &out, psl &pslNode)
	{
		out<<pslNode.match<<"\t";
		out<<pslNode.misMatch<<"\t";
		out<<pslNode.repMatch<<"\t";
		out<<pslNode.Ns<<"\t";
		out<<pslNode.qGapCount<<"\t";
		out<<pslNode.qGapBases<<"\t";
		out<<pslNode.tGapCount<<"\t";
		out<<pslNode.tGapBases<<"\t";
		out<<pslNode.strand<<"\t";
		out<<pslNode.qName<<"\t";
		out<<pslNode.qSize<<"\t";
		out<<pslNode.qStart<<"\t";
		out<<pslNode.qEnd<<"\t";
		out<<pslNode.tName<<"\t";
		out<<pslNode.tSize<<"\t";
		out<<pslNode.tStart<<"\t";
		out<<pslNode.tEnd<<"\t";
		out<<pslNode.blockCount<<"\t";
		for(std::vector<int>::iterator it = pslNode.blockSizes.begin(); it != pslNode.blockSizes.end(); ++it)
			out<<*it<<",";
		out<<"\t";
		for(std::vector<int>::iterator it = pslNode.qStarts.begin(); it != pslNode.qStarts.end(); ++it)
			out<<*it<<",";
		out<<"\t";
		for(std::vector<int>::iterator it = pslNode.tStarts.begin(); it != pslNode.tStarts.end(); ++it)
			out<<*it<<",";
		out<<"\t";
		out<<pslNode.percentIdentity<<"\t";
		out<<pslNode.score;
		out<<std::endl;
		return out;
	}

private :
	int round(double d)
	{
		if(d - floor(d) > 0.5)
			return ceil(d);
		else
			return floor(d);
	}

	int min(int a, int b)
	{
		return a < b ?a:b;
	}

	bool pslIsProtein(const psl *psl)
	{
		int lastBlock = psl->blockCount - 1;

		return  (((psl->strand[1] == '+' ) &&
			(psl->tEnd == psl->tStarts[lastBlock] + 3*psl->blockSizes[lastBlock])) 
			||
			((psl->strand[1] == '-') &&
			(psl->tStart == (psl->tSize-(psl->tStarts[lastBlock] + 
			3*psl->blockSizes[lastBlock])))));
	}

	int pslCalcMilliBad(bool isMrna)
	{
		int sizeMul = pslIsProtein(this) ? 3 : 1;
		int qAliSize, tAliSize, aliSize;
		int milliBad = 0;
		int sizeDif;
		int insertFactor;
		int total;

		qAliSize = sizeMul * (this->qEnd - this->qStart);
		tAliSize = this->tEnd - this->tStart;
		aliSize = min(qAliSize, tAliSize);
		if (aliSize <= 0)
			return 0;
		sizeDif = qAliSize - tAliSize;
		if (sizeDif < 0)
		{
			if (isMrna)
				sizeDif = 0;
			else
				sizeDif = -sizeDif;
		}
		insertFactor = this->qGapCount;
		if (!isMrna)
			insertFactor += this->tGapCount;

		total = (sizeMul * (this->match + this->repMatch + this->misMatch));
		if (total != 0)
			milliBad = (1000 * (this->misMatch*sizeMul + insertFactor + 
			round(3*log(1+sizeDif)))) / total;
		return milliBad;
	}	
public:
	double PercentIdentity()
	{
		percentIdentity =  100.0 - pslCalcMilliBad(true) * 0.1;
		return percentIdentity;
	}

	int pslScore()
	{
		int sizeMul = pslIsProtein(this) ? 3 : 1;

		score = sizeMul * (this->match + ( this->repMatch>>1)) -
			sizeMul * this->misMatch - this->qGapCount - this->tGapCount;
		return score;
	}

};

