#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <ctime>
#include <string.h>
#include <string>
#include <sys/time.h>
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <sys/stat.h>
#include <map>
#include <algorithm>

using namespace std;

struct PackSet{
	int el[5];
	PackSet(){
		int i;
		for (i = 0; i < 5; i++)
			el[i] = 0;
	}
};

struct Recomb{
	int pos;
	int len;
	int nPos;
	double d;
};

struct entryARG{
	bool type;//true - recombination, 
	int SNP;
	Recomb rec;
};

struct Branch{
	bool type; //false for single leaf, true for internal node
	int leaf; //leaf label
};

struct Recomb1{
	Branch movedBranch;
	Branch childBranch;
	int renameBr;
};

class Argentum{
	private:

		int M;
		
		//Current site (to be set through SetSiteNumber() if necessary - by default not set)
		int siteNumber;
		
		//PBWT variables
		int *y;
		int *a;
		int *b;
		int brInCol;
		
		//Local tree
	    double *d;
	    double *d_tmp;
			
		//Reduced tree
	    int rM;
	    int rM1;
	    int rM_tmp;
	    double *rd;
	    double *rd1;
	    int *ry;
	    PackSet *rPack;
	    PackSet *rPack1;
	    PackSet *rPack_tmp;
		
		
		void ReinitSite();
		void ReduceTree();
		void ParseInterval(int , int , int );
		void AddBranches();
		void FormBranch(int , int );
		void CopyBranchR(int , int=-1 , int=-1 , double=-1.0 , int=-1 , int=-1 );
		void RecombPBWT(bool = false);
		void RegisterRecomb();
		

	public:
		//int siteNumber;
		Argentum(int size){
			int i;
			M = size;
			if (M > 0){
				y = new int [M];
				a = new int [M];
				b = new int [M];
				// brTrack = new int [M];
				// recombPointer = new int [M];
				d = new double [M];
				d_tmp = new double [M];
				rd = new double [M];
				rd1 = new double [M];
				ry = new int [M];
				rPack = new PackSet [M];
				rPack1 = new PackSet [M];
				rPack_tmp = new PackSet [M];
				for (i = 0; i < M; i++){
					a[i] = i;
					d[i] = 1.0;
					d_tmp[i] = 1.0;
				}
			}
			std::cout << "New ARGentum session initialised with " << M << " haplotypes." << std::endl;
		}
		void FeedSite(std::vector<int>&, bool = false);
		void SetTree(std::vector<double>&);
		void PrintTree();
		void PrintTreeForTest(std::vector<int>&);
		void PrintReducedTree();
		void PrintReducedTree1();
		void DefaultControl();
		void RecombStrategy();
		void SetSiteNumber(int);
		int GetSize();

};

void ReadFile(char *filename);
