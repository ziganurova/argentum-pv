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

struct ARGEvent{
	bool type, side;
	int pos, len, nPos;
	double newD, old;
};

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
		int RecombStrat;
		std::map <std::string, std::string> RecombStrategies;

		int M;
		
		//PBWT variables
		int *y;
		int *a;
		int *b;
		int brInCol;
		
		//Unwind
		int *brTrack;
		int *recombPointer;
		std::vector<Branch> recombList;
		
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
		void CopyBranchR(int , int , int , double , int , int );
		void RecombPBWT();
		void RegisterRecomb();
		
		//read and write ARG
		std::ofstream OutARGFile;
		std::ifstream InARGFile;
		entryARG entry;
		
		ARGEvent event;
		std::vector<ARGEvent> events;
		eventPointer;

	public:
		Argentum(int size){
			int i;
			M = size;
			RecombStrat = 0;
			if (M > 0){
				y = new int [M];
				a = new int [M];
				b = new int [M];
				brTrack = new int [M];
				recombPointer = new int [M];
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
		void FeedSite(std::vector<int>&);
		void SetTree(std::vector<double>&);
		void PrintTree();
		void PrintReducedTree();
		void PrintReducedTree1();
		void DefaultControl();
		void RecombStrategy();
		int GetSize();
};

void Argetum::Reinit(size){
	this->Argentum(size);
}

class ControlPanel{
	public:
		std::map<std::string, std::string> commands;
		int lineLimit;
		int lineReport;
//		std::map<std::string, std::string> recombStrategies;
	    ControlPanel(){
	        commands["readVCF"] = "-r";
			commands["write"] = "-w";
			commands["test"] = "-test";
			commands["recombStr"] = "-rs";
			commands["lineLimit"] = "-ll";
			commands["lineReport"] = "-lr";
			LineLimit = -1;
			LineReport = 10000;
//	        recombStrategies['1'] = 'RecombTrOr';
//			recombStrategies['2'] = 'RecombPBWT';
		}
//		void ControlPanel::RecombStrategy(Argentum* , std::string , bool );
};


void UpString(char*, char** );
void ReadVcf (char* );