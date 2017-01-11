/*(c) Vladimir Shchur, Wellcome Trust Sanger Institute, 2015*/
#include "ARG.h"

void PerformTest(){
	int i;
//	int xa[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
	int xa[] = {1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1};
	std::vector<int> x;
//	double da[] = {0, 4, 3, 5, 2, 3, 5, 4, 1, 2};
	double da[] = {0, 5, 4, 3, 4, 5, 1, 5, 2, 4, 2, 3};
	std::vector<double> dFunc;
	int len = sizeof(xa)/sizeof(xa[0]);
    if (len != sizeof(da)/sizeof(da[0]) ){
        std::cout << "Test cannot be run: test arrays' lengths are non-equal." << std::endl;
		exit (EXIT_FAILURE);
	}
	for (i = 0; i < len; i++){
		x.push_back(xa[i]);
		dFunc.push_back(da[i]);
	}
//	int x[] = {0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1};
//	double dFunc[] = {0, 4, 4, 3, 4, 4, 2, 5, 4, 4, 5, 3, 5, 4, 5, 4, 5, 4, 5, 1, 4, 2, 5, 3, 6, 5, 3, 3, 5, 6, 3, 5, 6};
	std::cout << "Perform test." << std::endl;
    Argentum ARG (len);
	ARG.SetTree(dFunc);
    ARG.FeedSite(x);
	ARG.PrintTree();
}
    
void Help(){
    std::cout << "Help: to be developed." << std::endl;
}

int main(int argc, char *argv[]){
	std::cout << "Read VCF." << std::endl;
	char fn[] = "../aa_1k/chr20_aa_complete.vcf.gz";//INPUT FILE: edit here or add command line argument
	ReadFile(fn);
	return 1;
}
