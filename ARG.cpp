#include "ARG.h"
/*    #Recombination list
    self.recombList = []
    self.localRecomb = []
    for i in range(self.M):
        tmp = [0, 0, 0]
        self.localRecomb.append(tmp)
        del tmp
        
    self.brStack = []
    self.brP = 0
    for i in range(self.M):
        tmp = [0, 0, 0, 0]
        self.brStack.append(tmp)
        del tmp
}
*/

int Argentum::GetSize(){
	return M;
}

void Argentum::ReinitSite(){
	rM = 0;
	rM1 = 0;
	brInCol = 0;
//	brP = 0;
}
    
void Argentum::FeedSite(std::vector<int>& x){
	int i;
	if (x.size() != M){
		std::cerr << "Incorrect number of alleles: contains " << x.size() << ", expected " << M << '.' << std::endl;
		exit (EXIT_FAILURE);
	}
	for (i = 0; i < M; i++)
		y[a[i]] = x[i];
//	PrintTree();
	ReinitSite();
	ReduceTree();
//	PrintReducedTree();
	AddBranches();
//	PrintReducedTree1();
	RecombPBWT();
}

void Argentum::ReduceTree(){
    int L = 0, h, i;
	double H = -1.0;
    for (i = 0; i < M; i++){
        if (d[i] < H or H == -1.0){
            H = d[i];
            h = i;
		}
        if (i == M-1 or y[i] != y[i+1]){
            if (i == M-1 or d[i+1] < H)
                h = i + 1;
            ParseInterval(L, i+1, h);
            H = -1;
            L = i+1;
		}
	}
}
			
void Argentum::ParseInterval(int L, int R, int h){
	int allele = y[L], i, j;
	int hBr, HBr = -1;
	double DL, DR;
	char c;
    rM_tmp = 0;
    if (L == R-1){
        rd[rM] = d[L];
        rPack[rM].el[0] = L;
        rPack[rM].el[1] = R;
        rPack[rM].el[2] = allele;
		rPack[rM].el[3] = -1;
        rM ++;
        return ;
	}
    i = L;
    while (i < h){
        j = i + 1;
        DL = d[i];
		hBr = -1;
        while (j < M and d[j] > DL){
			if (hBr == -1 or HBr > d[j]){
				hBr = j;
				HBr = d[j];
			}
            j += 1;
		}
        rd[rM] = DL;
        rPack[rM].el[0] = i;
        rPack[rM].el[1] = j;
        rPack[rM].el[2] = allele;
		rPack[rM].el[3] = hBr;
        rM += 1;
        i = j;
	}
    i = R;
    while (i > h){
        j = i - 1;
        DR = d[i];
		hBr = -1;
        while (d[j] > DR){
			if (hBr == -1 or HBr > d[j]){
				hBr = j;
				HBr = d[j];
			}
    		j--;
		}
        rPack_tmp[rM_tmp].el[0] = j;
        rPack_tmp[rM_tmp].el[1] = i;
		rPack_tmp[rM_tmp].el[3] = hBr;
        rM_tmp++;
        i = j;
	}
    for (i = rM_tmp-1; i > -1; i--){
        L = rPack_tmp[i].el[0];
        R = rPack_tmp[i].el[1];
		hBr = rPack_tmp[i].el[3];
        rd[rM] = d[L];
        rPack[rM].el[0] = L;
        rPack[rM].el[1] = R;
        rPack[rM].el[2] = allele;
		rPack[rM].el[3] = hBr;
        rM++;
	}
}

void Argentum::AddBranches(){
    int i = 0, j, k, H, v = 0;
    while (i < rM){
        if (i == rM - 1){
			CopyBranchR(i);
            break;
		}
        if (rd[i] > rd[i+1]){
			CopyBranchR(i);
            i++;
            continue;
		}
        j = i;
        H = rd[i+1];
		v = 0;
        if (rPack[j].el[2] == 1)
            v++;
        while (j < rM - 1 and rd[j+1] == H){
            j++;
            if (rPack[j].el[2] == 1)
                v++;
		}
        if (j < rM - 1 and rd[j+1] > H){
            if (rPack[j].el[2] == 1)
                v--;
            j--;
		}
        if (v > 1){
			FormBranch(i, j+1);
            i = j + 1;
		}
        else{
            for (k = 0; k <j-i+1; k++)
                CopyBranchR(i+k);
            i = j+1;
		}
	}
}
					

void Argentum::CopyBranchR(int id, int L , int R , double newD , int allele, int h ){
	int i;
	int H;
    if (L == -1)
        L = rPack[id].el[0];
    if (R == -1)
        R = rPack[id].el[1];
    if (allele == -1)
        allele = rPack[id].el[2];
    if (newD == -1)
        newD = rd[id];
    if (allele == 1)
        brInCol++;
	if (h == -1){
		for (i = L+1; i < R; i++)
			if (h == -1 or d[i] < H){
				h = i;
				H = d[i];
			}
	}
    rPack1[rM1].el[0] = L;
    rPack1[rM1].el[1] = R;
    rPack1[rM1].el[2] = allele;
	rPack1[rM1].el[3] = h;
    rd1[rM1] = newD;
    rM1++;
}

             
void Argentum::FormBranch(int L, int R){ //this relies on the fact that there is at least one 0-branch and one 1-branch, otherwise the interval would be reduced at the previous step
    int p = 0, q = 0, D = -1, u;
	int i, j;
	int oneBr;
	int lastZero = L;
	double DL = rd[L];
/*	for (i = L; i < R; i++)
		if (rPack[i].el[2] == 1)
			numOnes += rPack[i].el[1] - rPack[i].el[0];*/
	for (i = L; i < R; i++){
		if (rPack[i].el[2] == 1){
			if (p == 0){
				oneBr = rM1;
				u = rPack[i].el[1];
				CopyBranchR(i);
				for (j = rPack[i].el[0]+1; j < rPack[i].el[1]; j++)
					d[j]++;
			}
			else{
				for (j = rPack[i].el[0]; j < rPack[i].el[1]; j++){
					d[u] = d[j]+1;
					a[u] = a[j];
					u++;
				}
			}
			p += rPack[i].el[1] - rPack[i].el[0];
		}
		else{
			if (i < R - 1 && rPack[i+1].el[2] == 1)
				lastZero = i;
			if (p == 0)
				CopyBranchR(i);
			else{
				CopyBranchR(i, rPack[i].el[0]-p, rPack[i].el[1]-p);
				for (j = rPack[i].el[0]; j < rPack[i].el[1]; j++){
					d_tmp[q] = d[j];
					b[q] = a[j];
					q++;
				}
			}
		}
	}
	rPack1[oneBr].el[1] = rPack1[i].el[0] + p;
	for (i = oneBr+1; i < rM1; i++){
	    rPack1[i].el[0] += p;
	    rPack1[i].el[1] += p;
	}
	if (oneBr+1 < rM1){
		j = 0;
    	for(i = rPack1[oneBr+1].el[0]; i < rPack1[lastZero].el[1]; i++){
			d[i] = d_tmp[j];
			a[i] = b[j];
			j++;
		}
	}
}

void Argentum::RecombPBWT(){
    int size = 0, brId;
    int num1 = 0; //total number of ones
    int num2 = 0; //number of ones after the stable position
    int u = 0, p, i, j, H;
    double DL = 0.0, D1;
    for (i = 0; i < rM1; i++){
        if (rPack1[i].el[2] == 0)
            continue;
        num1 += rPack1[i].el[1] - rPack1[i].el[0];
        num2 += rPack1[i].el[1] - rPack1[i].el[0];
        if (size < rPack1[i].el[1] - rPack1[i].el[0]){
            size = rPack1[i].el[1] - rPack1[i].el[0];
            brId = i;
            num2 = rPack1[i].el[1] - rPack1[i].el[0];
		}
	}
    p = rPack1[brId].el[0] - (num1 - num2);
    D1 = rd1[brId];
    for (i = 0; i < rM1; i++){
        if (rPack1[i].el[2] == 1){
            DL = std::min(rd1[i], DL);
            if (i == brId){
                u = rPack1[i].el[0] + num2;
				DL = -1;
			}
            d_tmp[p] = d[rPack1[brId].el[3]] + 1.0;
            b[p] = a[ rPack1[i].el[0] ];
            p++;
			if (rPack1[i].el[3] == -1)
				H = 0;
			else
				H = d[rPack1[i].el[3]];
            for (j = rPack1[i].el[0]+1; j < rPack1[i].el[1]; j++){
				if (i == brId)
					d_tmp[p] = d[j] + 1.0;
				else
					d_tmp[p] = d[j] + 2.0 + d[rPack1[brId].el[3]] - H;
                b[p] = a[j];
                p++;
			}
		}
        else{
            b[u] = a[ rPack1[i].el[0] ];
			if (DL == -1)
				d_tmp[u] = d[ rPack1[i].el[0] ];
			else
				d_tmp[u] = std::min(d[ rPack1[i].el[0] ], DL);
            u++;
            if (i < rM1 - 1)
                DL = rd1[i+1];
            for (j = rPack1[i].el[0]+1; j < rPack1[i].el[1]; j++){
                d_tmp[u] = d[j];
                b[u] = a[j];
                u++;
			}
		}
	}
    d_tmp[ rPack1[brId].el[0] - (num1 - num2) ] = D1;
/*	int swap1 = *a;
	*a = *b;
	*b = swap1;
	double swap2 = *d;
	*d = *d_tmp;
	*d_tmp = swap2;*/
	std::swap(a, b);
	std::swap(d, d_tmp);
}

void Argentum::SetTree(std::vector<double>& dFunc){
	int i;
	int arSize = dFunc.size();
	if (arSize != M){
		std::cerr << "Unable to set tree: incorrect array size: M = " << M << ", array size = " << arSize << std::endl;
		exit (EXIT_FAILURE);
	}
	if (*std::min_element(dFunc.begin(),dFunc.end()) < 0 ){
		std::cerr << "Unable to set tree: negative value of similarity function." << std::endl;
		exit (EXIT_FAILURE);
	}
	for (i = 0; i < M; i++)
		d[i] = dFunc[i];
	std::clog << "New tree is set." << std::endl;
}	

/*int Argentum::BranchRoot(int pos, int len){//Checks if interval [pos, pos+len) is a valid branch. If yes returns it root node, otherwise -1.
	if (pos < 0 or pos+len > M)
		return -1;
	int h = (int)(std::min_element(d+pos+1,d+pos+len) - d);
	if (h < d[pos] or (len+pos < M and h < d[len+pos]))
		return -1;
	else
		return h;
}

void Argentum::ApplyRecombination(Recomb* recomb, bool correct = false){
	int h = -1;
	double H;
	int i;
	int pos1, pos2, pos3;
	double d1, d2, d3;
	if (recomb.nPos < 0 or recomb.nPos > M - 1){
		std::cout << "Invalid new position." << std::endl;
		exit(EXIT_FAILURE);
	}
	if (recomb.d < d[nPos]){
		std::cout << "Invalid similarity value." << std::endl;
		exit(EXIT_FAILURE);
	}
	h = BranchRoot(recomb.pos, recomb.len);
	if (h == -1){
		std::cout << "The interval is not a valid branch." << std::endl;
		exit(EXIT_FAILURE);
	}
	H = d[h];
	if (h < recomb.d){
		if (!correct){
			std::cout << "Invalid recombination." << std::endl;
			exit(EXIT_FAILURE);
		}
		else{
			int cor = recomb.d - h + 1;
			for (i = pos+1; i < pos + len; i++)
				d[i] += cor;
		}
	}
	if (recomb.pos < recomb.nPos){
		d1 = min_element(d[recomb.pos], d[recomb.pos+recomb.len]);
		d2 = d[recomb.nPos]:;
		d3 = recomb.d;
		
		memcpy(b+recomb.nPos-recomb.len, a+recomb.pos, sizeof(int)*recomb.len);
		memcpy(d_tmp+recomb.nPos-recomb.len, d+recomb.pos, sizeof(double)*recomb.len);

		memcpy(b+recomb.pos, a+recomb.pos+recomb.len, sizeof(int)*(recomb.nPos - recomb.pos-recomb.len) );
		memcpy(d_tmp+recomb.pos, d+recomb.pos+recomb.len, sizeof(double)*(recomb.nPos - recomb.pos-recomb.len) );
		
		memcpy(a+recomb.pos, b+recomb.pos, sizeof(int)*(recomb.nPos - recomb.pos) );
		memcpy(d+recomb.pos, d_tmp+recomb.pos, sizeof(double)*(recomb.nPos - recomb.pos) );
		
		d[recomb.pos] = d1;
		d[recomb.pos+recomb.len] = d2;
		d[recomb.nPos] = d3;
	}
	else{
		d1 = d[nPos];
		d2 = recomb.d;
		if (recomb.pos+recomb.len < M)
			d3 = min_element(d[recomb.pos], d[recomb.pos+recomb.len]);
		
		memcpy(b+recomb.nPos, a+recomb.pos, sizeof(int)*recomb.len);
		memcpy(d_tmp+recomb.nPos, d+recomb.pos, sizeof(double)*recomb.len);

		memcpy(b+recomb.nPos+recomb.len, a+recomb.nPos, sizeof(int)*(recomb.nPos - recomb.pos) );
		memcpy(d_tmp+recomb.nPos+recomb.len, d+recomb.nPos, sizeof(double)*(recomb.nPos - recomb.pos) );
		
		memcpy(a+recomb.nPos, b+recomb.nPos, sizeof(int)*(recomb.pos+recomb.len-recomb.nPos) );
		memcpy(d+recomb.nPos, d_tmp+recomb.nPos, sizeof(double)*(recomb.pos+recomb.len-recomb.nPos) );
		
		d[recomb.nPos] = d1;
		d[recomb.nPos+recomb.len] = d2;
		if (recomb.pos+recomb.len < M)
			d[recomb.pos+recomb.len] = d3;
	}
}
*/
