#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>

#include "boxes-ref.dat"


// Dimension of input square matrix
#define N 4
#define SWAP(a, b) do { typeof(a) temp = a; a = b; b = temp; } while (0)
#define f_size 256 //size of the Galois field
#define fstar_size 255
#define m_size 4

//arranging elemwnts to form a 4X4 Hadamard Matrix
unsigned char ArrangeHadamard(unsigned char h[4], unsigned char H[4][4])
{
	int i, j;

	for(i = 0;i < 4;i++){
		for(j = 0;j < 4; j++){
			H[i][j] = h[i^j];
		}
	}
}

//addition in Galois field
unsigned char gadd(unsigned char a, unsigned char b) {
	return a ^ b;
}

unsigned char gmul(unsigned char a, unsigned char b) {
   /* multiply two elements of GF(2^m)
    * needed for MixColumn and InvMixColumn
    */
	if (a && b) return Alogtable[(Logtable[a] + Logtable[b])%255];
	else return 0;
}

//finding modular inverse of an element
unsigned char modInverse(unsigned char a) 
{ 
	int i,j;
	unsigned char res = 1;
	for(i=0;i<254;i++){ /*Multiplying a by itself 254 times since, a^254 = a^-1*/
		res = gmul(res,a);
	}
	return res;
}

//exponentiation in Galois field
unsigned char gpow(unsigned char a, int pow){
	unsigned char res = 1;
	for(int i = 0; i < pow; i++){
		res = gmul(res,a);
	}
	return res;
}

//division in Galois field (a*b^-1 = a/b)
unsigned char gdiv(unsigned char a, unsigned char b) {

	return gmul(a,modInverse(b));
}

// Function to get determinant of matrix with elements in GF(256)
unsigned char gdet(int n, unsigned char a[n][n]){
	int i, j, k;
	unsigned char temp[n][n], ratio, det=1;
	
	/*using Gauss Elimination
    Technique for transforming matrix to
    upper triangular matrix */

	if(n == 1){
		det = a[0][0];
	}
	else{
		for(i = 0; i < n; i++)
			for(j = 0; j < n; j++)
				temp[i][j] = a[i][j];
		
		for(i=0;i< n;i++){
			if(temp[i][i] != 0){
				for(j=i+1;j< n;j++){
					ratio = gdiv(temp[j][i],temp[i][i]);
					for(k=0;k< n;k++){
						temp[j][k] = gadd(temp[j][k],gmul(ratio,temp[i][k]));
					}
				}
			}
		}
		
		/* Finding determinant by multiplying
		elements in principal diagonal elements */
		for(i=0;i< n;i++)
		{
		    det = gmul(det,temp[i][i]);
		}
    }
    return det;
}

//multiply matrices
void multiplyMatrices(int r1, int c1, int c2, unsigned char first[r1][c1],
                      unsigned char second[c1][c2],
                      unsigned char result[r1][c2]) {

	// Initializing elements of matrix mult to 0.
	unsigned char temp = 0;
   
	for(int i = 0; i < r1; i++)
		for(int j = 0; j < c2; j++)
			result[i][j] = 0;
	// Multiplying first and second matrices and storing it in result
	for(int i = 0; i < r1; i++) {
		for(int j = 0; j < c2; j++) {
			for(int k = 0; k < c1; k++) {
				temp = temp^gmul(first[i][k],second[k][j]);
			}
			result[i][j] = temp;
			temp = 0;
		}
	}
}

void addMatrices(int r, int c, unsigned char first[r][c], unsigned char second[r][c], unsigned char result[r][c]){
	for(int i = 0; i < r; i++)
		for(int j = 0; i < c; j++)
			result[i][j] = gadd(first[i][j],second[i][j]);
}

int weight(unsigned char a[4][1]){
	int size;
	int weight = 0;
	for(size = 0; size < 4; size++){
		if(a[size][0] != 0) weight++;
	}
	
	return weight;
}

int main(int argc, char *argv[]){

	int ii, i, j, k, l, m, n, i1, j1, mixcount=0, ctr=0, ctr1=0, ctr3=0, ctr4=0, ctr5=0, index1 = 0, index2 = 0, related = 0, check = 0, bcount = 0, coeffctr = 0, rcount = 0, success = 0, totalMats = 0, elem;
	unsigned char M[4][4], h[4], Minv[4][4], b[4][1], c[4][1], c1[4][1], subarr1[144], subarr2[144], submat1[3][3], submat2[2][2], temp[1][1], Mres[4][4], coeffArray[16], coeffMat[4][4], Mdet;
	unsigned char d1[4][7140];
	unsigned char d2[4][7140];
	unsigned char d3[4][7140];
	unsigned char d4[4][7140];
	unsigned char d5[4][14280];
	unsigned char d6[4][14280];
	char rel[17][100] = {"a+b+c","a+b+d","a+c+d","b+c+d","a(a+b+c+d)","b(a+b+c+d)","c(a+b+c+d)","d(a+b+c+d)","ab+cd", "ac+bd", "ad+bc", "a^2+b^2", "a^2+c^2", "a^2+d^2", "b^2+c^2", "b^2+d^2", "c^2+d^2"};
	char relinv[17][100] = {"e+f+g","e+f+h","e+g+h","f+g+h","e(e+f+g+h)","f(e+f+g+h)","g(e+f+g+h)","h(e+f+g+h)","ef+gh", "eg+fh", "eh+fg", "e^2+f^2", "e^2+g^2", "e^2+h^2", "f^2+g^2", "f^2+h^2", "g^2+h^2"};
	char filename1[100] = {"noRelfile"};
	char filename2[100] = {".txt"};
	char filename3[100] = {"found"};
	unsigned char relctr[17], relctrinv[17];
	
	strcat(filename1,argv[1]);
	strcat(filename3,argv[1]);
	strcat(filename1,filename2);
	strcat(filename3,filename2);
	
	FILE *fptr, *fptr1;
	fptr = fopen(filename1, "r");
	fptr1 = fopen(filename3, "w");

	while(fscanf(fptr, "%d", &elem) != EOF){
	
		h[0] = elem;
		fscanf(fptr, "%d", &elem);
		h[1] = elem;
		fscanf(fptr, "%d", &elem);
		h[2] = elem;
		fscanf(fptr, "%d", &elem);
		h[3] = elem;
		totalMats++;

		ArrangeHadamard(h,M);
		
		Mdet = gdet(4,M);
		
		for(i=0;i<144;i++){
			subarr1[i]=0;
		}

		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				submat1[i][j]=0;
				
		for(i=0;i<4;i++){
			for(j=0;j<4;j++){
				for(k=0;k<4;k++){
					for(l=0;l<4;l++){
						if(i != k && j != l){
							subarr1[ctr] = M[k][l]; ctr++;
						}
						else continue;
					}
				}
			}
		}
		ctr = 0;
		
		for(i=0;i<144;i=i+9){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					submat1[j][k]=subarr1[3*j+k+i];
				}
			}
			if(gdet(3,submat1) != 0){
				coeffArray[ctr3] = gdet(3,submat1);
				ctr3++;
			}
		}
		ctr3 = 0;
		
		for(j = 0; j < 4; j++){
			for(k = 0; k <4; k++){
				coeffMat[j][k] = coeffArray[coeffctr];
				coeffctr++;
			}
		}
		coeffctr = 0;
		
		for(j = 0; j < 4; j++){
			for(k = 0; k < 4; k++){
				Minv[k][j] = gmul(coeffMat[j][k], modInverse(gdet(4,M)));
			}
		}
		
		for(i = 0; i < 4; i++){
			for(j = 0; j < 1; j++)
				b[i][j] = 0;
		}
		
		//code for generating b,c with hamming weights 1 and 4
		for(j = 0; j < 4; j++){
			for(k = 1; k < 256; k++){
				b[j][0] = k;
				multiplyMatrices(4,4,1,M,b,c);
						
				for(i1 = 0; i1 < 4; i1++){
					d1[i1][index1] = b[i1][0];
					d3[i1][index1] = b[i1][0];
				}
							
				for(i1 = 0; i1 < 4; i1++)
					d2[i1][index1] = c[i1][0];
							
				multiplyMatrices(4,4,1,Minv,b,c);
				for(i1 = 0; i1 < 4; i1++)
					d4[i1][index1] = c[i1][0];
						
				index1++;
			}
			b[j][0] = 0;
		}
		index2 = index1;

		for(i = 0; i < 3; i++){
			for(j = i+1; j < 4; j++){
				for(k = 1; k < 256; k++){
					b[i][0] = k;
							
					for(l = 1; l < 256; l++){
						b[j][0] = l;
						multiplyMatrices(4,4,1,M,b,c);
						multiplyMatrices(4,4,1,Minv,b,c1);
								
						if(weight(c) == 3){
							for(i1 = 0; i1 < 4; i1++)
								d1[i1][index1] = b[i1][0];
							
							for(i1 = 0; i1 < 4; i1++)
								d2[i1][index1] = c[i1][0];
										
							index1++;
						}
								
						if(weight(c1) == 3){
							for(i1 = 0; i1 < 4; i1++)
								d3[i1][index2] = b[i1][0];
						
							for(i1 = 0; i1 < 4; i1++)
								d4[i1][index2] = c1[i1][0];
										
							index2++;
						}
					}
					b[j][0] = 0;
				}
				b[i][0] = 0;
			}
		}
		index1 = 0;
		index2 = 0;
				
		for(i = 0; i < 4; i++){
			for(j = 0; j < 7140; j++){
				d5[i][j] = d1[i][j];
				d5[i][j+7140] = d4[i][j];
				d6[i][j] = d2[i][j];
				d6[i][j+7140] = d3[i][j]; 
			}
		}
		
		/*Store the vectors and compare all possible pairs of vector b,b' and check if they are related differences by the definition in paper. If b,b' are related, then check the same for corresponding b.M=c,b'.M=c'. If both pairs are related, they form related differentials (b,c) and (b',c')*/
		for(i = 0; i < 14279; i++){
			for(j = i+1; j < 14280; j++){
				k = 0;
				while((d5[k][i] == d5[k][j] || d5[k][i] == 0 || d5[k][j] == 0) && (d6[k][i] == d6[k][j] || d6[k][i] == 0 || d6[k][j] == 0)){
					k++;
					if(k == 4) break;
				}
				if(k == 4) check++;
			}
		}
		
		if(check == 0){ 
			success++;
		}

		else{
			for(i = 0; i < 4; i++)
				fprintf(fptr1, " %d ", h[i]);
			fprintf(fptr1, "\n");
			fprintf(fptr1, "Total related differentials =  %d\n", check);
		}
		check = 0;
	}
	
	if(success == totalMats){ 
		printf("Success in %s", filename1);
	}
	
	else printf("Counter examples stored in %s",filename3);
	fclose(fptr);
	fclose(fptr1);	
	return 0;
}
