//Shows example of MDS Hadamard matrices that allows Related Differentials and the relation(s) satisfied by the matrix elements 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <unistd.h>
#include <string.h>
#include "boxes-ref.dat"

//arranging elements to form a 4X4 Hadamard Matrix
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
	int i,j;
	unsigned char res = 1;
	
	for(i = 0; i < pow; i++){
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
		
		for(i = 0;i < n; i++){
			if(temp[i][i] != 0){
				for(j = i+1; j < n; j++){
					ratio = gdiv(temp[j][i],temp[i][i]);
					for(k = 0; k < n; k++){
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

int main(){

	int ii, i, j, k, l, m, n, i1, j1, mixcount=0, ctr=0, ctr1=0, ctr3=0, ctr4=0, ctr5=0, index1 = 0, index2 = 0, relatedb = 0, relatedc = 0, check = 0, coeffctr = 0;
	unsigned char h[4], M[4][4], Minv[4][4], b[4][1], c[4][1], c1[4][1], subarr1[144], subarr2[144], submat1[3][3], submat2[2][2], temp[1][1], Mres[4][4], GFMT[256][256], coeffArray[16], coeffMat[4][4], Mdet = 0;
	 
	unsigned char d1[4][7140];
	unsigned char d2[4][7140];
	unsigned char d3[4][7140];
	unsigned char d4[4][7140];
	unsigned char d5[4][14820];
	unsigned char d6[4][14820];
	char rel[17][100] = {"a+b+c","a+b+d","a+c+d","b+c+d","a(a+b+c+d)","b(a+b+c+d)","c(a+b+c+d)","d(a+b+c+d)","ab+cd", "ac+bd", "ad+bc", "a^2+b^2", "a^2+c^2", "a^2+d^2", "b^2+c^2", "b^2+d^2", "c^2+d^2"};
	char relinv[17][100] = {"e+f+g","e+f+h","e+g+h","f+g+h","e(e+f+g+h)","f(e+f+g+h)","g(e+f+g+h)","h(e+f+g+h)","ef+gh", "eg+fh", "eh+fg", "e^2+f^2", "e^2+g^2", "e^2+h^2", "f^2+g^2", "f^2+h^2", "g^2+h^2"};
	unsigned char relctr[17], relctrinv[17];
	
	
	srand(time(NULL));

	while(1){
	
		for(i = 0; i < 4; i++)
			h[i] = rand()%256;

		ArrangeHadamard(h,M);

		//determinant of original matrix
		Mdet = gdet(4,M);
		if(gdet(4,M) != 0){ mixcount++;}

		for(i = 0;i < 144; i++){
			subarr1[i] = 0; subarr2[i] = 0;
		}

		for(i = 0; i < 3; i++)
			for(j = 0; j < 3; j++)
				submat1[i][j] = 0;
		
		for(i = 0; i < 2; i++)
			for(j = 0; j < 2; j++)
				submat2[i][j] = 0;

		for(i = 0; i < 4; i++){
			for(j = 0; j < 4; j++){
				for(k = 0; k < 4; k++){
					for(l = 0; l < 4; l++){
						if(i != k && j != l){
							subarr1[ctr] = M[k][l]; ctr++;
						}
						else continue;
					}
				}
			}
		}

		//check determinants of submatrices of size 3
		for(i = 0; i < 144; i = i+9){
			for(j = 0; j < 3; j++){
				for(k = 0; k < 3; k++){
					submat1[j][k] = subarr1[3*j+k+i];
				}
			}
			if(gdet(3,submat1) != 0){
				coeffArray[ctr3] = gdet(3,submat1);
				ctr3++;
			}
		}

		for(j = 0; j < 4; j++){
			for(k = 0; k <4; k++){
				coeffMat[j][k] = coeffArray[coeffctr];
				coeffctr++;
			}
		}
		coeffctr = 0;
		//generate M^-1
		for(j = 0; j < 4; j++){
			for(k = 0; k < 4; k++){
				Minv[k][j] = gmul(coeffMat[j][k], modInverse(gdet(4,M)));
			}
		}

		for(i = 0; i < 3; i++){
			for(j = i+1; j < 4; j++){
				if(i != j){
					for(k = 0; k < 3; k++){
						for(l = k+1; l < 4; l++){
							if(k != l){
								for(m = 0; m < 4; m++){
									for(n = 0; n < 4; n++){
										if(m != i  && m != j && n != k && n != l){
										subarr2[ctr1] = M[m][n];
										ctr1++;
										}
										else continue;
									}
								}
							}
						}
					}
				}
			}
		}

		//check determinants of submatrices of size 2
		for(i = 0;i < 144;i = i+4){
			for(j = 0; j < 2; j++){
				for(k = 0; k < 2; k++){
					submat2[j][k] = subarr2[2*j+k+i];
				}
			}
			if(gdet(2,submat2) != 0) ctr4++;
		}
		
		//check determinants of submatrices of size 1
		
		for(i = 0; i < 4; i++){
			for(j = 0; j < 4; j++){
				temp[0][0] = M[i][j];
				if(gdet(1,temp) != 0) ctr5++;
			}
		}
		
		/*The counters mixcount, ctr3, ctr4 and ctr5 are used to keep track if the determinants of all the submatrices of the random Haramard Matrix M are nonzero. There can be total 16 submatrices of size 3, and total 4C2 X 4C2 = 36 submatrices of size 2. */
		
		if(mixcount == 1 && ctr3 == 16 && ctr4 == 36 && ctr5 == 16){
			for(i = 0; i < 4; i++){
				for(j = 0; j < 1; j++)
					b[i][j] = 0;
			}
			//code for generating b,c with hamming weights (1,4)/(4,1)
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
			//code for generating b,c with hamming weights (2,3)/(3,2)
			for(i = 0; i < 3; i++){
				for(j = i+1; j < 4; j++){
					for(k = 1; k < 256; k++){
						b[i][0] = k;
						
						for(l = 1; l < 256; l++){
							b[j][0] = l;
							multiplyMatrices(4,4,1,M,b,c);
							multiplyMatrices(4,4,1,Minv,b,c1);
							
							if(weight(c) <= 3){
								for(i1 = 0; i1 < 4; i1++)
									d1[i1][index1] = b[i1][0];
						
								for(i1 = 0; i1 < 4; i1++)
									d2[i1][index1] = c[i1][0];
									
								index1++;
							}
							
							if(weight(c1) <= 3){
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
			
			/*Store the vectors and compare all possible pairs of vector b,b' and check if they are related differences by the definition in paper. If b,b' are related, then check the same for corresponding M.b=c,M.b'=c'. If both pairs are related, they form related differentials (b,c) and (b',c')*/
			
			for(i = 0; i < 14279; i++){
				for(j = i+1; j < 14280; j++){
					k = 0;
					while((d5[k][i] == d5[k][j] || d5[k][i] == 0 || d5[k][j] == 0) && (d6[k][i] == d6[k][j] || d6[k][i] == 0 || d6[k][j] == 0)){
						k++;
						if(k == 4) break;
					}
					if(k == 4){
						check++;
						printf("b =");
						for(k = 0; k < 4; k++){
							printf(" %x",d5[k][i]);
						}
						printf("\n");
						printf("b' =");
						for(k = 0; k < 4; k++){
							printf(" %x ", d5[k][j]);
						}
						printf("\n");
						printf("c =");
						for(k = 0; k < 4; k++){
							printf(" %x ", d6[k][i]);
						}
						printf("\n");
						printf("c' =");
						for(k = 0; k < 4; k++){
							printf(" %x ", d6[k][j]);
						}
						printf("\n\n");
					}
					relatedb = 0;
					relatedc = 0;
				}
			}
			
			if(check > 0){
				printf("\n\n Related differentials found\n\n");
				printf("\n\nFor the Hadamard MDS matrix\n\n");
				for(i = 0; i < 4; i++){
					for(j = 0; j < 4; j++){
						printf(" %x ", M[i][j]);
					}
					printf("\n");
				}
				
				printf("\n\nAnd its inverse matrix\n\n");
				for(i = 0; i < 4; i++){
					for(j = 0; j < 4; j++){
						printf(" %x ", Minv[i][j]);
					}
					printf("\n");
				}
			
				multiplyMatrices(4,4,4,M,M,Mres);	
				printf("\n M X M is \n");
				for(i = 0; i < 4; i++){
					for(j = 0; j < 4; j++){
						printf(" %x ", Mres[i][j]);
					}
					printf("\n");
				}
				
				relctr[0] = M[0][0]^M[0][1]^M[0][2];
				relctr[1] = M[0][0]^M[0][1]^M[0][3];
				relctr[2] = M[0][0]^M[0][2]^M[0][3];
				relctr[3] = M[0][1]^M[0][2]^M[0][3];
				relctr[4] = gmul(M[0][0],M[0][0]^M[0][1]^M[0][2]^M[0][3]);
				relctr[5] = gmul(M[0][1],M[0][0]^M[0][1]^M[0][2]^M[0][3]);
				relctr[6] = gmul(M[0][2],M[0][0]^M[0][1]^M[0][2]^M[0][3]);
				relctr[7] = gmul(M[0][3],M[0][0]^M[0][1]^M[0][2]^M[0][3]);
				relctr[8] = gadd(gmul(M[0][0],M[0][1]),gmul(M[0][2],M[0][3]));
				relctr[9] = gadd(gmul(M[0][0],M[0][2]),gmul(M[0][1],M[0][3]));
				relctr[10] = gadd(gmul(M[0][0],M[0][3]),gmul(M[0][1],M[0][2]));
				relctr[11] = gadd(gpow(M[0][0],2),gpow(M[0][1],2));
				relctr[12] = gadd(gpow(M[0][0],2),gpow(M[0][2],2));
				relctr[13] = gadd(gpow(M[0][0],2),gpow(M[0][3],2));
				relctr[14] = gadd(gpow(M[0][1],2),gpow(M[0][2],2));
				relctr[15] = gadd(gpow(M[0][1],2),gpow(M[0][3],2));
				relctr[16] = gadd(gpow(M[0][2],2),gpow(M[0][3],2));
	
				relctrinv[0] = Minv[0][0]^Minv[0][1]^Minv[0][2];
				relctrinv[1] = Minv[0][0]^Minv[0][1]^Minv[0][3];
				relctrinv[2] = Minv[0][0]^Minv[0][2]^Minv[0][3];
				relctrinv[3] = Minv[0][1]^Minv[0][2]^Minv[0][3];
				relctrinv[4] = gmul(Minv[0][0],Minv[0][0]^Minv[0][1]^Minv[0][2]^Minv[0][3]);
				relctrinv[5] = gmul(Minv[0][1],Minv[0][0]^Minv[0][1]^Minv[0][2]^Minv[0][3]);
				relctrinv[6] = gmul(Minv[0][2],Minv[0][0]^Minv[0][1]^Minv[0][2]^Minv[0][3]);
				relctrinv[7] = gmul(Minv[0][3],Minv[0][0]^Minv[0][1]^Minv[0][2]^Minv[0][3]);
				relctrinv[8] = gadd(gmul(Minv[0][0],Minv[0][1]),gmul(Minv[0][2],Minv[0][3]));
				relctrinv[9] = gadd(gmul(Minv[0][0],Minv[0][2]),gmul(Minv[0][1],Minv[0][3]));
				relctrinv[10] = gadd(gmul(Minv[0][0],Minv[0][3]),gmul(Minv[0][1],Minv[0][2]));
				relctrinv[11] = gadd(gpow(Minv[0][0],2),gpow(Minv[0][1],2));
				relctrinv[12] = gadd(gpow(Minv[0][0],2),gpow(Minv[0][2],2));
				relctrinv[13] = gadd(gpow(Minv[0][0],2),gpow(Minv[0][3],2));
				relctrinv[14] = gadd(gpow(Minv[0][1],2),gpow(Minv[0][2],2));
				relctrinv[15] = gadd(gpow(Minv[0][1],2),gpow(Minv[0][3],2));
				relctrinv[16] = gadd(gpow(Minv[0][2],2),gpow(Minv[0][3],2));
				
				
				printf("\n\nRelations in M\n\n");
				for(i = 0; i < 4; i++){
					if(relctr[i] == 0){
						printf("\n%s = %x\n", rel[i],relctr[i]);
					}
				}
				
				for(i = 4; i < 8; i++){
					for(j = 8; j < 17; j++){
						if(relctr[i] == relctr[j]){
							printf("\n%s = %s = %x\n", rel[i],rel[j],relctr[i]);
						}
					}
				}
				
				for(i = 8; i < 11; i++){
					for(j = 11; j < 17; j++){
						if(relctr[i] == relctr[j]){
							printf("\n%s = %s = %x\n", rel[i],rel[j],relctr[i]);
						}
					}
				}
				
				printf("\n\nRelations in Minv\n\n");
				for(i = 0; i < 4; i++){
					if(relctrinv[i] == 0){
						printf("\n%s = %x\n", relinv[i],relctrinv[i]);
					}
				}
				
				for(i = 4; i < 8; i++){
					for(j = 8; j < 17; j++){
						if(relctrinv[i] == relctrinv[j]){
							printf("\n%s = %s = %x\n", relinv[i],relinv[j],relctrinv[i]);
						}
					}
				}
				
				for(i = 8; i < 11; i++){
					for(j = 11; j < 17; j++){
						if(relctrinv[i] == relctrinv[j]){
							printf("\n%s = %s = %x\n", relinv[i],relinv[j],relctrinv[i]);
						}
					}
				}
				
				printf("\nTotal related differentials =  %d ", check);
				printf("\n\n======================================\n\n");
				system("read -p 'Press Enter to continue...' var");
			}
			check = 0;
		}
		
		mixcount = 0;
		ctr = 0;
		ctr1 = 0;
		ctr3 = 0;
		ctr4 = 0;
		ctr5 = 0;
	}
	return 0;
}
