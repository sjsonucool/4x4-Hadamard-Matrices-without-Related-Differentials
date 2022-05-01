#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include "boxes-ref.dat"

#define   SIZE   4

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
	int i;
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

//arranging elemwnts to form a 4X4 Hadamard Matrix
unsigned char ArrangeHadamard(unsigned char h[4], unsigned char H[4][4])
{
	int i, j;

	for(i = 0; i < 4; i++){
		for(j = 0; j < 4; j++){
			H[i][j] = h[i^j];
		}
	}
}

int main(int argc, char *argv[])
{
	int ii, i, j, k, l, m, n, i1, i2, i3, i4, matcount = 0, mdscount = 0, mixcount=0, ctr=0, ctr1=0, ctr3=0, ctr4=0, ctr5=0, rcount = 0, index1 = 0, coeffctr = 0;
	char *p;
	char filename1[100] = {"noRelfile"};
	char filename2[100] = {".txt"};
	unsigned char relctr[17], relctrinv[17];
	unsigned char det=1;
	unsigned char M[4][4], Minv[4][4], h[4];
	unsigned char subarr1[144], subarr2[144], submat1[3][3], submat2[2][2], temp[1][1], coeffArray[16], coeffMat[4][4];

	strcat(filename1,argv[1]);
	strcat(filename1,filename2);
	FILE *fp;
	fp = fopen(filename1, "w");
	h[0] = strtol(argv[1], &p, 10);

	for(i2 = 1; i2 < 256; i2++){
		h[1] = i2;
		for(i3 = 1; i3 < 256; i3++){
			h[2] = i3;
			for(i4 = 1; i4 < 256; i4++){
				h[3] = i4;
				ArrangeHadamard(h,M);

				if(gdet(4,M) != 0){ mixcount++;}
				
				for(i = 0; i < 144; i++){
					subarr1[i] = 0;subarr2[i] = 0;
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
				for(i = 0; i < 144;i = i+4){	
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
				
				if(mixcount == 1 && ctr3 == 16 && ctr4 == 36 && ctr5 == 16){
					mdscount++;
					
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
					
					rcount = 0;
					
					//new code
					
					if((relctr[0] == 0 || relctr[1] == 0 || relctr[2] == 0 || relctr[3] == 0) && (relctrinv[0] == 0 || relctrinv[1] == 0 || relctrinv[2] == 0 || relctrinv[3] == 0))
						rcount++;
						
					if((relctr[4] == relctr[11] || relctr[4] == relctr[12] || relctr[4] == relctr[13]) && (relctrinv[4] == relctrinv[11] || relctrinv[4] == relctrinv[12] || relctrinv[4] == relctrinv[13]))
						rcount++;
						
					if((relctr[5] == relctr[11] || relctr[5] == relctr[14] || relctr[5] == relctr[15]) && (relctrinv[5] == relctrinv[11] || relctrinv[5] == relctrinv[14] || relctrinv[5] == relctrinv[15]))
						rcount++;
						
					if((relctr[6] == relctr[12] || relctr[6] == relctr[14] || relctr[6] == relctr[16]) && (relctrinv[6] == relctrinv[12] || relctrinv[6] == relctrinv[14] || relctrinv[6] == relctrinv[16]))
						rcount++;
						
					if((relctr[7] == relctr[13] || relctr[7] == relctr[15] || relctr[7] == relctr[16]) && (relctrinv[7] == relctrinv[13] || relctrinv[7] == relctrinv[15] || relctrinv[7] == relctrinv[16]))
						rcount++;
						
					if((relctr[8] == relctr[12] || relctr[8] == relctr[13] || relctr[8] == relctr[14] || relctr[8] == relctr[15]) && (relctrinv[8] == relctrinv[12] || relctrinv[8] == relctrinv[13] || relctrinv[8] == relctrinv[14] || relctrinv[8] == relctrinv[15]))
						rcount++;
						
					if((relctr[9] == relctr[11] || relctr[9] == relctr[13] || relctr[9] == relctr[14] || relctr[9] == relctr[16]) && (relctrinv[9] == relctrinv[11] || relctrinv[9] == relctrinv[13] || relctrinv[9] == relctrinv[14] || relctrinv[9] == relctrinv[16]))
						rcount++;
						
					if((relctr[10] == relctr[11] || relctr[10] == relctr[12] || relctr[10] == relctr[15] || relctr[10] == relctr[16]) && (relctrinv[10] == relctrinv[11] || relctrinv[10] == relctrinv[12] || relctrinv[10] == relctrinv[15] || relctrinv[10] == relctrinv[16]))
						rcount++;
					
					//end
					
					if(rcount == 0){
						for(i = 0; i < 4; i++)
							fprintf(fp, " %d ", h[i]);
						fprintf(fp, "\n");
					}
				}
				mixcount = 0;
				ctr = 0;
				ctr1 = 0;
				ctr3 = 0;
				ctr4 = 0;
				ctr5 = 0;
			}
		}
	}

	printf("\n %d total mds matrices\n", mdscount);
	fclose(fp);
	return 0;
}

