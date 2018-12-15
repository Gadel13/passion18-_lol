#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
void mult(double** &A, double** &B, double** &C, int nA, int mA, int nB, int mB)
{
	if (mA != nB)
		cout << "Block matrix sizes error" << endl;

	for(int i = 0; i < nA;i++)
		for (int j = 0; j < mB; j++)
			for(int k = 0; k < nB;k++)
				C[i][j] += A[i][k]*B[k][j];
}

int main(int argc, char** argv) // A*B = C    ./main <namefile A> <namefile B> <namefile C> 
{
	int numprocs;
	int myid;
	double ** bufA;
	double ** bufB;
	double** bufC;

	int rc;
	if (rc= MPI_Init(&argc, &argv)) 
 	{ 
  		cout << "ERROR! MPI_Init " << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
  	} 

	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 

  	int n = pow(numprocs, 1.0/3);

  	MPI_Comm cube, level, linei, linej, stlb;
	int reorder = 0; 
	int dimcube[3] = {n,n,n};
	int period[3] = {0,0,0};
	MPI_Cart_create (MPI_COMM_WORLD, 3,dimcube,period,reorder,&cube); // NEW Comm CUBE;


	int remain_dims_level[3] = {1,1,0};
	MPI_Cart_sub(cube, remain_dims_level, &level); // SUB LEVELS;

	int remain_dims_linei[3] = {1,0,0};
	MPI_Cart_sub(cube, remain_dims_linei, &linei); // SUB Linei;

	int remain_dims_linej[3] = {0,1,0};
	MPI_Cart_sub(cube, remain_dims_linej, &linej); // SUB Linej;

	int remain_dims_stlb[3] = {0,0,1};
	MPI_Cart_sub(cube, remain_dims_stlb, &stlb);

	MPI_Comm_rank(cube,&myid); 

	char typeA,typeB;
	int nA,mA,nB,mB;

	MPI_File A;
	MPI_File B;
	MPI_Status status;
	MPI_File_open(cube, argv[1],MPI_MODE_RDONLY,MPI_INFO_NULL, &A);
	MPI_File_open(cube, argv[2],MPI_MODE_RDONLY,MPI_INFO_NULL, &B);

	MPI_File_read(A, (char*)&typeA,1, MPI_CHAR, &status);
	MPI_File_read(B, (char*)&typeB,1, MPI_CHAR, &status);

	MPI_File_read(A, (char*)&nA,1, MPI_INT, &status);
	MPI_File_read(A, (char*)&mA,1, MPI_INT, &status);

	MPI_File_read(B, (char*)&nB,1, MPI_INT, &status);
	MPI_File_read(B, (char*)&mB,1, MPI_INT, &status);

	if (mA != nB || nA != mA || nB != mB)
	{
		cout << "ERROR check matrix sizes" << endl;
		return 1;
	}

	int block = nA/n;
	int osta = nA%n;

	int cordinate[3];
	MPI_Cart_coords(cube,myid,3,cordinate);

	int bn = block,bm = block;
	int Bbn = block, Bbm = block;

	if(cordinate[0] == n-1)
		bn = block+osta;
	if(cordinate[2] == n-1)
		bm = block+osta;

	double** blA;
	blA = new double*[bn];
		for(int i = 0; i < bn; i++)
			blA[i] = new double[bm]; 

	if(cordinate[2] == n-1)
		Bbn = block+osta;
	if(cordinate[1] == n-1)
		Bbm = block+osta;

	double** blB;
	blB = new double*[Bbn];
		for(int i = 0; i < Bbn; i++)
			blB[i] = new double[Bbm]; 

	for(int i = 0; i<bn;i++)
		for(int j = 0; j < bm; j++)
			blA[i][j] == 0;
	for(int i = 0; i<Bbn;i++)
		for(int j = 0; j < Bbm; j++)
			blB[i][j] == 0;

	if (cordinate[2] == 0) // level 0
	{	

		if( cordinate[0] == n-1 && cordinate[1] == n-1 )
		{
			bufA = new double*[block+osta];
			for (int i = 0; i < block+osta; i++)
				bufA[i] = new double[block + osta];
			MPI_File_seek(A,sizeof(char) + 2*sizeof(int) + (nA-block-osta)*mA*sizeof(double), MPI_SEEK_SET);
			for(int i = 0; i < block+osta; i++)
			{
				MPI_File_seek(A,(mA-block-osta)*sizeof(double) , MPI_SEEK_CUR);
				for(int j = 0; j < block+osta; j++)
					MPI_File_read(A, (char*)&bufA[i][j], sizeof(double), MPI_CHAR, &status);
			}

				bufB = new double*[block+osta];
			for (int i = 0; i < block+osta; i++)
				bufB[i] = new double[block + osta];
			MPI_File_seek(B,sizeof(char) + 2*sizeof(int) + (nB-block-osta)*mA*sizeof(double), MPI_SEEK_SET);
			for(int i = 0; i < block+osta; i++)
			{
				MPI_File_seek(B,(mB-block-osta)*sizeof(double) , MPI_SEEK_CUR);
				for(int j = 0; j < block+osta; j++)
					MPI_File_read(B, (char*)&bufB[i][j], sizeof(double), MPI_CHAR, &status);
			}
		} else
			if( cordinate[0] == n-1 )
			{
				bufA = new double*[block+osta];
				for (int i = 0; i < block+osta; i++)
					bufA[i] = new double[block];
				MPI_File_seek(A,sizeof(char) + 2*sizeof(int) + (nA-block-osta)*mA*sizeof(double), MPI_SEEK_SET);
				for(int i = 0; i < block+osta; i++)
				{
					MPI_File_seek(A, cordinate[1]*block*sizeof(double) , MPI_SEEK_CUR);
					for(int j = 0; j < block; j++)
						MPI_File_read(A, (char*)&bufA[i][j], sizeof(double), MPI_CHAR, &status);
					MPI_File_seek(A, (mA - (cordinate[1]+1)*block)*sizeof(double) , MPI_SEEK_CUR);

				}

				bufB = new double*[block+osta];
				for (int i = 0; i < block+osta; i++)
					bufB[i] = new double[block];
				MPI_File_seek(B,sizeof(char) + 2*sizeof(int) + (nB-block-osta)*mA*sizeof(double) , MPI_SEEK_SET);
				for(int i = 0; i < block+osta; i++)
				{
					MPI_File_seek(B, cordinate[1]*block*sizeof(double) , MPI_SEEK_CUR);
					for(int j = 0; j < block; j++)
						MPI_File_read(B, (char*)&bufB[i][j], sizeof(double), MPI_CHAR, &status);
					MPI_File_seek(B, (mB - (cordinate[1]+1)*block)*sizeof(double) , MPI_SEEK_CUR);

				}
			} else
				if( cordinate[1] == n-1 )
				{
					bufA = new double*[block];
					for (int i = 0; i < block; i++)
						bufA[i] = new double[block+osta];
					MPI_File_seek(A,sizeof(char) + 2*sizeof(int) + cordinate[0]*mA*block*sizeof(double), MPI_SEEK_SET);
					for(int i = 0; i < block; i++)
					{
						MPI_File_seek(A, cordinate[1]*block*sizeof(double) , MPI_SEEK_CUR);
						for(int j = 0; j < block+osta; j++)
							MPI_File_read(A, (char*)&bufA[i][j], sizeof(double), MPI_CHAR, &status);
					}


					bufB = new double*[block];
					for (int i = 0; i < block; i++)
						bufB[i] = new double[block+osta];
					MPI_File_seek(B,sizeof(char) + 2*sizeof(int) + cordinate[0]*mB*block*sizeof(double), MPI_SEEK_SET);
					for(int i = 0; i < block; i++)
					{
						MPI_File_seek(B, cordinate[1]*block*sizeof(double), MPI_SEEK_CUR);
						for(int j = 0; j < block+osta; j++)
							MPI_File_read(B, (char*)&bufB[i][j], sizeof(double), MPI_CHAR, &status);
					}
				} else
					{
						bufA = new double*[block];
						for (int i = 0; i < block; i++)
							bufA[i] = new double[block];

						MPI_File_seek(A,sizeof(char) + 2*sizeof(int) + cordinate[0]*mA*block*sizeof(double), MPI_SEEK_SET);
						for(int i = 0; i < block; i++)
						{
							MPI_File_seek(A, cordinate[1]*block*sizeof(double) , MPI_SEEK_CUR);
							for(int j = 0; j < block; j++)
								MPI_File_read(A, (char*)&bufA[i][j], sizeof(double), MPI_CHAR, &status);

							MPI_File_seek(A, (mA - (cordinate[1]+1)*block)*sizeof(double) , MPI_SEEK_CUR);
						}


						bufB = new double*[block];
						for (int i = 0; i < block; i++)
							bufB[i] = new double[block];

						MPI_File_seek(B,sizeof(char) + 2*sizeof(int) + cordinate[0]*mB*block*sizeof(double), MPI_SEEK_SET);
						for(int i = 0; i < block; i++)
						{
							MPI_File_seek(B, cordinate[1]*block*sizeof(double) , MPI_SEEK_CUR);
							for(int j = 0; j < block; j++)
								MPI_File_read(B, (char*)&bufB[i][j], sizeof(double), MPI_CHAR, &status);

							MPI_File_seek(B, (mB - (cordinate[1]+1)*block)*sizeof(double) , MPI_SEEK_CUR);
						}
					}	

		bn = block;
		bm = block;
		Bbn = block;
		Bbm = block;
		if(cordinate[0] == n-1)
			bn = block+osta;
		if(cordinate[1] == n-1)
			bm = block+osta;

		if(cordinate[0] == n-1)
			Bbn = block+osta;
		if(cordinate[1] == n-1)
			Bbm = block+osta;

		if(cordinate[1] == 0) 
		{
			for(int i = 0; i < bn; i++)
				for(int j = 0; j < bm; j++)
					blA[i][j] = bufA[i][j];
		} else
		{
			int adress[] = {cordinate[0],cordinate[1],cordinate[1]};
			int rank;
			MPI_Cart_rank(cube,adress,&rank);

			MPI_Send(&bn,(int)1,MPI_INT, rank,1,cube);
			MPI_Send(&bm,(int)1,MPI_INT, rank,1,cube);
			
			for(int i = 0; i < bn; i++)
				MPI_Send(bufA[i],bm,MPI_DOUBLE, rank,1,cube);
		}

		if(cordinate[0] == 0) 
		{

			for(int i = 0; i < Bbn; i++)
				for(int j = 0; j < Bbm; j++)
					blB[i][j] = bufB[i][j];
		} else
		{
			int adress[] = {cordinate[0],cordinate[1],cordinate[0]};
			int rank;
			MPI_Cart_rank(cube,adress,&rank);

			MPI_Send(&Bbn,(int)1,MPI_INT, rank,2,cube);
			MPI_Send(&Bbm,(int)1,MPI_INT, rank,2,cube);
			
			for(int i = 0; i < Bbn; i++)
				MPI_Send(bufB[i],Bbm,MPI_DOUBLE, rank,2,cube);
		}
	} else
	{
		if( cordinate[1] == cordinate[2] )
		{
			MPI_Recv(&bn,(int)1,MPI_INT, MPI_ANY_SOURCE,1,cube, &status);
			MPI_Recv(&bm,(int)1,MPI_INT, MPI_ANY_SOURCE,1,cube, &status);
			for(int i = 0; i < bn; i++)
				MPI_Recv(blA[i],bm,MPI_DOUBLE, MPI_ANY_SOURCE,1,cube, &status);
		}

		if( cordinate[0] == cordinate[2] ) 
		{
			MPI_Recv(&Bbn,(int)1,MPI_INT, MPI_ANY_SOURCE,2,cube, &status);
			MPI_Recv(&Bbm,(int)1,MPI_INT, MPI_ANY_SOURCE,2,cube, &status);
			for(int i = 0; i < Bbn; i++)
				MPI_Recv(blB[i],Bbm,MPI_DOUBLE, MPI_ANY_SOURCE,2,cube, &status);
		}
	}

	MPI_File_close(&A);
	MPI_File_close(&B);

	int rank;
	int adress[1];


	bn = block;
	bm = block;
	Bbn = block;
	Bbm = block;
	if(cordinate[0] == n-1)
		bn = block+osta;
	if(cordinate[2] == n-1)
		bm = block+osta;

	if(cordinate[2] == n-1)
		Bbn = block+osta;
	if(cordinate[1] == n-1)
		Bbm = block+osta;

 	adress[0] = cordinate[2];

 	MPI_Cart_rank(linej,adress,&rank);
	for(int i = 0; i < bn; i++)
		MPI_Bcast(blA[i],bm,MPI_DOUBLE,rank,linej);

 	MPI_Cart_rank(linei,adress,&rank);
	for(int i = 0; i < Bbn; i++)
		MPI_Bcast(blB[i],Bbm,MPI_DOUBLE,rank,linei);

	double** blC;
	blC = new double*[bn];
	for(int i = 0; i < bn; i++)
		blC[i] = new double[Bbm];

	mult(blA, blB, blC, bn,bm,Bbn,Bbm);

	if (cordinate[2] == 0)
	{
		bufC = new double*[bn];
		for(int i = 0; i < bn; i++)
			bufC[i] = new double[Bbm];

		for(int i = 0; i < bn; i++)
			for(int j = 0; j < Bbm; j++)
				bufC[i][j] = 0;
	}

	adress[0] = 0;
	MPI_Cart_rank(stlb,adress,&rank);

	for(int i = 0; i < bn; i++)
		MPI_Reduce(blC[i], bufC[i], Bbm, MPI_DOUBLE, MPI_SUM, rank, stlb);

	MPI_File C;	
	MPI_File_open(cube, argv[3],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL, &C);

	char c = 'd';
	MPI_File_write(C,(char*)&c, sizeof(char),MPI_CHAR,&status);
	MPI_File_write(C,(char*)&nA, sizeof(int),MPI_CHAR,&status);
	MPI_File_write(C,(char*)&mB, sizeof(int),MPI_CHAR,&status);

	if (cordinate[2] == 0) // level 0 WRITE
	{	
		if( cordinate[0] == n-1 && cordinate[1] == n-1 )
		{
			MPI_File_seek(C,sizeof(char) + 2*sizeof(int) + (nA-block-osta)*mB*sizeof(double), MPI_SEEK_SET);
			for(int i = 0; i < block+osta; i++)
			{
				MPI_File_seek(C,(mA-block-osta)*sizeof(double) , MPI_SEEK_CUR);
				for(int j = 0; j < block+osta; j++)
					MPI_File_write(C, (char*)&bufC[i][j], sizeof(double), MPI_CHAR, &status);
			}

		} else
			if( cordinate[0] == n-1 )
			{
				MPI_File_seek(C,sizeof(char) + 2*sizeof(int) + (nA-block-osta)*mA*sizeof(double), MPI_SEEK_SET);
				for(int i = 0; i < block+osta; i++)
				{
					MPI_File_seek(C, cordinate[1]*block*sizeof(double) , MPI_SEEK_CUR);
					for(int j = 0; j < block; j++)
						MPI_File_write(C, (char*)&bufC[i][j], sizeof(double), MPI_CHAR, &status);
					MPI_File_seek(C, (mA - (cordinate[1]+1)*block)*sizeof(double) , MPI_SEEK_CUR);
				}
			} else
				if( cordinate[1] == n-1 )
				{
					MPI_File_seek(C,sizeof(char) + 2*sizeof(int) + cordinate[0]*mA*block*sizeof(double), MPI_SEEK_SET);
					for(int i = 0; i < block; i++)
					{
						MPI_File_seek(C, cordinate[1]*block*sizeof(double) , MPI_SEEK_CUR);
						for(int j = 0; j < block+osta; j++)
							MPI_File_write(C, (char*)&bufC[i][j], sizeof(double), MPI_CHAR, &status);
					}
				} else
					{
						MPI_File_seek(C,sizeof(char) + 2*sizeof(int) + cordinate[0]*mA*block*sizeof(double), MPI_SEEK_SET);
						for(int i = 0; i < block; i++)
						{
							MPI_File_seek(C, cordinate[1]*block*sizeof(double) , MPI_SEEK_CUR);
							for(int j = 0; j < block; j++)
								MPI_File_write(C, (char*)&bufC[i][j], sizeof(double), MPI_CHAR, &status);

							MPI_File_seek(C, (mA - (cordinate[1]+1)*block)*sizeof(double) , MPI_SEEK_CUR);
						}
					}	
	}



	MPI_File_close(&C);

	MPI_Finalize();
	return 0;
}
