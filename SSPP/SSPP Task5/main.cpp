#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void mult(double* A, double* B, double* &C, int nA, int mA, int nB, int mB)
{
	if (mA != nB)
		cout << "Block matrix sizes error" << endl;

	for(int i = 0; i < nA;i++)
		for (int j = 0; j < mB; j++)
			for(int k = 0; k < nB;k++)
				C[i*mB + j] += A[i*mA + k]*B[k*mB + j];
}

int main(int argc, char** argv) // A*B = C    ./main <namefile A> <namefile B> <namefile C> 
{
	int numprocs;
	int myid;
	double* bufA;
	double* bufB;
	double* bufC;

	double maxtime = 0.0;
  	double maxtimeIO = 0.0;

  	double time = 0.0;
  	double timeIO = 0.0;


	int rc;
	if (rc= MPI_Init(&argc, &argv)) 
 	{ 
  		cout << "ERROR! MPI_Init " << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
  	} 

	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 

	int n = 1;
	while (n * n * n != numprocs)
		++n;

	MPI_Datatype filetype;


  	MPI_Comm cube, linei, linej, stlb, level;
	int reorder = 0; 
	int dimcube[3] = {n,n,n};
	int period[3] = {0,0,0};
	MPI_Cart_create (MPI_COMM_WORLD, 3,dimcube,period,reorder,&cube); // NEW Comm CUBE;

	MPI_Comm_rank(cube,&myid);

	int remain_dims_level[3];
	remain_dims_level[0] = 1;
	remain_dims_level[1] = 1;
	remain_dims_level[2] = 0;
	MPI_Cart_sub(cube, remain_dims_level, &level); // SUB LEVELS;

	char typeA,typeB;
	int nA,mA,nB,mB;

	MPI_File A;
	MPI_File B;
	MPI_File C;	
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

	int cordinate[3];
	MPI_Cart_coords(cube,myid,3,cordinate);


	int block = nA/n;
	int osta = nA%n;
	int bn = block, bm = block, Bbn = block, Bbm = block;

	if(cordinate[0] == n-1)
		bn = block+osta;
	if(cordinate[2] == n-1)
		bm = block+osta;

	if(cordinate[2] == n-1)
		Bbn = block+osta;
	if(cordinate[1] == n-1)
		Bbm = block+osta;

	double* blA;
	blA = new double[bn*bm];

	double* blB;
	blB = new double[Bbn*Bbm];

	for(int i = 0; i<bn*bm;i++)
		blA[i] = 0;
	for(int i = 0; i<Bbn*Bbm;i++)
		blB[i] = 0;
	

	MPI_File_close(&A);
	MPI_File_close(&B);

	if (cordinate[2] == 0) // level 0
	{	
		timeIO -= MPI_Wtime();
		MPI_File_open(level, argv[1],MPI_MODE_RDONLY,MPI_INFO_NULL, &A);
		MPI_File_open(level, argv[2],MPI_MODE_RDONLY,MPI_INFO_NULL, &B);
		MPI_File_open(level, argv[3],MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &C);
		int start[2];
		if( cordinate[0] == n-1 && cordinate[1] == n-1 )
		{
			char c = 'd';
			MPI_File_write(C,(char*)&c, sizeof(char),MPI_CHAR,&status);
			MPI_File_write(C,(char*)&nA, sizeof(size_t),MPI_CHAR,&status);
			MPI_File_write(C,(char*)&mB, sizeof(size_t),MPI_CHAR,&status);
			bn = block+osta;
			bm = block+osta;
			Bbn = block+osta;
			Bbm = block+osta;
			start[0] = nA - block - osta;
			start[1] = mA - block - osta;

		} else
			if( cordinate[0] == n-1 )
			{
				bn = block+osta;
				bm = block;
				Bbn = block+osta;
				Bbm = block;
				start[0] = nA - block - osta;
				start[1] = cordinate[1]*block;

			} else
				if( cordinate[1] == n-1 )
				{
					bn = block;
					bm = block+osta;
					Bbn = block;
					Bbm = block+osta;
					start[0] = cordinate[0]*block;
					start[1] = mA - block - osta;
				} else
					{
						bn = block;
						bm = block;
						Bbn = block;
						Bbm = block;
						start[0] = cordinate[0]*block;
						start[1] = cordinate[1]*block;
					}	

		bufA = new double[bn*bm];
		bufB = new double[Bbn*Bbm];

		int gsize[2];
		gsize[0] = nA,
		gsize[1] = mA;
		int lsize[2];
		lsize[0] = bn;
		lsize[1] = bm;

		MPI_Type_create_subarray(2, gsize, lsize, start, MPI_ORDER_C, MPI_DOUBLE, &filetype);
        MPI_Type_commit(&filetype);


		MPI_File_set_view(A, sizeof(char)+2*sizeof(int),MPI_DOUBLE,filetype, "native", MPI_INFO_NULL);
        MPI_File_set_view(B, sizeof(char)+2*sizeof(int),MPI_DOUBLE,filetype, "native", MPI_INFO_NULL);



        MPI_File_read(A, bufA, bn*bm, MPI_DOUBLE, &status);
        MPI_File_read(B, bufB, Bbn*Bbm, MPI_DOUBLE, &status);


        MPI_File_close(&A);
		MPI_File_close(&B);
		
		timeIO += MPI_Wtime();

		bn = block;
		bm = block;
		Bbn = block;
		Bbm = block;

		if(cordinate[0] == n-1)
		{
			bn = block + osta;
			Bbn = block + osta;
		}
		if(cordinate[1] == n-1)
		{
			bm = block + osta;
			Bbm = block + osta;
		}

		if(cordinate[1] == 0) 
		{
			for(int i = 0; i < bn*bm; i++)
				blA[i] = bufA[i];
		} else
		{
			int rank;
			int adress[3];
			adress[0] = cordinate[0];
			adress[1] = cordinate[1];
		 	adress[2] = cordinate[1];

			MPI_Cart_rank(cube,adress,&rank);

			MPI_Send(&bn,(int)1,MPI_INT, rank,1,cube);
			MPI_Send(&bm,(int)1,MPI_INT, rank,1,cube);
			
			MPI_Send(bufA,bn*bm,MPI_DOUBLE, rank,1,cube);
		}

		if(cordinate[0] == 0) 
		{

			for(int i = 0; i < Bbn*Bbm; i++)
				blB[i] = bufB[i];
		} else
		{
			int rank;
			int adress[3];
			adress[0] = cordinate[0];
			adress[1] = cordinate[1];
		 	adress[2] = cordinate[0];

			MPI_Cart_rank(cube,adress,&rank);

			MPI_Send(&Bbn,(int)1,MPI_INT, rank,2,cube);
			MPI_Send(&Bbm,(int)1,MPI_INT, rank,2,cube);
			
			MPI_Send(bufB,Bbn*Bbm,MPI_DOUBLE, rank,2,cube);
		}


	} else
	{
		if( cordinate[1] == cordinate[2] )
		{
			MPI_Recv(&bn,(int)1,MPI_INT, MPI_ANY_SOURCE,1,cube, &status);
			MPI_Recv(&bm,(int)1,MPI_INT, MPI_ANY_SOURCE,1,cube, &status);
			MPI_Recv(blA,bn*bm,MPI_DOUBLE, MPI_ANY_SOURCE,1,cube, &status);
		}

		if( cordinate[0] == cordinate[2] ) 
		{
			MPI_Recv(&Bbn,(int)1,MPI_INT, MPI_ANY_SOURCE,2,cube, &status);
			MPI_Recv(&Bbm,(int)1,MPI_INT, MPI_ANY_SOURCE,2,cube, &status);
			MPI_Recv(blB,Bbn*Bbm,MPI_DOUBLE, MPI_ANY_SOURCE,2,cube, &status);
		}
	}

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

	MPI_Barrier(level);

	int rank;
	int adress[3];

	adress[0] = cordinate[0];
	adress[1] = cordinate[2];
 	adress[2] = cordinate[2];

 	int remain_dims_linei[3];
 	remain_dims_linei[0] = false;
 	remain_dims_linei[1] = true;
 	remain_dims_linei[2] = false;
	MPI_Cart_sub(cube, remain_dims_linei, &linei); // SUB Linei;

 	
 	MPI_Cart_rank(linei,&cordinate[2],&rank);
		MPI_Bcast(blA,bn*bm,MPI_DOUBLE,rank,linei);

	MPI_Comm_free(&linei);

	adress[0] = cordinate[2];
	adress[1] = cordinate[1];
 	adress[2] = cordinate[2];


	int remain_dims_linej[3];
	remain_dims_linej[0] = true;
	remain_dims_linej[1] = false;
	remain_dims_linej[2] = false;

	MPI_Cart_sub(cube, remain_dims_linej, &linej); // SUB Linej;

 	MPI_Cart_rank(linej,&cordinate[2],&rank);
		MPI_Bcast(blB,Bbn*Bbm,MPI_DOUBLE,rank,linej);

	MPI_Comm_free(&linej);




	double* blC;
	blC = new double[bn*Bbm];

	for(int i = 0; i < bn*Bbm;i++)
		blC[i] = 0;

	MPI_Barrier(cube);

	time -= MPI_Wtime();
	mult(blA, blB, blC, bn,bm,Bbn,Bbm);
	time += MPI_Wtime();




	MPI_Barrier(cube);

	MPI_Reduce(&time, &maxtime, (int)1, MPI_DOUBLE, MPI_MAX, 0, cube);

	if (cordinate[2] == 0)
	{
		delete[] bufA;
		delete[] bufB;
		bufC = new double[bn*Bbm];

		for(int i = 0; i < bn*Bbm; i++)
			bufC[i] = 0;
	}

	int remain_dims_stlb[3] = {false,false,true};
	MPI_Cart_sub(cube, remain_dims_stlb, &stlb);

	adress[0] = cordinate[0];
	adress[1] = cordinate[1];
	adress[2] = 0;
	MPI_Cart_rank(stlb,&adress[2],&rank);

	MPI_Barrier(cube);
	
	MPI_Reduce(blC, bufC, bn*Bbm, MPI_DOUBLE, MPI_SUM, rank, stlb);


	MPI_Comm_free(&stlb);

	MPI_Barrier(cube);

	if (cordinate[2] == 0) // level 0
	{	

		timeIO -= MPI_Wtime();

		if( cordinate[0] == n-1 && cordinate[1] == n-1 )
		{
			bn = block+osta;
			bm = block+osta;
			Bbn = block+osta;
			Bbm = block+osta;
		} else
			if( cordinate[0] == n-1 )
			{
				bn = block+osta;
				bm = block;
				Bbn = block+osta;
				Bbm = block;
			} else
				if( cordinate[1] == n-1 )
				{
					bn = block;
					bm = block+osta;
					Bbn = block;
					Bbm = block+osta;
				} else
					{
						bn = block;
						bm = block;
						Bbn = block;
						Bbm = block;
					}	

		MPI_File_set_view(C, sizeof(char)+2*sizeof(int),MPI_DOUBLE,filetype, "native", MPI_INFO_NULL);

        MPI_File_write(C, bufC, bn*Bbm, MPI_DOUBLE, &status);
		
		MPI_File_close(&C);
		timeIO += MPI_Wtime();
	}


	MPI_Comm_free(&level);
	MPI_Comm_free(&cube);

	if (cordinate[2] == 0)
		delete[] bufC;

	delete[] blA;
	delete[] blB;
	delete[] blC;


	MPI_Reduce(&timeIO, &maxtimeIO, (int)1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	
	if (myid == 0)
	{
		ofstream MTIO("MTIO.txt", ios_base::out|ios::app);
		ofstream MT("MT.txt", ios_base::out|ios::app);
		MTIO << nA << " " << numprocs << " " << maxtimeIO << endl;
		MT << nA << " " << numprocs << " " << maxtime << endl;
		MTIO.close();
		MT.close();
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}



