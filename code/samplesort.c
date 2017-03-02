#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

/* Total number of elements is determined by MAX.	*/
int MAX=10000000;
/* set 1 to debug. otherwise set it to 0.	*/
int DEBUG_FLAG=0;

/*	PrefixSum is used for calculating bucket offest from bucket size.	*/
void prefixSum(int *arr,int low,int high,int *sum)
{
	int l = (high - low ) + 1;
	int i;
	sum[0] = 0;
	for(i=1; i < l; ++i)
	{
		sum[i] = sum[i-1] + arr[i-1];
	}
}

/*	randomSequence is used for generating a random sequence which is used as
	input.	*/
void randomSequence(unsigned int *seedp,int *arr,int size)
{
	int i=0;
	for(i=0; i<size; ++i)
	{
    	rand_r(seedp);
    	arr[i] = *seedp;
	}
}

/*	randomGenerator is used for genearting random number for partitioning.	*/
void randomGenerator(unsigned int *seedp)
{
    rand_r(seedp);
}

/*	Partition subroutine for quicksort algorithm.	*/
int partition(int arr[], int low, int high,unsigned int *seedp)
{
    randomGenerator(seedp);
    unsigned int randomNumber = *seedp;
    int l = (high - low + 1);
    int pivotIndex = low + (randomNumber%l); //generates a random number as a pivot
    int pivot;
    int i = low - 1;
    int j;
    pivot = arr[pivotIndex];
    int temp;
    temp = arr[pivotIndex];
    arr[pivotIndex] = arr[high];
    arr[high] = temp;
    for (j = low; j < high; j++)
    {
        if (arr[j] < pivot)
        {
            i++;
            temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
 
    }
    temp = arr[i+1];
    arr[i+1] = arr[high];
    arr[high] = temp;
    return i + 1;
}

/*	quicksort algorithm for serial local sorting.	*/
void quicksort(int arr[], int p, int q,unsigned int *seedp)
{
    int j;
    if (p < q)
    {
        j = partition(arr, p, q,seedp);
        quicksort(arr, p, j-1,seedp);
        quicksort(arr, j+1, q,seedp);
    }
}

/*	This is the function where all the magic happens.	*/
int main(int argc, char *argv[])
{
	/*	Some initial variables for calculation.	*/
	/*	numOfprocs = number of all the participating processors(in MPI sense, "nodes")	*/
	int numOfprocs;
	/*	rank of individual node.	*/
	int myrank;
	/*	n = total number of	elements to sort.	*/
	int n;
	/*	elementsPerProc = number of elements that a node have at starting.
		elementsPerProc = n/p 	*/
	int elementsPerProc;


	/*	MPI initialisation	*/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	/*	initialising MPI dependent variables.	*/
	n = MAX;
	unsigned int seed = myrank;
	int i,j;

	/*	Stage 1: Generating input	*/
	/*	Genearating input	*/
	elementsPerProc = n/numOfprocs;
	/*	arr = input array.	*/
	int *arr = (int *) malloc(sizeof(int)*elementsPerProc);
	double t1 = MPI_Wtime();
	/*	Actually filling input array with random integers.	*/
	randomSequence(&seed,arr,elementsPerProc);
	double t2 = MPI_Wtime();
	fprintf(stderr,"Time taken for generation of input = %e from rank = %d\n",t2-t1,myrank);

	/*	Storing input to a file. [optional].
		To debug and check correctness of input data.	*/
	if(DEBUG_FLAG == 1)
	{
		char in1[10] = "in";
		char in2[10];
		sprintf(in2, "%d", myrank);
		char in3[10] = ".txt";
		char *inputFile = (char *)malloc(sizeof(char)*(strlen(in1) + strlen(in2) + strlen(in3) + 1));
		sprintf(inputFile,"%s%s%s",in1,in2,in3);
		FILE *inputFP = fopen(inputFile,"w");
		for(i=0; i<elementsPerProc; ++i)
			fprintf(inputFP, "%d\n", arr[i]);

		fclose(inputFP);
	}

	/*	Stage 2: Generation of Local splitters	*/
	/*	Real computation starts from here.	*/
	t1 = MPI_Wtime();
	/*	Sorting locally in order to generate splitters in next stage 
		and to generate buckets in later stages.	*/
	quicksort(arr, 0, elementsPerProc-1, &seed);	

	int numOfSplitters = numOfprocs - 1;	/*	number of splitters = (p-1)	*/
	/*	Number of elements per bloc = number of elements in a bloc in a processor = roughly (n/p)/p 	*/
	int elementsPerBloc = elementsPerProc/numOfprocs;

	int *localSplitters = (int *) malloc(sizeof(int)*numOfSplitters);
	
	/*	Choosing local splitters	*/
	for(j=0;j<numOfSplitters;j++)
	{
		localSplitters[j] = arr[elementsPerBloc*(j+1)];
	}

	/*	Stage 3: generation of global splitters */
	/*	Gathering global splitters	*/
	int *GlobalSplitters = (int *) malloc(sizeof(int)*numOfSplitters);;
	int *CombinedSplitters;
	if(myrank == 0)
	{
		CombinedSplitters = (int *) malloc(numOfprocs * sizeof(int)*numOfSplitters);
	}
	MPI_Gather( localSplitters , numOfSplitters , MPI_INT, CombinedSplitters , numOfSplitters , MPI_INT, 0 , MPI_COMM_WORLD );

	/*	Selecting global splitters	*/
	if(myrank == 0)
	{
		quicksort(CombinedSplitters, 0, (numOfSplitters*numOfprocs)-1, &seed);

		for(j=0; j<numOfSplitters; j++)
		{
			GlobalSplitters[j] = CombinedSplitters[(numOfprocs-1)*(j+1)];
		}
	}

	/*	Scattering global Splitters to each node.	*/
	MPI_Bcast( GlobalSplitters , numOfSplitters , MPI_INT , 0 , MPI_COMM_WORLD);

	/*	Stage 4: Distribution of buckets to nodes where they belong.	*/
	/*	Choosing buckets and bucketsizes	*/
	int sendBucketSize[numOfprocs];
	int sendBucketOffset[numOfprocs];

	j = 0;
	for(i=0; i<numOfprocs ; ++i)
		sendBucketSize[i] = 0;

	for(i=0; i<elementsPerProc; i++)
	{
		if(j < (numOfprocs - 1))
		{
			if(arr[i] < GlobalSplitters[j])
			{
				sendBucketSize[j]++;
			}
			else
			{
				j++;
				i--;
			}
		}
		else
		{
			sendBucketSize[(numOfprocs - 1)]++;
		}
	}

	/*	Offset for sending buckets	*/
	prefixSum(sendBucketSize,0,numOfprocs-1,sendBucketOffset);

	/*	Storing splitters to a file. [optional].
		To debug and check correctness of splitters.	*/
	if(DEBUG_FLAG == 1)
	{
		for(i=0;i<numOfSplitters; ++i)
			printf("%d s from %d\n",GlobalSplitters[i],myrank);
		for(i=0;i<elementsPerProc; ++i)
			printf("%d e from %d\n",arr[i],myrank);
	}

	/*	Passing bucketsizes	*/
	int recvBucketSize[numOfprocs];
	int recvBucketOffset[numOfprocs];
	/*	Sending bucketsizes so that receiving node can allocate space for them.	*/
	MPI_Alltoall( sendBucketSize , 1 , MPI_INT , recvBucketSize , 1 , MPI_INT , MPI_COMM_WORLD );

	prefixSum(recvBucketSize,0,numOfprocs-1,recvBucketOffset);

	/*	ready for passing buckets	*/
	int recvArrSize = recvBucketSize[numOfprocs-1]+recvBucketOffset[numOfprocs-1];
	int *recvArr = (int *)malloc(sizeof(int)*recvArrSize);

	/*	Actually sending and receiving variable length buckets.	*/
	MPI_Alltoallv( arr , sendBucketSize , sendBucketOffset , MPI_INT , recvArr , recvBucketSize , recvBucketOffset , MPI_INT , MPI_COMM_WORLD );

	/*	To debug and check correctness of receiving buckets.	*/
	if(DEBUG_FLAG == 1)
	{
		for(i=0;i<recvArrSize; ++i)
			printf("%d r from %d\n",recvArr[i],myrank);
	}

	/*	Stage 5: Generating local sorted array.	*/
	quicksort(recvArr, 0, recvArrSize-1, &seed);
	t2 = MPI_Wtime();

	/*	Storing output to a file. [optional].
		To debug and check correctness of output sorted data.	*/
	if(DEBUG_FLAG == 1)
	{
		char out1[10] = "out";
		char out2[10];
		sprintf(out2, "%d", myrank);
		char out3[10] = ".txt";
		char *outputFile = (char *)malloc(sizeof(char)*(strlen(out1) + strlen(out2) + strlen(out3) + 1));
		sprintf(outputFile,"%s%s%s",out1,out2,out3);
		FILE *outputFP = fopen(outputFile,"w");
		for(i=0; i<recvArrSize; ++i)
			fprintf(outputFP, "%d\n", recvArr[i]);

		fclose(outputFP);
	}
	free(arr);
	free(localSplitters);
	free(GlobalSplitters);
	free(recvArr);
	
	MPI_Finalize();
	fprintf(stderr,"Time taken for output = %e from rank = %d\n",t2-t1,myrank);
	return 0;
}