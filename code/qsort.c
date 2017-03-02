#include <stdio.h>
#include <stdlib.h>
#define MAX 10000000

void randomSequence(unsigned int *seedp,int *arr,int size)
{
    int i=0;
    for(i=0; i<size; ++i)
    {
        rand_r(seedp);
        arr[i] = *seedp;
    }
}

void randomGenerator(unsigned int *seedp)
{
    rand_r(seedp);
}

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
 
void quick_sort(int arr[], int p, int q,unsigned int *seedp)
{
    int j;
    if (p < q)
    {
        j = partition(arr, p, q,seedp);
        quick_sort(arr, p, j-1,seedp);
        quick_sort(arr, j+1, q,seedp);
    }
}

int main()
{
    int i;
    int *arr = (int *) malloc(sizeof(int)*MAX);
    unsigned int seed = 5;
    struct timeval start, end;
    time_t t1 = time(NULL);
    randomSequence(&seed,arr,MAX);
    double t2 = time(NULL);
    fprintf(stderr,"Time taken for generation of input = %e\n",t2-t1);
    t1 = time(NULL);
    gettimeofday(&start, NULL);
    quick_sort(arr, 0, MAX-1,&seed);
    gettimeofday(&end, NULL);
    t2 = time(NULL);
    double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    fprintf(stderr,"Time taken for output = %e .\n",delta);
    return 0;
}
