#include <stdlib.h>
#include <malloc.h>

unsigned int memoryBytesUsed = 0;
void *gtrack_alloc(unsigned int numElements, unsigned int sizeInBytes)
{
	memoryBytesUsed += numElements*sizeInBytes;
    return calloc(numElements, sizeInBytes);
}
void gtrack_free(void *pFree, unsigned int sizeInBytes)
{
	memoryBytesUsed -= sizeInBytes;
    free(pFree);
}