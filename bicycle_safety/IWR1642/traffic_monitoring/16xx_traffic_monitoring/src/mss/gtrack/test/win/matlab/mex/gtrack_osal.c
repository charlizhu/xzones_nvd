#include <stdlib.h>
#include <malloc.h>

char myArray[1024*1024*16];
unsigned int memIndex;

void *gtrack_alloc(unsigned int numElements, unsigned int sizeInBytes)
{
	void *ptr = &myArray[memIndex];
	memIndex += (numElements*sizeInBytes + 0x100) & 0xFFFFFFFF00;
    return ptr;
}
void gtrack_free(void *ptr, unsigned int sizeInBytes)
{
}
