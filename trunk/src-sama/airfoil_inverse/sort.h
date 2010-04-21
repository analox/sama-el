/*
 * HEADER:  sort.h
 * AUTHOR:  Ivan Griffin (ivan.griffin@ul.ie)
 * DATE:    2 December 1997
 *
 * PURPOSE: Generic sort interface for sort implementation modules
 *          (bubble.c, select.c and qsort.c)
 *
 */

#ifndef _sort_h
#define _sort_h

enum _SortOrder
{
    ASCENDING, DESCENDING
};
typedef enum _SortOrder SortOrder;

void sort(double array[], double indexarray[], int arraySize, SortOrder sortOrder);

#endif
