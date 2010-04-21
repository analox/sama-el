/*
 * MODULE:  qsort.c
 * AUTHOR:  Ivan Griffin (ivan.griffin@ul.ie)
 * DATE:    16 December 1997
 *
 * PURPOSE: Implementation of the Quick Sort
 */

#include "sort.h"

static void quicksort(int left, int right, double *array, double *indexarray, int order);

void sort(double array[], double indexarray[], int size, SortOrder order)
{
    quicksort(0, size-1, array, indexarray, order);
}

static void quicksort(int left, int right, double *array, double *indexarray, int order)
{
    double target, swap_temp, swap_temp_index;
    int comparison, leftscan, rightscan;

    if (right > left)
    {
        target = array[left];
        
        leftscan = left;
        rightscan = right;

        while (rightscan > leftscan)
        {
            for (rightscan = rightscan;
                 rightscan > leftscan;
                 rightscan --)
            {
                comparison = 0;

                if (order == ASCENDING)
                {
                    comparison = (array[rightscan] < target);
                }
                else
                {
                    comparison = (array[rightscan] > target);
                }
 
                if (0 != comparison)
                {
                    swap_temp = array[rightscan];
                    swap_temp_index = indexarray[rightscan];
                    
                    array[rightscan] = array[leftscan];
                    indexarray[rightscan] = indexarray[leftscan];
                    
                    array[leftscan] = swap_temp;
                    indexarray[leftscan] = swap_temp_index;
                    break;
                }
            }

            for (leftscan = leftscan+1;
                 leftscan < rightscan;
                 leftscan ++)
            {
                comparison = 0;

                if (order == ASCENDING)
                    comparison = (array[leftscan] > target);
                else
                    comparison = (array[leftscan] < target);

                if (0 != comparison)
                {
                    swap_temp = array[leftscan];
                    swap_temp_index = indexarray[leftscan];                    
                    
                    array[leftscan] = array[rightscan];
                    indexarray[leftscan] = indexarray[rightscan];
                                        
                    array[rightscan] = swap_temp;
                    indexarray[rightscan] = swap_temp_index;                    
                    break;
                }
            }
        }

        quicksort(left, leftscan - 1, array, indexarray, order);
        quicksort(rightscan + 1, right, array, indexarray, order);
    }
}
