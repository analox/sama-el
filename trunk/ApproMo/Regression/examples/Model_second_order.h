#include "Array.h"
#include "nrutil.h"
#include "stdio.h"
#include "ludcmp.h"
#include "lubksb.h"
#include "sqr.h"

int Model_second_order(Array<double> InputData,
		       Array<double> TargetData,
		       int start,
		       int ende,
		       Array<double> &coefficient);




int Model_second_order_interaction(Array<double> InputData,
				   Array<double> TargetData,
				   int start,
				   int ende,
				   Array<double> &coefficient);
