#include <Rcpp.h>
using namespace Rcpp;

#include "std_utils.h"

void fold_hc(NumericVector x1, NumericVector y1,
	NumericVector x2, NumericVector y2,
	NumericVector x3, NumericVector y3,
	NumericVector x4, NumericVector y4) {
	
	double offset = size(x1);

	hflip(x1, y1);
	turn(x1, y1, -90);
	vmove(x2, y2, offset);

	// move to top right
	move(x3, y3, offset, offset);

	hflip(x4, y4);
	turn(x4, y4, 90);
	hmove(x4, y4, offset);

	return;
}

void fold_moore(NumericVector x1, NumericVector y1,
	NumericVector x2, NumericVector y2,
	NumericVector x3, NumericVector y3,
	NumericVector x4, NumericVector y4) {

	double offset = size(x1);

	turn(x1, y1, 90);

	turn(x2, y2, 90);
	vmove(x2, y2, offset);

	turn(x3, y3, -90);
	move(x3, y3, offset, offset);

	turn(x4, y4, -90);
	hmove(x4, y4, offset);
	
	return;
}


void fold_l1(NumericVector x1, NumericVector y1,
	NumericVector x2, NumericVector y2,
	NumericVector x3, NumericVector y3,
	NumericVector x4, NumericVector y4) {

	double offset = size(x1);

	turn(x1, y1, -180);

	vmove(x2, y2, offset);

	move(x3, y3, offset, offset);

	turn(x4, y4, 180);
	hmove(x4, y4, offset);
	
	return;
}



void fold_l2(NumericVector x1, NumericVector y1,
	NumericVector x2, NumericVector y2,
	NumericVector x3, NumericVector y3,
	NumericVector x4, NumericVector y4) {

	double offset = size(x1);

	vflip(x1, y1);

	turn(x2, y2, 90);
	vmove(x2, y2, offset);

	turn(x3, y3,-90);
	move(x3, y3, offset, offset);

	vflip(x4, y4);
	hmove(x4, y4, offset);
	
	return;
}



void fold_l3(NumericVector x1, NumericVector y1,
	NumericVector x2, NumericVector y2,
	NumericVector x3, NumericVector y3,
	NumericVector x4, NumericVector y4) {

	double offset = size(x1);

	hflip(x1, y1);
	turn(x1, y1, -90);

	vmove(x2, y2, offset);

	move(x3, y3, offset, offset);

	turn(x4, y4, -180);
	hmove(x4, y4, offset);
	
	return;
}


void fold_l4(NumericVector x1, NumericVector y1,
	NumericVector x2, NumericVector y2,
	NumericVector x3, NumericVector y3,
	NumericVector x4, NumericVector y4) {

	double offset = size(x1);

	vflip(x1, y1);

	turn(x2, y2, 90);
	vmove(x2, y2, offset);

	turn(x3, y3, -90);
	move(x3, y3, offset, offset);

	turn(x4, y4, -90);
	hmove(x4, y4, offset);
	
	return;
}



// [[Rcpp::export]]
List hilbert_curve_cpp(int level, int type = 1) {

	// left, bottom, bottom, right
	if(level >= 2) {
		List pos = hilbert_curve_cpp(level - 1);
		NumericVector x = pos[0];
		NumericVector y = pos[1];

		NumericVector x1 = clone(x);
		NumericVector y1 = clone(y);
		NumericVector x2 = clone(x);
		NumericVector y2 = clone(y);
		NumericVector x3 = clone(x);
		NumericVector y3 = clone(y);
		NumericVector x4 = clone(x);
		NumericVector y4 = clone(y);
		
		if(type == 1) {
			fold_hc(x1, y1, x2, y2, x3, y3, x4, y4);
		} else if(type == 2) {
			fold_moore(x1, y1, x2, y2, x3, y3, x4, y4);
		} else if(type == 3) {
			fold_l1(x1, y1, x2, y2, x3, y3, x4, y4);
		} else if(type == 4) {
			fold_l2(x1, y1, x2, y2, x3, y3, x4, y4);
		} else if(type == 5) {
			fold_l3(x1, y1, x2, y2, x3, y3, x4, y4);
		} else if(type == 6) {
			fold_l4(x1, y1, x2, y2, x3, y3, x4, y4);
		}

		NumericVector x_combine = c_vec(x1, x2, x3, x4);
		NumericVector y_combine = c_vec(y1, y2, y3, y4);

		List pos2 = List::create(x_combine, y_combine);
		return pos2;

	} else {
		// start phass: facing bottom, orientation from left to right
		NumericVector x = NumericVector::create(0, 0, 1, 1);
		NumericVector y = NumericVector::create(0, 1, 1, 0);
		List pos = List::create(x, y);

		return pos;
	}
}
