#include <Rcpp.h>
using namespace Rcpp;

double size(NumericVector x) {
	return sqrt(x.size());
}


NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3) {
	int n1 = x1.size();
	int n2 = x2.size();
	int n3 = x3.size();

	int n = n1 + n2 + n3;

	NumericVector x(n);
	int offset = 0;
	for(int i = 0; i < n1; i ++) {
		x[i] = x1[i];
	}

	offset = n1;
	for(int i = 0; i < n2; i ++) {
		x[offset + i] = x2[i];
	}

	offset = n1 + n2;
	for(int i = 0; i < n3; i ++) {
		x[offset + i] = x3[i];
	}

	return x;
}


NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4) {
	int n1 = x1.size();
	int n2 = x2.size();
	int n3 = x3.size();
	int n4 = x4.size();

	int n = n1 + n2 + n3 + n4;

	NumericVector x(n);
	int offset = 0;
	for(int i = 0; i < n1; i ++) {
		x[i] = x1[i];
	}

	offset = n1;
	for(int i = 0; i < n2; i ++) {
		x[offset + i] = x2[i];
	}

	offset = n1 + n2;
	for(int i = 0; i < n3; i ++) {
		x[offset + i] = x3[i];
	}

	offset = n1 + n2 + n3;
	for(int i = 0; i < n4; i ++) {
		x[offset + i] = x4[i];
	}

	return x;
}


NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector x5) {
	int n1 = x1.size();
	int n2 = x2.size();
	int n3 = x3.size();
	int n4 = x4.size();
	int n5 = x5.size();

	int n = n1 + n2 + n3 + n4 + n5;

	NumericVector x(n);
	int offset = 0;
	for(int i = 0; i < n1; i ++) {
		x[i] = x1[i];
	}

	offset = n1;
	for(int i = 0; i < n2; i ++) {
		x[offset + i] = x2[i];
	}

	offset = n1 + n2;
	for(int i = 0; i < n3; i ++) {
		x[offset + i] = x3[i];
	}

	offset = n1 + n2 + n3;
	for(int i = 0; i < n4; i ++) {
		x[offset + i] = x4[i];
	}

	offset = n1 + n2 + n3 + n4;
	for(int i = 0; i < n5; i ++) {
		x[offset + i] = x5[i];
	}

	return x;
}


NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector x5,
	NumericVector x6, NumericVector x7) {
	int n1 = x1.size();
	int n2 = x2.size();
	int n3 = x3.size();
	int n4 = x4.size();
	int n5 = x5.size();
	int n6 = x6.size();
	int n7 = x7.size();

	int n = n1 + n2 + n3 + n4 + n5 + n6 + n7;

	NumericVector x(n);
	int offset = 0;
	for(int i = 0; i < n1; i ++) {
		x[i] = x1[i];
	}

	offset = n1;
	for(int i = 0; i < n2; i ++) {
		x[offset + i] = x2[i];
	}

	offset = n1 + n2;
	for(int i = 0; i < n3; i ++) {
		x[offset + i] = x3[i];
	}

	offset = n1 + n2 + n3;
	for(int i = 0; i < n4; i ++) {
		x[offset + i] = x4[i];
	}

	offset = n1 + n2 + n3 + n4;
	for(int i = 0; i < n5; i ++) {
		x[offset + i] = x5[i];
	}

	offset = n1 + n2 + n3 + n4 + n5;
	for(int i = 0; i < n6; i ++) {
		x[offset + i] = x6[i];
	}

	offset = n1 + n2 + n3 + n4 + n5 + n6;
	for(int i = 0; i < n7; i ++) {
		x[offset + i] = x7[i];
	}

	return x;
}


NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, 
	NumericVector x5, NumericVector x6, NumericVector x7, NumericVector x8, NumericVector x9) {
	int n1 = x1.size();
	int n2 = x2.size();
	int n3 = x3.size();
	int n4 = x4.size();
	int n5 = x5.size();
	int n6 = x6.size();
	int n7 = x7.size();
	int n8 = x8.size();
	int n9 = x9.size();

	int n = n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8 + n9;

	NumericVector x(n);
	int offset = 0;
	for(int i = 0; i < n1; i ++) {
		x[i] = x1[i];
	}

	offset = n1;
	for(int i = 0; i < n2; i ++) {
		x[offset + i] = x2[i];
	}

	offset = n1 + n2;
	for(int i = 0; i < n3; i ++) {
		x[offset + i] = x3[i];
	}

	offset = n1 + n2 + n3;
	for(int i = 0; i < n4; i ++) {
		x[offset + i] = x4[i];
	}

	offset = n1 + n2 + n3 + n4;
	for(int i = 0; i < n5; i ++) {
		x[offset + i] = x5[i];
	}

	offset = n1 + n2 + n3 + n4 + n5;
	for(int i = 0; i < n6; i ++) {
		x[offset + i] = x6[i];
	}

	offset = n1 + n2 + n3 + n4 + n5 + n6;
	for(int i = 0; i < n7; i ++) {
		x[offset + i] = x7[i];
	}

	offset = n1 + n2 + n3 + n4 + n5 + n6 + n7;
	for(int i = 0; i < n8; i ++) {
		x[offset + i] = x8[i];
	}

	offset = n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8;
	for(int i = 0; i < n9; i ++) {
		x[offset + i] = x9[i];
	}

	return x;
}

NumericVector c_vec(NumericVector x1, NumericVector x2) {
	int n1 = x1.size();
	int n2 = x2.size();

	int n = n1 + n2;

	NumericVector x(n);
	int offset = 0;
	for(int i = 0; i < n1; i ++) {
		x[i] = x1[i];
	}

	offset = n1;
	for(int i = 0; i < n2; i ++) {
		x[offset + i] = x2[i];
	}

	return x;
}

void move(NumericVector x, NumericVector y, double h = 0, double v = 0) {
	x = x + h;
	y = y + v;
	return;
}

void hmove(NumericVector x, NumericVector y, double h = 0) {
	x = x + h;
	return;
}

void vmove(NumericVector x, NumericVector y, double v = 0) {
	y = y + v;
	return;
}

void hflip(NumericVector x, NumericVector y) {
	double offset = (size(x)-1.0)/2;
	hmove(x, y, -offset);
	x = 0 - x;
	hmove(x, y, offset);
	return;
}

void vflip(NumericVector x, NumericVector y) {
	double offset = (size(x)-1.0)/2;
	vmove(x, y, -offset);
	y = 0 - y;
	vmove(x, y, offset);
	return;
}

// clockwise
void turn(NumericVector x, NumericVector y, int angle = 90) {
	double offset = (size(x) - 1.0)/2;

	NumericVector z;
	if(angle == 90) {
		move(x, y, -offset, -offset);
		z = 0 + x;
		x = 0 + y;
		y = 0 + z;
		x = 0 - x;
		move(x, y, offset, offset);

	} else if(angle == -90) {
		move(x, y, -offset, -offset);
		z = 0 + x;
		x = 0 + y;
		y = 0 + z;
		y = 0 - y;
		move(x, y, offset, offset);

	} else if(angle == -180 || angle == 180) {
		move(x, y, -offset, -offset);
		x = 0 - x;
		y = 0 - y;
		move(x, y, offset, offset);

	}
	return;
}

void reverse(NumericVector x, NumericVector y) {
	NumericVector x2 = clone(x);
	NumericVector y2 = clone(y);

	x = rev(x2);
	y = rev(y2);
	return;
}


void move(List pos, double h_offset, double v_offset) {
	NumericVector x = pos[0];
	NumericVector y = pos[1];

	x = x + h_offset;
	y = y + v_offset;

	return;
}

// 1: bottomleft-topright
// 2: topleft-bottomright
void diag_flip(NumericVector x, NumericVector y, int type = 1) {
	if(type == 1) {
		hflip(x, y);
		turn(x, y, -90);
	} else {
		hflip(x, y);
		turn(x, y, 90);
	}
	return;
}


List c_list(List pos1, List pos2, List pos3) {
	NumericVector x1 = pos1[0];
	NumericVector x2 = pos2[0];
	NumericVector x3 = pos3[0];

	NumericVector x = c_vec(x1, x2, x3);

	NumericVector y1 = pos1[1];
	NumericVector y2 = pos2[1];
	NumericVector y3 = pos3[1];

	NumericVector y = c_vec(y1, y2, y3);

	List pos = List::create(x, y);
	return pos;
}

List c_list(List pos1, List pos2, List pos3, List pos4) {
	NumericVector x1 = pos1[0];
	NumericVector x2 = pos2[0];
	NumericVector x3 = pos3[0];
	NumericVector x4 = pos4[0];

	NumericVector x = c_vec(x1, x2, x3, x4);

	NumericVector y1 = pos1[1];
	NumericVector y2 = pos2[1];
	NumericVector y3 = pos3[1];
	NumericVector y4 = pos4[1];

	NumericVector y = c_vec(y1, y2, y3, y4);

	List pos = List::create(x, y);
	return pos;
}


List c_list(List pos1, List pos2, List pos3, List pos4, List pos5) {
	NumericVector x1 = pos1[0];
	NumericVector x2 = pos2[0];
	NumericVector x3 = pos3[0];
	NumericVector x4 = pos4[0];
	NumericVector x5 = pos5[0];

	NumericVector x = c_vec(x1, x2, x3, x4, x5);

	NumericVector y1 = pos1[1];
	NumericVector y2 = pos2[1];
	NumericVector y3 = pos3[1];
	NumericVector y4 = pos4[1];
	NumericVector y5 = pos5[1];

	NumericVector y = c_vec(y1, y2, y3, y4, y5);

	List pos = List::create(x, y);
	return pos;
}


List c_list(List pos1, List pos2, List pos3, List pos4, List pos5, List pos6, List pos7) {
	NumericVector x1 = pos1[0];
	NumericVector x2 = pos2[0];
	NumericVector x3 = pos3[0];
	NumericVector x4 = pos4[0];
	NumericVector x5 = pos5[0];
	NumericVector x6 = pos6[0];
	NumericVector x7 = pos7[0];

	NumericVector x = c_vec(x1, x2, x3, x4, x5, x6, x7);

	NumericVector y1 = pos1[1];
	NumericVector y2 = pos2[1];
	NumericVector y3 = pos3[1];
	NumericVector y4 = pos4[1];
	NumericVector y5 = pos5[1];
	NumericVector y6 = pos6[1];
	NumericVector y7 = pos7[1];

	NumericVector y = c_vec(y1, y2, y3, y4, y5, y6, y7);

	List pos = List::create(x, y);
	return pos;
}


bool is_bit_one(int x, int pos) {
	bool l = x & (1 << (pos - 1));
	return l;
}

