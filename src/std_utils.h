#ifndef __STD_UTILS__
#define __STD_UTILS__

double size(NumericVector x);
NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4);
NumericVector c_vec(NumericVector x1, NumericVector x2);
NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3);
NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector x5, NumericVector x6, NumericVector x7, NumericVector x8, NumericVector x9);
NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector x5);
NumericVector c_vec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector x5, NumericVector x6, NumericVector x7);
void move(NumericVector x, NumericVector y, double h = 0, double v = 0);
void hmove(NumericVector x, NumericVector y, double h = 0);
void vmove(NumericVector x, NumericVector y, double v = 0);
void hflip(NumericVector x, NumericVector y);
void vflip(NumericVector x, NumericVector y);
void turn(NumericVector x, NumericVector y, int angle = 90);
void reverse(NumericVector x, NumericVector y);
void diag_flip(NumericVector x, NumericVector y, int type = 1);
void move(List pos, double h_offset, double v_offset);
List c_list(List pos1, List pos2, List pos3, List pos4);
List c_list(List pos1, List pos2, List pos3, List pos4, List pos5);
List c_list(List pos1, List pos2, List pos3);
List c_list(List pos1, List pos2, List pos3, List pos4, List pos5, List pos6, List pos7);
bool is_bit_one(int x, int pos);

#endif
