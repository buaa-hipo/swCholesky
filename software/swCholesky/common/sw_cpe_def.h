#ifndef SW_CPE_DEF
#define SW_CPE_DEF

#define SLAVES 64
#define COLS 8
#define ROWS 8
#define COL(x) (x & 0x07)
#define ROW(x) ((x & 0x38) >> 3)


#endif
