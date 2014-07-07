#ifndef _ALGEBRAOP_H_
#define _ALGEBRAOP_H_

struct omxMatrix;
typedef void (*algebra_op_t)(FitContext *fc, int want, struct omxMatrix**, int, struct omxMatrix*);

#endif
