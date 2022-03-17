#include "SVHSS.h"

int main() {
    BV_Para *bvPara = (BV_Para *) (malloc(sizeof(*bvPara)));
    BV_SK *bvSk = (BV_SK *) (malloc(sizeof(*bvSk)));
    BV_PK *bvPk = (BV_PK *) (malloc(sizeof(*bvPk)));
    BV_Gen(bvPara, bvSk, bvPk);

    hPara pa;
    int repeatTime = 1, task;
    printf("Choose a computing task： 1. compute a degree-d term  2. compute a  degree-d polynomial with n variables\n");
    printf("Enter the index (1 or 2) of the task:\n");
    scanf("%d", &task);
    if (task == 1) {
        printf("Enter d, t, k, m, e.g., 2 1 1 2\n");
        scanf("%d %d %d %d", &pa.d, &pa.t, &pa.k, &pa.m);
        pa.varNum = 1;
    }else if (task == 2){
        printf("Enter d, n, t, k, m, e.g., 2 10 1 1 2\n");
        scanf("%d %d %d %d %d", &pa.d, &pa.varNum, &pa.t, &pa.k, &pa.m);
    }else {
        printf("Error! The input is not in {1, 2}\n");
    }
    pa.sNumPerVar = getComNum(pa.m, pa.t);

    basicFunction(&pa, bvPara, task);
    SVHSS(repeatTime, pa, bvPara, bvPk, bvSk);

    return 0;
}
