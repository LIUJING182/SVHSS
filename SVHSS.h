#ifndef SVHSS_SVHSS_H
#define SVHSS_SVHSS_H

#include "split.h"
#include "compute.h"

 /*!
  * Generate the basic polynomial according parameters
  * @param pa: the parameters from SVHSS
  * @param bvPara: the parameters from cryptosystem
  * @param type: the type of computation task
  */
void basicFunction(hPara *pa, BV_Para *bvPara, int type);

/*!
 * The running process of SVHSS scheme
 * @param repeatTime: the number of times the algorithm is repeated
 * @param pa: the parameters from SVHSS
 * @param bvPara: the parameters from cryptosystem
 * @param bvPk: the public key from cryptosystem
 * @param bvSk: the secret key from cryptosystem
 */
void SVHSS(int repeatTime, hPara pa, BV_Para* bvPara, BV_PK* bvPk, BV_SK* bvSk);
#endif //SVHSS_SVHSS_H
