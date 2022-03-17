#include "SVHSS.h"

func aFunc;
fmpz_t *vars;
long long allSTermNum=0;

void basicFunction(hPara *pa, BV_Para *bvPara, int type){
    vars = malloc(sizeof(fmpz_t) * pa->varNum);
    flint_rand_t state;
    flint_randinit(state);
    for(int i=0; i<pa->varNum;i++)
    {
        fmpz_init(vars[i]);
        fmpz_randm(vars[i], state,bvPara->msg);
    }
    /*********************prepare the computation task*********************/
    if (type == 1){
        aFunc.tNum=1;
        aFunc.tDegree=(int *)malloc(sizeof(int)*aFunc.tNum);
        aFunc.tSymbol=(int *)malloc(sizeof(int)*aFunc.tNum);
        aFunc.tCon=(int **)malloc(sizeof(int*)*aFunc.tNum);
        aFunc.tDegree[0]=pa->d;
        aFunc.tSymbol[0]=1;
        aFunc.tCon[0]=(int *)malloc(sizeof(int)*aFunc.tDegree[0]);
        for(int j=0; j<aFunc.tDegree[0]; j++)
        {
            aFunc.tCon[0][j]=0;
        }
        allSTermNum+=pow(pa->sNumPerVar, aFunc.tDegree[0]);

    }else if (type == 2){
        int *tNumByD = malloc(sizeof(int)*(pa->d+1));
        int *tNumByDSum = malloc(sizeof(int)*(pa->d+1));
        tNumByDSum[0]=0;
        for(int i=1; i<=pa->d; i++){
            tNumByD[i]=getComNum(pa->varNum+i-1, i);
            tNumByDSum[i] = tNumByDSum[i-1] + tNumByD[i];
        }

        aFunc.tNum= getComNum(pa->d+pa->varNum, pa->d)-1;//choose maxDg from maxDg+pa.varNum
        aFunc.tDegree=(int *)malloc(sizeof(int)*aFunc.tNum);
        aFunc.tSymbol=(int *)malloc(sizeof(int)*aFunc.tNum);
        aFunc.tCon=(int **)malloc(sizeof(int*)*aFunc.tNum);

        int de=1;
        for(int i=0; i<aFunc.tNum; i++)
        {
            if (i+1>tNumByDSum[de]){de++;}
            aFunc.tDegree[i]=de;
            aFunc.tSymbol[i]=1;
            if (DEBUG){printf("d=%d\n", aFunc.tDegree[i]);}
            aFunc.tCon[i]=(int *)malloc(sizeof(int)*aFunc.tDegree[i]);
            allSTermNum+=pow(pa->sNumPerVar, aFunc.tDegree[i]);
        }
        //find all terms by degree and variable
        allTermsByDV(aFunc, tNumByD, pa->varNum, pa->d);

        if (DEBUG){
            //check aFunc.tCon
            printf("aFunc.tCon\n");
            for(int i=0; i<aFunc.tNum; i++){
                for(int j=0; j<aFunc.tDegree[i]; j++){
                    printf("%d", aFunc.tCon[i][j]);
                }
                printf("\n");
            }
        }
    }
}

void SVHSS(int repeatTime, hPara pa, BV_Para* bvPara, BV_PK* bvPk, BV_SK* bvSk){
    /**************prepareSharesIndex**************/
    int *tmpArr= malloc(sizeof(int)*pa.m);
    for(int i=0;i<pa.m;i++)
    {
        tmpArr[i]=i+1;
    }
    //combination before merging
    int *data = malloc(sizeof(int)*pa.t);
    pa.sIndex = (int **)(malloc(sizeof(int *) * pa.sNumPerVar));
    for(int i=0; i<pa.sNumPerVar; i++)
    {
        pa.sIndex[i]= ((int*) malloc(sizeof(int) * pa.t));
    }
    combinationUtil(tmpArr, pa.m, pa.t, 0, data, 0, pa.sIndex);

    //free memory
    free(tmpArr);
    tmpArr=NULL;
    free(data);
    data=NULL;
    /********************************************split********************************************/
    //[serverIndex][varIndex][shareIndex]
    share ***serShares = malloc(sizeof(share **) * pa.m);
    share ***alphaSerShares = malloc(sizeof(share **) * pa.m);
    for(int i=0; i<pa.m; i++)
    {
        serShares[i] = malloc(sizeof(share *) * pa.varNum);
        alphaSerShares[i] = malloc(sizeof(share *) * pa.varNum);
        for(int j=0; j<pa.varNum; j++)
        {
            serShares[i][j] = malloc(sizeof(share) * pa.sNumPerVar);
            alphaSerShares[i][j] = malloc(sizeof(share) * pa.sNumPerVar);
        }
    }

    fmpz_t *alphaVars = malloc(sizeof(fmpz_t) * pa.varNum);
    fmpz_t alpha;
    split(serShares, alphaSerShares, vars, alphaVars, alpha, pa, bvPara, bvPk);
    //direct computation
    fmpz_t directRes;
    fmpz_init(directRes);
    fmpz_zero(directRes);
    directCompute(directRes, aFunc, vars, bvPara);

    if (DEBUG){
        for(int i=0; i<pa.varNum;i++)
        {
            printf("vars[%d]:", i);
            fmpz_print(vars[i]);
            printf("\n");
        }
        printf("alpha:");
        fmpz_print(alpha);
        printf("\n");
    }

    //free memory
    free(vars);
    vars=NULL;
    free(alphaVars);
    alphaVars=NULL;
    /********************************************Compute********************************************/
    //Calculate all possible combinations ahead of time
    long long *allResSetNum = malloc(sizeof(long long)*(pa.d+1));
    int ***allResSet = malloc(sizeof(int **)*(pa.d+1));
    int **countRes = malloc(sizeof(int*)*(pa.d+1));//Count the number of occurrences of "share" and record the minimum number
    for (int de = 1; de <= pa.d ; de++) {
        allResSetNum[de] = pow(pa.sNumPerVar, de);//the number of share terms of degree de
        allResSet[de] = malloc(sizeof(int *)*allResSetNum[de]);
        countRes[de] = malloc(sizeof(int *)*allResSetNum[de]);
        for(long long i=0; i<allResSetNum[de]; i++)
        {
            allResSet[de][i] = malloc(sizeof(int)*de);
        }
    }
    shareCombine(allResSet, pa.d, pa);

    if(DEBUG){
        //check allResSet
        for (int de = 1; de <= pa.d ; de++) {
            printf("\n%lld terms combination\n", allResSetNum[de]);
            for(int i=0; i<allResSetNum[de]; i++)
            {
                for(int j=0; j<de; j++)
                {
                    printf("%d", allResSet[de][i][j]);
                }
                printf("\n");
            }
        }
    }

    //compute the state (encrypted or not encrypted)
    for (int de = 1; de <= pa.d ; de++) {
        for(long long i=0; i<allResSetNum[de];i++)
        {
            int sCnt[pa.m];
            for (int j = 0; j < pa.m; j++) {//virtual server index
                sCnt[j]=0;
                for(int k=0; k<de; k++)
                {
                    sCnt[j]=sCnt[j]+serShares[j][0][allResSet[de][i][k]].state;
                }
            }

            countRes[de][i] = findMinVIdx(sCnt, pa.m);

            if(DEBUG){
                //check sCnt
                for (int j = 0; j < pa.m; j++) {
                    printf("%d-%d,", j, sCnt[j]);
                }

                printf("minVIndx:%d\n", countRes[de][i]);
            }
        }
    }

    //alloc computation for the servers
    long long sTermsIndex = 0;
    sTerm* sTerms = malloc(sizeof(sTerm)*allSTermNum);
    for(int i=0; i<aFunc.tNum; i++)
    {
        //according degree to get the share term by spliting
        for(long long j=0; j<allResSetNum[aFunc.tDegree[i]]; j++)
        {
            sTerms[sTermsIndex].conLen = aFunc.tDegree[i];
            sTerms[sTermsIndex].symbol = aFunc.tSymbol[i];
            sTerms[sTermsIndex].tCon = malloc(sizeof(int)*sTerms[sTermsIndex].conLen);
            sTerms[sTermsIndex].sCon = malloc(sizeof(int)*sTerms[sTermsIndex].conLen);
            for(int k=0; k<sTerms[sTermsIndex].conLen; k++)
            {
                sTerms[sTermsIndex].tCon[k] =  aFunc.tCon[i][k];//the variable information
                sTerms[sTermsIndex].sCon[k] = allResSet[aFunc.tDegree[i]][j][k];//the share information
            }
            sTerms[sTermsIndex].ser = countRes[aFunc.tDegree[i]][j];//not real server index, the order of left servers
            sTermsIndex++;
        }
    }

    if(DEBUG){
        //check sTerms
        printf("allSTermNum:%lld\n", allSTermNum);
        for(long long i=0; i<allSTermNum;i++)
        {
            printf("conLen:%d\n", sTerms[i].conLen);
            printf("symbol:%d\n", sTerms[i].symbol);
            for (int j = 0; j < sTerms[i].conLen; j++)
            {
                printf("v%d-s%d  ", sTerms[i].tCon[j], sTerms[i].sCon[j]);
            }
            printf("===server:%d\n", sTerms[i].ser);
        }
    }

    //free momery
    for (int de = 1; de <= pa.d ; de++) {
        for(long long i=0; i<allResSetNum[de]; i++)
        {
            free(allResSet[de][i]);
            allResSet[de][i]=NULL;
        }
        free(allResSet[de]);
        allResSet[de]=NULL;
        free(countRes[de]);
        countRes[de]=NULL;
    }
    free(allResSet);
    allResSet=NULL;
    free(countRes);
    countRes=NULL;
    free(allResSetNum);
    allResSetNum=NULL;

    //Close different servers and compute
    clock_t start, end;
    double time = 0;
    cypher *res = malloc(sizeof(cypher)*pa.m);
    cypher *alphaRes = malloc(sizeof(cypher)*pa.m);

    for(int i=0; i<pa.m; i++){
        res[i] = (cypher)(malloc(sizeof(fmpz_mod_poly_t) * (pa.k+1)));
        alphaRes[i] = (cypher)(malloc(sizeof(fmpz_mod_poly_t) * (pa.k+1)));
        for(int j=0; j<=pa.k; j++){
            fmpz_mod_poly_init(res[i][j], bvPara->ctx_q);
            fmpz_mod_poly_init(alphaRes[i][j], bvPara->ctx_q);
        }
    }
    start = clock();
    for (int i = 0; i < repeatTime; i++) {
        compute(res, serShares, sTerms, allSTermNum, pa, bvPara, bvPk);
        alphaCompute(alphaRes, serShares, alphaSerShares, sTerms, allSTermNum, pa, bvPara, bvPk);
    }
    end = clock();
    time = (double) ((end - start)/repeatTime )/ CLOCKS_PER_SEC;
    printf("server running time: %f ms\n", time * 1000);

    //free memory
    free(sTerms);
    sTerms=NULL;
    for(int i=0; i<pa.m; i++)
    {
        for(int j=0; j<pa.varNum; j++)
        {
            free(serShares[i][j]);
            serShares[i][j]=NULL;

            free(alphaSerShares[i][j]);
            alphaSerShares[i][j]=NULL;
        }
        free(serShares[i]);
        serShares[i]=NULL;

        free(alphaSerShares[i]);
        alphaSerShares[i]=NULL;
    }
    free(serShares);
    serShares=NULL;

    free(alphaSerShares);
    alphaSerShares=NULL;

    /******************************************Dec and Verification*********************************/
    start = clock();
    for (int re = 0; re < repeatTime; re++) {
        cypher fRes, fAlphaRes;

        fRes = (cypher)(malloc(sizeof(fmpz_mod_poly_t) * (pa.k+1)));
        fAlphaRes = (cypher)(malloc(sizeof(fmpz_mod_poly_t) * (pa.k+1)));
        for(int i=0; i<=pa.k; i++){
            fmpz_mod_poly_init(fRes[i], bvPara->ctx_q);
            fmpz_mod_poly_zero(fRes[i], bvPara->ctx_q);

            fmpz_mod_poly_init(fAlphaRes[i], bvPara->ctx_q);
            fmpz_mod_poly_zero(fAlphaRes[i], bvPara->ctx_q);
        }
        for(int i=0; i<pa.m; i++){
            BV_Add(fRes, fRes, res[i], pa.k, bvPara);
            BV_Add(fAlphaRes, fAlphaRes, alphaRes[i], pa.k, bvPara);
        }

        fmpz_t decRes, alphaDecRes, tmp;
        fmpz_init(decRes);
        fmpz_init(alphaDecRes);

        BV_Dec(decRes, fRes, pa.k, bvPara, bvSk);
        BV_Dec(alphaDecRes, fAlphaRes, pa.k, bvPara, bvSk);

        fmpz_mod(decRes, decRes, bvPara->msg);
        fmpz_mod(alphaDecRes, alphaDecRes, bvPara->msg);

        fmpz_mod_mul(tmp, alpha, decRes, bvPara->ctx_msg);

        if (DEBUG){
            //check tmp
            printf("decRes:\n");
            fmpz_print(decRes);
            printf("\n");
            printf("tmp:\n");
            fmpz_print(tmp);
            printf("\n");
            printf("alphaDecRes:\n");
            fmpz_print(alphaDecRes);
            printf("\n");
        }

        if (fmpz_equal(tmp, alphaDecRes)){
            printf("Verification Successes!!!\n");
        }else{
            printf("Verification Fails!!\n");
        }
    }
    end = clock();
    time = (double) ((end - start)/repeatTime) / CLOCKS_PER_SEC;
    printf("client running time: %f ms\n", time * 1000);
    if (DEBUG){
        printf("directRes:");
        fmpz_print(directRes);
        printf("\n");
    }
}

