#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define INFINITE 99999

    double *costCoeff; // Cost Coefficient Matrix
    double **ConstraintCoeff;// Constraint Coefficient Matrix
    double *bValues;// B Matrix


    double *CVector;// C Vector
    double **ABVector;// AB Combined Matrix

    int numberOfCostCoeff;// No. of Cost Coefficients
    int numberOfConstraints;// No. of Constraints imposed
    int numberOfVariables;// No. of Variables

    // To keep track of basic variables.... 0 for basic..... 1 for non-basic
    int *basicIdentifier;

    // To check whether values obtained from current iteration is optimal or not....
    int checkForOptimality()
    {
        for(int i = 0; i < numberOfVariables; i++)
        {
            if(CVector[i] < 0)
                return 1;

        }
            return 0;
    }

    // To get the index of the entering variable.....
    int getEnteringIndex()
    {
        int minIndex = 0;

        for(int i = 1; i < numberOfVariables; i++)
        {
            if(CVector[minIndex] > CVector[i])
                minIndex = i;
        }
        return minIndex;//  Min. and negative index
    }

int main()
{

    int i,j;
    basicIdentifier = (int*)malloc(numberOfVariables*sizeof(int));

    printf("**********ENTER DATA**********\n");
    printf("Enter no. of Cost Coefficients : ");
    scanf("%d",&numberOfCostCoeff);
    printf("\n");

    costCoeff = (double*)malloc(numberOfCostCoeff * sizeof(double));

    for( i = 0; i < numberOfCostCoeff; i++){
        printf("ENTER COST COEFFICIENT c[%d] : ",i);
        scanf("%lf",&costCoeff[i]);
    }

    printf("\n\n");

    printf("Enter no. of Constraints : ");
    scanf("%d",&numberOfConstraints);
    printf("\n");

    bValues = (double*)malloc(numberOfConstraints * sizeof(double));

    for(i = 0; i < numberOfConstraints; i++){
        printf("ENTER B VALUE b[%d] : ",i);
        scanf("%lf",&bValues[i]);
    }
    printf("\n\n");


    numberOfVariables = numberOfConstraints + numberOfCostCoeff;

    ConstraintCoeff = (double**)malloc(numberOfConstraints * sizeof(double*));

    for(i = 0; i < numberOfConstraints; i++)
    {
        ConstraintCoeff[i] = (double*)malloc(numberOfVariables * sizeof(double));
    }

    for(i = 0; i < numberOfConstraints; i++)
    {
        for(j = 0; j < numberOfVariables; j++)
        {
            printf("ENTER CONSTRAINT COEFFICIENTS a[%d][%d] : ",i,j);
            scanf("%lf",&ConstraintCoeff[i][j]);
        }
    }
    printf("\n\n");


    printf("***********LPP IN STANDARD FORM***********\n");
    printf("MINIMIZATION FOR Z : \n\n");
    printf("\t");
    for(i = 0; i < numberOfCostCoeff; i++)
    {
        if(i != numberOfCostCoeff - 1)
            printf("(%.2lf)X%d +",costCoeff[i],i+1);
        else
            printf("(%.2lf)X%d\n",costCoeff[i],i+1);
    }
    printf("\n");
    printf("Constraints imposed on Z : \n");
    printf("\n\n");
    for(i = 0; i < numberOfConstraints; i++)
    {
        printf("\t");
        for(j = 0; j <= numberOfVariables; j++)
        {
            if(j == numberOfVariables)
                printf(" = %.2lf\n",bValues[i]);
            else
            {
                if(j < numberOfVariables - 1)
                    printf("(%.2lf)X%d +",ConstraintCoeff[i][j],j+1);
                else
                    printf("(%.2lf)X%d",ConstraintCoeff[i][j],j+1);
            }
        }
    }

    CVector = (double*)malloc(numberOfVariables*sizeof(double));
    for( i = 0; i < numberOfVariables; i++)
    {
        if(i < numberOfCostCoeff)
            CVector[i] = costCoeff[i];
        else
            CVector[i] = 0;
    }

    ABVector = (double**)malloc(numberOfConstraints*sizeof(double*));
    for(i = 0; i < numberOfConstraints; i++)
    {
        ABVector[i] = (double*)malloc((numberOfVariables+1)*sizeof(double));
    }

    for(i = 0; i < numberOfConstraints; i++)
    {
        for( j = 0; j < numberOfVariables+1; j++)
        {
            if(j != numberOfVariables)
                ABVector[i][j] = ConstraintCoeff[i][j];
            else
                ABVector[i][j] = bValues[i];
        }
    }

    // Solution Array
    double *solution = (double*)malloc(numberOfVariables * sizeof(double));

    // Basic Feasible Solution
    for( i = 0; i < numberOfVariables; i++)
    {
        if(i < numberOfCostCoeff){
            solution[i] = 0;
            basicIdentifier[i] = 0;
        }
        else{
            solution[i] = bValues[i - numberOfCostCoeff];
            basicIdentifier[i] = 1;
        }
    }

    double Min = 0;


    int enteringBasisVariableIndex = -1;

    int steps = 0;

    printf("------------------------------------------------------------------\n");
    printf("\n\nSTEP : 0\n\n");

    printf("\n\nBASIC:\n");
        for(i = 0; i < numberOfVariables; i++)
        {
            printf("%d ",basicIdentifier[i]);
        }
        printf("\n");

    printf("\nABVector : \n");
        for( i = 0; i < numberOfConstraints; i++){
            for(j = 0; j <= numberOfVariables; j++)
            {
                printf("%.2lf ",ABVector[i][j]);
            }
            printf("\n");
        }

    printf("CVector : \n");
        for(j = 0; j < numberOfVariables; j++)
            {
                printf("%.2lf ",CVector[j]);
            }
            printf("\n");


    printf("------------------------------------------------------------------\n");

    while(checkForOptimality())
    {

        printf("------------------------------------------------------------------\n");
        printf("\n\nSTEP : %d\n",++steps);
        printf("\n");
        enteringBasisVariableIndex = getEnteringIndex();  // pivot element column

        double *tempArray;
        tempArray = (double*)malloc(numberOfConstraints*sizeof(double));

        for(i = 0; i < numberOfConstraints; i++)
        {
            if(ABVector[i][enteringBasisVariableIndex] <= 0)
                tempArray[i] = INFINITE;
            else
                tempArray[i] = ABVector[i][numberOfVariables]/ABVector[i][enteringBasisVariableIndex];
        }

        int acount = 0;
        for(i = 0; i < numberOfConstraints; i++)
        {
            if(tempArray[i] == INFINITE)
                acount++;
        }


        int minTempIndex = 0;
        for( i = 1; i < numberOfConstraints; i++)
        {
            if(tempArray[minTempIndex] > tempArray[i])
                minTempIndex = i;
        }


        double pivot = ABVector[minTempIndex][enteringBasisVariableIndex];

        printf("%.2lf is pivot\n",pivot);

        for(i = 0; i < numberOfVariables; i++)
        {
            if(ABVector[minTempIndex][i] == 1 && basicIdentifier[i] == 1)
            {
                basicIdentifier[i] = 0;
                basicIdentifier[enteringBasisVariableIndex] = 1;
            }
        }

        printf("\n\nBASIC:\n");
        for(i = 0; i < numberOfVariables; i++)
        {
            printf("%d ",basicIdentifier[i]);
        }
        printf("\n");

        for(i = 0; i <= numberOfVariables; i++)
            ABVector[minTempIndex][i] = ABVector[minTempIndex][i]/pivot;

        double t=0;

        for(i = 0; i < numberOfConstraints; i++)
        {
            t = ABVector[i][enteringBasisVariableIndex];
            if(i == minTempIndex)
                    continue;
            for(j = 0; j <= numberOfVariables; j++)
            {
                if(t < 0)
                    ABVector[i][j] = ABVector[i][j] + ABVector[minTempIndex][j]*fabs(t);
                else
                    ABVector[i][j] = ABVector[i][j] - ABVector[minTempIndex][j]*t;
            }
        }

        printf("\nABVector : \n");
        for( i = 0; i < numberOfConstraints; i++){
            for(j = 0; j <= numberOfVariables; j++)
            {
                printf("%.2lf ",ABVector[i][j]);
            }
            printf("\n");
        }

        t = CVector[enteringBasisVariableIndex];
        for(j = 0; j < numberOfVariables; j++)
            {
                if(t < 0)
                    CVector[j] = CVector[j] + ABVector[minTempIndex][j]*fabs(t);
                else
                    CVector[j] = CVector[j] - ABVector[minTempIndex][j]*t;
            }

            printf("CVector : \n");
        for(j = 0; j < numberOfVariables; j++)
            {
                printf("%.2lf ",CVector[j]);
            }
            printf("\n");


    printf("------------------------------------------------------------------\n");

            if(acount == numberOfConstraints)
           {
            printf("Solution is UNBOUNDED.\n");

            Min = 0;
            for(i = 0; i < numberOfCostCoeff; i++)
                Min = Min + solution[i]*costCoeff[i];

            printf("MINIMUM VALUE : %.2lf",Min);

            return 0;
        }

    }

    printf("\n\nSOLUTION:\n");
    for( i = 0; i < numberOfConstraints; i++){
        for(j = 0; j < numberOfVariables; j++)
        {
            if(ABVector[i][j] == 1 && basicIdentifier[j] == 1)
            {
                solution[j] = ABVector[i][numberOfVariables];
            }
        }
    }

    for(j = 0; j < numberOfVariables; j++)
    {
        if(basicIdentifier[j] == 0)
            solution[j] = 0;
    }

    for(j = 0; j < numberOfCostCoeff; j++)
        printf("%.2lf ",solution[j]);
    printf("\n\n");

    Min = 0;
    for(i = 0; i < numberOfCostCoeff; i++)
        Min = Min + solution[i]*costCoeff[i];

    printf("MINIMUM VALUE : %.2lf\n",Min);

    return 0;
}
