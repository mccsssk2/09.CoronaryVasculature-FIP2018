#include <stdio.h>
#include <stdlib.h>
void doit(int ** s , int row , int col ) 
{
    for(int i =0 ; i <col;i++){
            for(int j =0 ; j <col;j++)
                printf("%d ",s[i][j]);
            printf("\n");
    }
}
int main()
{
    int row =10 , col=10;
    int ** c = (int**)malloc(sizeof(int*)*row);
    for(int i =0 ; i <col;i++)
        *(c+i) = (int*)malloc(sizeof(int)*row);
    for(int i =0 ; i <col;i++)
            for(int j =0 ; j <col;j++)
                c[i][j]= 1; // i*j;
    doit(c,row,col);
}
