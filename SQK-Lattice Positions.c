#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{   
    int i, j, L;        
    float lc ; // lc -> lattice constant
    FILE *file01;

    file01 = fopen("/home/sumank/project codes/square_kagome_positions.txt", "w+");
    printf("Enter the lattice dimension : ");
    scanf("%d", &L);

    float posA[L][L][2], posB[L][L][2], posC[L][L][2], posD[L][L][2], posE[L][L][2], posF[L][L][2];
    
    lc = 1.0 + sqrt(3.0);
    
    for(j=0; j<L; j++)
    {
        for(i=0; i<L; i++)                     // Note the i and j order, it is such to match the indexing of Newman book
        {
            posA[i][j][0] = (float)i*lc + 0.5*sqrt(3.0);
            posA[i][j][1] = (float)j*lc + 0.5*sqrt(3.0);
            
            posB[i][j][0] = (float)i*lc + 0.5*sqrt(3.0) + 1.0;
            posB[i][j][1] = (float)j*lc + 0.5*sqrt(3.0);
            
            posC[i][j][0] = (float)i*lc + 0.5*sqrt(3.0) + 1.0;
            posC[i][j][1] = (float)j*lc + 0.5*sqrt(3.0) + 1.0;
            
            posD[i][j][0] = (float)i*lc + 0.5*sqrt(3.0);
            posD[i][j][1] = (float)j*lc + 0.5*sqrt(3.0) + 1.0;
            
            posE[i][j][0] = (float)i*lc + 0.5*sqrt(3.0) + 0.5;
            posE[i][j][1] = (float)j*lc + 0.0;
            
            posF[i][j][0] = (float)i*lc + sqrt(3.0) + 1.0;
            posF[i][j][1] = (float)j*lc + 0.5*sqrt(3.0) + 0.5;

            fprintf(file01, " %f \t %f \n", posA[i][j][0], posA[i][j][1]);
            fprintf(file01, " %f \t %f \n", posB[i][j][0], posB[i][j][1]);
            fprintf(file01, " %f \t %f \n", posC[i][j][0], posC[i][j][1]);
            fprintf(file01, " %f \t %f \n", posD[i][j][0], posD[i][j][1]);
            fprintf(file01, " %f \t %f \n", posE[i][j][0], posE[i][j][1]);
            fprintf(file01, " %f \t %f \n", posF[i][j][0], posF[i][j][1]);
        }

    }
    fclose(file01);
    return 0;
}
