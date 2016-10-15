//
//  ElectroMagneticWave.c
//  
//
//  Created by 原口俊樹 on 2016/09/08.
//
//



#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>

#define NX 230               /* 空間セル数 X [pixels] */
#define NY 102				/* 空間セル数 Y [pixels] */
#define NZ 200            /* 空間セル数 Z [pixels] */

#define dx 0.0001			/* 空間刻み [m] */
#define dy 0.0001				/* 空間刻み [m] */
#define dz 100e-6				/* 空間刻み [m] */
#define dt 1.92e-13		/* 時間刻み [s] */

#define Nstep 300			/* 計算ステップ数 [step] */
#define Fileoutstep 10		/* 音圧分布ファイル出力ステップ数 [step] */

#define freq 9.4e9			/* 初期波形の周波数 [Hz] */

#define ips 8.854e-12				/* 空気中の誘電率 */
#define mu 4*M_PI*1e-7                /* 空気中の透磁率 */
#define sigma   0                  /*導電率*/
#define sigma2 0                   /*磁気導電率*/
#define v 418616961.7
#define a 22.86e-3
#define c 3e8


double Ex[NX+1][NY+1][NZ+1];		/* x方向電界 [V/m] */
double Ey[NX+1][NY+1][NZ+1];		/* y方向電界 [V/m] */
double Ez[NX+1][NY+1][NZ+1];        /* z方向電界 [V/m] */

double oldEy0[NX][NY];
double oldEy1[NX][NY];
double oldEyNZ1[NX][NY];
double oldEyNZ2[NX][NY];

double Hx[NX+1][NY+1][NZ+1];			/* x方向磁界 [A/m] */
double Hy[NX+1][NY+1][NZ+1];			/* y方向磁界 [A/m] */
double Hz[NX+1][NY+1][NZ+1];			/* z方向磁界 [A/m] */
double ip[NX][NY][NZ];

double seedips = 90.0;

void UpdateE(),UpdateH(),Mur();

double B = ((mu/dt-sigma2/2)/(mu/dt+sigma2/2));
float Amp = 1;


/* メイン関数 ************************************************************/

int main(void)
{
    
   
    int i,j,k,m,l,AAA;
    int n;
    double sig;
    FILE *waveformfile, *fieldfile;
    char fieldfilename[30];
    time_t t_stert,t_end;
    n=0;
    
    /* 事前準備 *********************************************************/
    if((waveformfile = fopen("waveform.txt","w"))==NULL){
        printf("open error [waveform.txt]\n");
        return(1);
    }
    
    /* 粒子速度分布・音圧分布を初期化 *********************************************************/
    
    



    
    for(i=0;i<NX;i++){
        for(j=0;j<NY;j++){
            for(k=0;k<NZ;k++){
                Ex[i][j][k] = 0.0;
                Ey[i][j][k] = 0.0;
                Ez[i][j][k] = 0.0;
                Hx[i][j][k] = 0.0;
                Hy[i][j][k] = 0.0;
                Hz[i][j][k] = 0.0;
                ip[i][j][k] = ips;
                
            }
            oldEy0[i][j] = 0;
            oldEy1[i][j] = 0;
            }
    }
   /* for(i=109;i<119;i++){
        for(j=0;j<10;j++){
            for(k=495;k<505;k++){
                ip[i][j][k] = seedips * ips;
            }
        }
    }*/
    

     t_stert = time(NULL);
    /* メインループ *********************************************************/
    while(n<=Nstep){
        
        /* 更新（ここが FDTD の本丸!!） */
       
        UpdateE();

        UpdateH();
        
     /* 初期波形を準備（正弦波×１波 with ハン窓）*/
        if(n  < 10)
            for(i=1;i<NX-1;i++)
                for(j=1;j<NY-1;j++)
                    Ey[i][j][100] = Amp*sin((3.14/a)*i*dx)*sin(2*3.14*freq*dt*n) + Ey[i][j][100];


        /* 音源 */
    
        /* 波形ファイル出力（時刻, 音源, 中央点の音圧） */
        fprintf(waveformfile,"%e\t%e\t%e\n", dt*n, sig, Ey[(int)(50)][(int)(50)][(int)(50)]);
        fprintf(stderr,"%5d / %5d\r", n, Nstep );
        
        /* 音圧分布ファイル出力 */
        if( n % Fileoutstep == 0 ){
            sprintf(fieldfilename, "field%.6d.txt",n);
            if((fieldfile = fopen(fieldfilename,"w"))==NULL){
                printf("open error [field.txt]\n");
                return(1);
            }
            for(i=0; i<NX; i++){
                for(j=0; j<NZ; j++){
                   /* if(Ey[i][50][j]>-1 && Ey[i][50][j]<1){
                        Ey[i][50][j] = 0.0;
                    }*/
                    fprintf(fieldfile, "%e\t",Ey[i][50][j]);
                }
                fprintf(fieldfile, "\n");
            }
            fclose(fieldfile);
        }
        n++;
    }
    /* 事後処理 *********************************************************/
    fclose(waveformfile);
    t_end = time(NULL);
    printf("Succeed! %d \n",(int)(t_end - t_stert));
    return(0);
}


/* サブルーチン ************************************************************/

void UpdateE()
{
    int i,j,k;
    
    #pragma omp parallel for private(i,j,k)
    for(i=0;i<NX;i++){
        for(j=1;j<NY;j++){
            for(k=1;k<NZ;k++){
                Ex[i][j][k] = ((ip[i][j][k]/dt-sigma/2)/(ip[i][j][k]/dt+sigma/2))*Ex[i][j][k] + (Hz[i][j][k]-Hz[i][j-1][k])/(dy*(ip[i][j][k]/dt+sigma/2)) - (Hy[i][j][k]-Hy[i][j][k-1])/(dz*(ip[i][j][k]/dt+sigma/2));
                
            }
        }
    }

    #pragma omp parallel for private(i,j,k)
    
    for(i=1;i<NX;i++){
        for(j=0;j<NY-1;j++){
            for(k=1;k<NZ;k++){
                Ey[i][j][k] = ((ip[i][j][k]/dt-sigma/2)/(ip[i][j][k]/dt+sigma/2))*Ey[i][j][k] + (Hx[i][j][k]-Hx[i][j][k-1])/(dz*(ip[i][j][k]/dt+sigma/2)) - (Hz[i][j][k]-Hz[i-1][j][k])/(dx*(ip[i][j][k]/dt+sigma/2));

                    }
            }
    }

    #pragma omp parallel for private(i,j,k)
    
    for(i=1;i<NX;i++){
        for(j=1;j<NY;j++){
            for(k=0;k<NZ-1;k++){
                Ez[i][j][k] = ((ip[i][j][k]/dt-sigma/2)/(ip[i][j][k]/dt+sigma/2))*Ez[i][j][k] + (Hy[i][j][k]-Hy[i-1][j][k])/(dx*(ip[i][j][k]/dt+sigma/2)) - (Hx[i][j][k]-Hx[i][j-1][k])/(dy*(ip[i][j][k]/dt+sigma/2));
            }
        }
    }
       Mur();
}





void Mur()
{
    int i,j,k;
    double k0 = 2 * M_PI * freq * sqrt(mu * ips);

    
    for(i=1;i<NX-1;i++)
        for(j=0;j<NY-1;j++){
                Ey[i][j][0] = oldEy1[i][j] - ((dz - v*dt)/(dz + v*dt))*(Ey[i][j][1] - oldEy0[i][j]);
                Ey[i][j][NZ-1] = oldEyNZ2[i][j] - ((dz - v*dt)/(dz + v*dt))*(Ey[i][j][NZ-2] - oldEyNZ1[i][j]);
            }
    for(i=1;i<NX-1;i++){
        for(j=0;j<NY-1;j++){
                oldEy0[i][j] = Ey[i][j][0];
                oldEy1[i][j] = Ey[i][j][1];
                oldEyNZ1[i][j] = Ey[i][j][NZ-1];
                oldEyNZ2[i][j] = Ey[i][j][NZ-2];
            }
    }
    
}



/* 音圧の更新 */
void UpdateH()
{
    int i,j,k;
    
    #pragma omp parallel for private(i,j,k)
    
    for(i=0;i<NX;i++){
        for(j=0;j<NY-1;j++){
            for(k=0;k<NZ-1;k++){
                Hx[i][j][k] = B*Hx[i][j][k] - (Ez[i][j+1][k]-Ez[i][j][k])/(dy*(mu/dt+sigma2/2)) + (Ey[i][j][k+1]-Ey[i][j][k])/(dz*(mu/dt+sigma2/2));
            }
        }
    }
    
    #pragma omp parallel for private(i,j,k)
    
    for(i=0;i<NX-1;i++){
        for(j=0;j<NY;j++){
            for(k=0;k<NZ-1;k++){
                Hy[i][j][k] = B*Hy[i][j][k] - (Ex[i][j][k+1]-Ex[i][j][k])/(dz*(mu/dt+sigma2/2)) + (Ez[i+1][j][k]-Ez[i][j][k])/(dx*(mu/dt+sigma2/2));
            }
        }
    }
    
    #pragma omp parallel for private(i,j,k)
    
    
    for(i=0;i<NX-1;i++){
        for(j=0;j<NY-1;j++){
            for(k=0;k<NZ;k++){
                Hz[i][j][k] = B*Hz[i][j][k] - (Ey[i+1][j][k]-Ey[i][j][k])/(dx*(mu/dt+sigma2/2)) + (Ex[i][j+1][k]-Ex[i][j][k])/(dy*(mu/dt+sigma2/2));
            }
        }
    }
    


}
