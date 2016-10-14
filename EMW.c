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

#define Nstep 200			/* 計算ステップ数 [step] */
#define Fileoutstep 10		/* 音圧分布ファイル出力ステップ数 [step] */

#define freq 9.4e9			/* 初期波形の周波数 [Hz] */

#define ips 8.854e-12				/* 空気中の誘電率 */
#define mu 4*M_PI*1e-7                /* 空気中の透磁率 */
#define sigma   0                  /*導電率*/
#define sigma2 0                   /*磁気導電率*/
#define v 418616961.7
#define a 22.86e-3
#define Mode 2



double Ex[NX][NY][NZ];		/* x方向電界 [V/m] */
double Ey[NX][NY][NZ];		/* y方向電界 [V/m] */
double Ez[NX][NY][NZ];        /* z方向電界 [V/m] */
double oldEx[NX][NY][4];
double oldEy[NX][NY][4];
double oldEz[NX][NY][4];

double Hx[NX][NY][NZ];			/* x方向磁界 [A/m] */
double Hy[NX][NY][NZ];			/* y方向磁界 [A/m] */
double Hz[NX][NY][NZ];			/* z方向磁界 [A/m] */
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
            oldEx[i][j][0] = 0;
            oldEy[i][j][0] = 0;
            oldEz[i][j][0] = 0;
            oldEx[i][j][1] = 0;
            oldEy[i][j][1] = 0;
            oldEz[i][j][1] = 0;
            oldEx[i][j][2] = 0;
            oldEy[i][j][2] = 0;
            oldEz[i][j][2] = 0;
            oldEx[i][j][3] = 0;
            oldEy[i][j][3] = 0;
            oldEz[i][j][3] = 0;
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
    
    for(i=1;i<NX-1;i++){
        for(j=1;j<NY-1;j++){
            for(k=1;k<NZ-1;k++){
                    if(j==1 || j == NY-2 )
                        Ex[i][j][k] = 0;
                    else
                        Ex[i][j][k] = ((ip[i][j][k]/dt-sigma/2)/(ip[i][j][k]/dt+sigma/2))*Ex[i][j][k] + (Hz[i][j+1][k]-Hz[i][j-1][k])/(dy*(ip[i][j][k]/dt+sigma/2)) - (Hy[i][j][k+1]-Hy[i][j][k-1])/(dz*(ip[i][j][k]/dt+sigma/2));
                
                    //printf("Ex[%d][%d][%d] = %5.2lf\n",i,j,k,Ex[i][j][k]);
            }
        }
    }

    #pragma omp parallel for private(i,j,k)
    
    for(i=1;i<NX-1;i++){
        for(j=1;j<NY-1;j++){
            for(k=1;k<NZ-1;k++){
                    if(i==1 || i== NX-2)
                        Ey[i][j][k] = 0;
                    else
                        Ey[i][j][k] = ((ip[i][j][k]/dt-sigma/2)/(ip[i][j][k]/dt+sigma/2))*Ey[i][j][k] + (Hx[i][j][k+1]-Hx[i][j][k-1])/(dz*(ip[i][j][k]/dt+sigma/2)) - (Hz[i+1][j][k]-Hz[i-1][j][k])/(dx*(ip[i][j][k]/dt+sigma/2));
                    }
            }
    }

    #pragma omp parallel for private(i,j,k)
    
    
        //  printf("Eyの演算完了！\n");
    for(i=1;i<NX-1;i++){
        for(j=1;j<NY-1;j++){
            for(k=1;k<NZ-1;k++){
                if(i==1 || i== NX-2 || j==1 || j == NY-2 )
                    Ez[i][j][k] = 0;
                else
                    Ez[i][j][k] = ((ip[i][j][k]/dt-sigma/2)/(ip[i][j][k]/dt+sigma/2))*Ez[i][j][k] + (Hy[i+1][j][k]-Hy[i-1][j][k])/(dx*(ip[i][j][k]/dt+sigma/2)) - (Hx[i][j+1][k]-Hx[i][j-1][k])/(dy*(ip[i][j][k]/dt+sigma/2));
                
                                          //  printf("Ez[%d][%d][%d] = %5.2lf\n",i,j,k,Ez[i][j][k]);
            }
        }
    }
  //  printf("Ezの演算完了！\n");
       Mur();
}

void Mur()
{
 double k0 = 2 * M_PI * freq * sqrt(mu * ips);
}



/* 音圧の更新 */
void UpdateH()
{
    int i,j,k;
    
    #pragma omp parallel for private(i,j,k)
    
    for(i=1;i<NX-1;i++){
        for(j=1;j<NY-1;j++){
            for(k=1;k<NZ-1;k++){
                Hx[i][j][k] = B*Hx[i][j][k] - (Ez[i][j+1][k]-Ez[i][j-1][k])/(dy*(mu/dt+sigma2/2)) + (Ey[i][j][k+1]-Ey[i][j][k-1])/(dz*(mu/dt+sigma2/2));
                            //  printf("Hx[%d][%d][%d] = %5.2lf\n",i,j,k,Hx[i][j][k]);
            }
        }
    }
    
    #pragma omp parallel for private(i,j,k)
    
    for(i=1;i<NX-1;i++){
        for(j=1;j<NY-1;j++){
            for(k=1;k<NZ-1;k++){
                Hy[i][j][k] = B*Hy[i][j][k] - (Ex[i][j][k+1]-Ex[i][j][k-1])/(dz*(mu/dt+sigma2/2)) + (Ez[i+1][j][k]-Ez[i-1][j][k])/(dx*(mu/dt+sigma2/2));
                            //   printf("Hy[%d][%d][%d] = %5.2lf\n",i,j,k,Hy[i][j][k]);
            }
        }
    }
    
    #pragma omp parallel for private(i,j,k)
    
    
    for(i=1;i<NX-1;i++){
        for(j=1;j<NY-1;j++){
            for(k=1;k<NZ-1;k++){
                Hz[i][j][k] = B*Hz[i][j][k] - (Ey[i+1][j][k]-Ey[i-1][j][k])/(dx*(mu/dt+sigma2/2)) + (Ex[i][j+1][k]-Ex[i][j-1][k])/(dy*(mu/dt+sigma2/2));
                              //  printf("Hz[%d][%d][%d] = %5.2lf\n",i,j,k,Hz[i][j][k]);
            }
        }
    }
    


}
