#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PI 3.1415
#define direc "E:\\workstation\\workstation\\sumProductDecoding\\H.txt"
#define NUM_1 12
#define NUM_2 20
#define NUM_3 200
#define NUM_4 330
#define CODE_LEN 155
#define SIGMA 1
#define B_SIZE 3
#define A_SIZE 5
#define CHECK_SIZE 93
#define MIN_SNR 2.0
#define MAX_SNR 4.0
#define ERR_SUM 30.0

const int m[CODE_LEN] = {0} ;
double m_changed[CODE_LEN] ;
double r[CODE_LEN];
double Initi_M[CODE_LEN][CHECK_SIZE];//Initiazation
int Z[CODE_LEN];//Test_z

int N , M ;      /* size of the matrix */
int biggest_num_n,biggest_num_m;
int B[CODE_LEN][B_SIZE],A[CHECK_SIZE][A_SIZE];
int H[CODE_LEN][CHECK_SIZE]={0};
double E[CHECK_SIZE][CODE_LEN] = {0};//矩阵，CheckMessage

//matrix = (alist_matrix *)malloc(sizeof(alist_matrix));

void write_alist();
double Gauss();
double e(int *B,double *);
double** CheckMessage(int *a,double **e,int num);
//int* Test_z(**E,*r,**A,N);

int main()
{
    int CNT = 0;//信噪比取值计数
    int cnt[30] = {0};//每一个信噪比对应的总数
    double snr = MIN_SNR,gauss,stdev;
	write_alist();
	int i,j;
/*    for(int i = 0;i < M;i++)
	{
		for(int j = 0;j < N;j++)
			printf("%d ",H[i][j]);
		printf("\n");
	}
	*/
    while(snr <= MAX_SNR)
    {
        int I = 0;//取30个错误
		int Err = 0;
        long sum_Single = 0;//每30次错误运行结束后的总运行次数
		//根据snr 计算标准差  snr=10log10(1/2*r*stdev^2)
		stdev=sqrt(1.0/((2.0*(N - M)/N)*pow(10,snr/10)));
		printf("stdev is %lf\n",stdev);
        while(Err < 30)
        {
			for(i = 0;i < CODE_LEN;i++)//通过高斯信道得到码字m_changed,以及Initialization中的R
			{
				gauss = stdev * Gauss();
				m_changed[i] = m[i] + gauss +1;
				r[i] = 4 * m_changed[i] * snr;
			}
			for(i = 0;i < N;i++)//Initialization中的M
            {
                for(int j = 0;j < M;j++)
                    Initi_M[j][i] = r[i];
            }

			/*///printf Initialization中的M
			printf("Initi_M is\n");
			for(i = 0;i < M;i++)
            {
                for(int j = 0;j < N;j++)
                    printf("%.4lf ",Initi_M[i][j]);
				printf("\n");
            }*/
			//----------------
			int flag = 0;//标志有错时一直执行到该次译码成功
			while(!flag)
			{
				for(j = 0;j < M;j++)//CheckMessage中的E
				{
					for(i = 0;i < B_SIZE;i++)
					{
						double t = 1.0;
						for(int z = 0;z < B_SIZE;z++)
						{
							if(i == z)continue;
							t *= tanh(Initi_M[j][B[j][z]]/2);
						}
						E[j][B[j][i]] = log((1 + t)/(1 - t));
					}
				}
				/*/printf E
				printf("test E\n");
				for(int i = 0;i < M/2;i++)//CheckMessage中的E
				{
					for(int j = 0;j < N/2;j++)
					{
						printf("%lf ",E[i][j]);
					}
					printf("\n");
				}*/
				//------------------------
				for(i = 0;i < N;i++)//求得Test中Z
				{
					double t = 0.0;
					for(j = 0;j < A_SIZE;j++)
						t += E[A[i][j]][i];
					Z[i] = (t + r[i] <= 0)?1:0;
				}

				//printf Z
				printf("\nZ[i]\n");
				for(i = 0;i < N;i++)
				{
					printf("%d ",Z[i]);
				}
				//------------------------
				int tmp_Err = 0;
				for(i = 0;i < N;i++)
				{
					if(Z[i] != 0)
					{
						tmp_Err++;
					}
				}
				if(tmp_Err != 0)
				{
					Err++;
					if(Err == 30)break;
				}
				else flag = 1;
				if(!flag)
				{
					for(i = 0;i < N;i++)
					{
						for(int z = 0;z < A_SIZE;z++)
						{
							double t = 0.0;
							for(j = 0;j < A_SIZE;j++)
							{
								if(j == z)continue;
								t += E[A[i][j]][i];
							}
							Initi_M[A[i][z]][i] = t + r[i];
						}
					}
					/*///printf Initialization中的M
					printf("Initi_M is\n");
					for(i = 0;i < M;i++)
					{
						for(int j = 0;j < N;j++)
							printf("%.4lf ",Initi_M[i][j]);
						printf("\n");
					 }*/
				//----------------
				}
				sum_Single ++;
			}
        }
		double err_rate = ERR_SUM/sum_Single;
        printf("When snr = %lf,errorRate = %lf\n",snr,err_rate);
		CNT++;
        snr += 0.1;
    }
    return 0;
}


/*int* Test_z(**E,*r,**A,int N)
{
    int t = 0,*L;
    int i,j;
    for(i = 0;i < N;i++)
    {
        for(j = 0;j < A_SIZE;j++)
        {
            t += E[A[i][j]][i];
        }
        *(L + i) = (int)malloc(sizeof(int));
        *(L + i) = t + *(r + i);
        *(L + i) = (*(L + i) <= 0?1:0);
    }
    return L;
}
*/
/*double** CheckMessage(int *a,double **e,int num)
{
    double ret[CODE_LEN][CODE_LEN];
    int i,j;
    for(j = 0;j < num;j++)
    {
        double x = 1.0;
        for(;i = *(a + j);i = a[j++])
        {
            x *= tanh(e[j][i] / 2);
        }
        ret[j][i] = log((1 + x)/(1-x));
    }
    return **ret;
}
*/
void write_alist ()
 {
  /* this assumes that A and B have the form of a rectangular
     matrix in the file; if lists have unequal lengths, then the
     entries should be present (eg zero values) but are ignored
     */
    FILE *fp;
    fp = fopen(direc,"r");
    if(!fp)
    {
        printf("open file failed!");
        return;
    }
    char tmp[NUM_4];
	memset(tmp,0,sizeof(tmp));
    fgets(tmp,NUM_1,fp);
    int i=0,t = 10;
    M = N = 0;
    while(tmp[i] != ' ')
    {
        N = N * t + (tmp[i] - '0');
        ++i;
    }
    ++i;
    while((tmp[i] != ' ')&&(tmp[i] != '\n'))
    {
        M = M * t + (tmp[i] - '0');
        ++i;
    }
    memset(tmp,0,sizeof(tmp));
    fgets(tmp,NUM_1,fp);
    i = 0;
	biggest_num_n = biggest_num_m = 0;
    while(tmp[i] != ' ')
    {
        biggest_num_n = biggest_num_n * t + (tmp[i] - '0');
        ++i;
    }
    ++i;
    while(tmp[i] != ' ')
    {
        biggest_num_m = biggest_num_m * t + (tmp[i] - '0');
        ++i;
    }
    memset(tmp,0,sizeof(tmp));
    fgets(tmp,NUM_4,fp);
    fgets(tmp,NUM_4,fp);
    memset(tmp,0,sizeof(tmp));

    int j = 0,k,_tmp;
    int tmp_N = N,tmp_M = M;
    //matrix->B = (int **)malloc(sizeof(int) * )
    while(tmp_N--)//存储论文中的B数组
    {
        i = 0;
        k = 0;
        fgets(tmp,NUM_3,fp);
        while(tmp[i] != '\n')
        {
			_tmp = 0;
            while(tmp[i] != ' ')
            {
                _tmp = _tmp * t + (tmp[i] - '0');
                ++i;
            }
            B[j][k] = _tmp-1;
            H[j][_tmp-1] = 1;
			++k;
            ++i;
        }
        ++j;
    }
    j=0;
    while(tmp_M--)//存储论文中的A数组
    {
        i = 0;
        k = 0;
        fgets(tmp,NUM_2,fp);
        while(tmp[i] != '\n' && tmp[i] != '\0')
        {
			_tmp = 0;
            while(tmp[i] != ' ')
            {
                _tmp = _tmp * t + (tmp[i] - '0');
                ++i;
            }
            A[j][k++] = _tmp-1;
            H[_tmp-1][j] = 1;
            ++i;
        }
        ++j;
    }
    fclose(fp);
}

double Gauss()
{
	double ret;
	double UA,UB;
	static double U1,U2;
	double s;
	static double fac;
	static int phase = 0;

	if(phase == 0)
	{
		do
		{
			UA = (float)rand()/(RAND_MAX);
			UB = (float)rand()/(RAND_MAX);
			U1 = 1-2*UA;
			U2 = 1-2*UB;
			s = U1*U1 + U2*U2;
		}while(s>=1 || s<=0);  //s<=0 for avoiding log(s) to be -inf.
		fac = sqrt(-2.0*SIGMA*SIGMA*log(s)/s);
		ret = U1*fac;
	}
	else
	{
		ret = U2*fac;
	}
	phase = 1-phase;
	return ret;
}
