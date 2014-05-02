#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.1415
#define direc "E:\\workstation\\workstation\\sumProductDecoding\\H.txt"
#define I 30
#define NUM_1 12
#define NUM_2 20
#define NUM_3 200
#define NUM_4 330
#define CODE_LEN 1008

const int m[CODE_LEN] = {0} ;
int m_change[CODE_LEN] ;

typedef struct {
  int N , M ;      /* size of the matrix */
  int **mlist;     /* list of integer coordinates in the m direction where the non-zero entries are */
  int **nlist;     /* list of integer coordinates in the n direction where the non-zero entries are */
  int *num_mlist;  /* weight of each row, m */
  int *num_nlist;  /* weight of each column n */
  int *l_up_to ;
  int *u_up_to ;
  int *norder ;
  int biggest_num_m ;       /* actual biggest sizes */
  int biggest_num_n ;
  int biggest_num_m_alloc ; /* sizes used for memory allocation */
  int biggest_num_n_alloc ;
  int tot ;
  int same_length ;  /* whether all vectors in mlist and nlist have same length */
} alist_matrix ;

void write_alist(FILE *,alist_matrix *);
enum P{p1 = 5,p2,p3,p4,p5};
double Gauss();
double e(int *B,double *);

int main()
{

    return 0;
}

void write_alist ( FILE *fp , alist_matrix *a )
 {
  /* this assumes that mlist and nlist have the form of a rectangular
     matrix in the file; if lists have unequal lengths, then the
     entries should be present (eg zero values) but are ignored
     */
    if(!(fp=fopen(direc,"r")))
       return;
    char tmp[NUM_4];
    fgets(tmp,NUM_1,fp);
    int i=0,t = 1;
    a->M = a->N = 0;
    while(tmp[i] != ' ')
    {
        (a->N) = (a->N) * t + (tmp[i] - '0');
        t *= 10;
        ++i;
    }
    ++i;
    t=1;
    while(tmp[i] != ' ')
    {
        (a->M) = (a->M) * t + (tmp[i] - '0');
        t *= 10;
        ++i;
    }
    memset(tmp,0,sizeof(tmp));
    fgets(tmp,NUM_1,fp);
    i = 0;
    t = 1;
    while(tmp[i] != ' ')
    {
        (a->biggest_num_n) = (a->biggest_num_n) * t + (tmp[i] - '0');
        t *= 10;
        ++i;
    }
    ++i;
    t = 1;
    while(tmp[i] != ' ')
    {
        (a->biggest_num_m) = (a->biggest_num_m) * t + (tmp[i] - '0');
        t *= 10;
        ++i;
    }
    memset(tmp,0,sizeof(tmp));
    fgets(tmp,NUM_4,fp);
    fgets(tmp,NUM_4,fp);
    memset(tmp,0,sizeof(tmp));

    int j = 0,k;
    int tmp_N = a->N,tmp_M = a->M;
    while(tmp_N--)//存储论文中的B数组
    {
        i = 0;
        t = 1;
        _tmp = 0;
        k = 0;
        fgets(tmp,NUM_1,fp);
        while(tmp[i] != '\n')
        {
            while(tmp[i] != ' ')
            {
                _tmp = _tmp * t + (tmp[i] - '0');
                t *= 10;
                ++i;
            }
            a->nlist[j][k++] = _tmp;
            ++i;
        }
        ++j;
    }
    j=0;
    while(tmp_M--)//存储论文中的A数组
    {
        i = 0;
        t = 1;
        _tmp = 0;
        k = 0;
        fgets(tmp,NUM_2,fp);
        while(tmp[i] != '\n')
        {
            while(tmp[i] != ' ')
            {
                _tmp = _tmp * t + (tmp[i] - '0');
                t *= 10;
                ++i;
            }
            a->mlist[j][k++] = _tmp;
            ++i;
        }
        ++j;
    }

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
