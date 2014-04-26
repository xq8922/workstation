#include <stdio.h>
#include <math.h>
#define PI 3.1415
#define direc "E:\\workstation\\workstation\\sumProductDecoding\\H.txt"
#define I 30

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
double AWGN(double *,double *,double *);
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
     while(fp != EOF)
     {
         fprintf(fp,"%d %d\n",a->N,a->M);

     }
}

double AWGN(double x,double everage,double sigama)
{
    return pow(M_E,-1*(x - everage)*(x - everage)/(2*sigama*sigama))/(sqrt(2*PI)*sigama);
}
