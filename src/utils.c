
/* 
 * File: utils.cpp
 * Author: Miquel Orobitg, Fernando Cores, Fernando Guirado
 *
 */

#include "utils.h"


double alloc_mem;
double max_mem;
double tot_mem;
Memcontrol *memlast;

/*
 * Open a file
 */

FILE *openfile(const char *filename, const char *mode){

    FILE *fp;
    
    switch (mode[0]){
        case 'r':
            if((fp = fopen(filename, mode)) == NULL){
                fprintf(stderr, "\nERROR - Opening the file %s\n", filename);
                exit(-1);
            }
            break;
        case 'w':
            if((fp = fopen(filename, mode)) == NULL){
                fprintf(stderr, "\nERROR - Opening the file %s\n", filename);
                exit(-1);
            }
            break;
        default:
            fprintf(stderr, "\nERROR - There is some error in the call of openfile. Unknown mode\n");
            exit(-1);
    }

    return fp;
}

Fname *declare_fname (int size){

    Fname *F;

    size+=1000+FILENAMELEN+1;

    F=vcalloc ( 1, sizeof (Fname));
    F->name = (char *) vcalloc(size, sizeof (char));
    F->path = (char *) vcalloc(size, sizeof (char));
    F->suffix= (char *) vcalloc(size, sizeof (char));
    F->full= (char *) vcalloc(size, sizeof (char));
    
    return F;
}

Fname *free_fname(Fname *F){
    vfree(F->name);
    vfree(F->path);
    vfree(F->suffix);
    vfree(F->full);
    return NULL;
}

Paths *declare_paths(){
    
    Paths *PB;
    
    PB=vcalloc ( 1, sizeof (Paths));
    
    PB->mta_home= PB->tcoffee_bin = (char *) vcalloc(ALLPATH, sizeof (char));
    PB->tcoffee_bin = (char *) vcalloc(ALLPATH, sizeof (char));
    PB->clustalw_bin = (char *) vcalloc(ALLPATH, sizeof (char));
    PB->clustalo_bin = (char *) vcalloc(ALLPATH, sizeof (char));
    PB->mafft_bin = (char *) vcalloc(ALLPATH, sizeof (char));
    PB->normd_bin = (char *) vcalloc(ALLPATH, sizeof (char));
    PB->strike_bin = (char *) vcalloc(ALLPATH, sizeof (char));
    PB->tc_score_bin = (char *) vcalloc(ALLPATH, sizeof (char));
    
    PB->nwtomafft = (char *) vcalloc(ALLPATH, sizeof (char));
    
    return PB;    
}

Paths *free_paths(Paths *PB){
    //vfree(PB->tcoffee_bin);
    //vfree(PB->clustalw_bin);
    //vfree(PB->clustalo_bin);
    //vfree(PB->normd_bin);
    //vfree(PB->strike_bin);
    //vfree(PB->tc_score_bin);
    
    return NULL;
}

/*
 * Parse the file names to the name structures

*/

Fname *parse_fname(char* array){

    int l;
    Fname *F;
    F = declare_fname(sizeof(array));
    sprintf(F->full, "%s", array);
    sprintf(F->path, "%s", array);
    l=strlen (array);
    while (l!=-1 && (F->path)[l]!='/')(F->path)[l--]='\0';

    sprintf(F->name, "%s", array+l+1);
    l=strlen(F->name);
    while(l!=-1){
        if((F->name)[l]=='.'){
            F->name[l]='\0';
            sprintf (F->suffix, "%s", F->name+l+1);
            break;
        }
        else{ l--;
        }
    }
    
    return F;
}

Paths *parse_paths(){
    
    Paths *PB;
    
    PB=declare_paths();
    
    PB->mta_home = getenv("MTA_HOME");
    printf("%s\n", PB->mta_home);
    PB->tc_score_bin=getenv("TCOFFEE_BIN");
    if(PB->tc_score_bin==NULL)
        fprintf(stderr, "Warning: Triplet o Coffee path not found. Scores disabled\n");
    //sprintf(PB->tc_score_bin, "%s/bin/plugins/t_coffee", PB->tc_score_bin);
    
    PB->tcoffee_bin=getenv("TCOFFEE_BIN");
    if(PB->tcoffee_bin==NULL)
        fprintf(stderr, "Warning: T-Coffee path not found. T-coffee method disabled\n");
    
    PB->clustalw_bin=getenv("CLUSTALW_BIN");
    if(PB->clustalw_bin==NULL)
        fprintf(stderr, "Warning: ClustalW path not found. ClustalW method disabled\n");
    
    PB->clustalo_bin=getenv("CLUSTALO_BIN");
    if(PB->clustalo_bin==NULL)
        fprintf(stderr, "Warning: ClustalO path not found. ClustalO method disabled\n");
    
    PB->mafft_bin=getenv("MAFFT_BIN");
    if(PB->mafft_bin==NULL)
        fprintf(stderr, "Warning: Mafft path not found. Mafft method disabled\n");
    
    PB->normd_bin=getenv("NORMD_BIN");
    if(PB->normd_bin==NULL)
        fprintf(stderr, "Warning: NorMD path not found. NorMD score disabled\n");
    
    PB->strike_bin=getenv("STRIKE_BIN");
    if(PB->strike_bin==NULL)
        fprintf(stderr, "Warning: STRIKE path not found. STRIKE score disabled\n");
    
    sprintf(PB->nwtomafft, "%s/src/scripts/newick2mafft.rb", PB->mta_home);
    printf("%s\n", PB->mta_home);
    
    //printf("%s\n%s\n%s\n%s\n%s\n%s\n", PB->tcoffee_bin,  PB->clustalw_bin,  PB->clustalo_bin, PB->normd_bin, PB->strike_bin, PB->tc_score_bin);
    
    return PB;
}

void remove_file(char* file){
    
    FILE *fp;
    
    if((fp=fopen(file, "r")) != NULL){
        if(remove(file)==-1){
            fprintf(stderr, "\nWARNING: File %s not deleted.", file);
        }
        fclose(fp);
    }
}

FILE *print_array_char (FILE *out, char **array, int n, char *sep){
    int a;
    if ( array==NULL || read_size_char (array,sizeof (char*))<n){
        fprintf( stderr, "\nORB in print_array_char [FATAL]\n");
        exit(EXIT_FAILURE);
        crash("");
    }
    for(a=0; a< n; a++){
        fprintf(out, "%s%s", array[a],sep);
    }
    return out;
}

 /* Skip the spaces of a file*/
FILE *skip_space(FILE *fp){
    int   c;

    do {
     	c = getc(fp);
    } while(isspace(c));
    if ( c==EOF){
        fprintf ( stderr, "\nEOF");
        exit(EXIT_FAILURE);
    }
    ungetc(c, fp);
    return fp;
}

/*
 * Count the number of characters of a file

int count_n_char_in_file(const char *name){

    int c, n;
    FILE *fp;

    n = 0;
    openfile(fp, name, "r");

    while ((c = fgetc(fp)) != EOF){
        n++;
    }
    fclose(fp);

    return n;
}


/*
 * Obtain the time
 */


double get_time(){
    struct timeval t;
    gettimeofday(&t, 0);
    return t.tv_sec + (0.000001 * t.tv_usec);
}

int name_is_in_list(char *name, char **name_list, int n_name, int len){
    int a;
    int pos=-1;
	/*Note: RETURNS THE Offset of the LAST Occurence of name in name_list*/

    if ( name_list==NULL || name ==NULL)return -1;

    for ( a=0; a< n_name; a++){
        if ( name_list[a]==NULL);
        else if ( len!=-1){
            if(strncmp(name, name_list[a], len)==0){
                pos=a;
            }
        }
        else if(strm( name, name_list[a])){
            pos=a;
        }
    }
    return pos;
}


/*
 * Memory utils
 */

int vstrcmp(const char *s1, const char *s2) {
    if(!s1 && !s2) return 0;
    else if (!s1 || !s2) return 1;
    else return strcmp(s1, s2);
}

int vstrncmp(const char *s1, const char *s2, int n)
{
    if(!s1 && !s2) return 0;
    else if(!s1 || !s2) return 1;
    else return strncmp(s1, s2, n);
}

#define READ_ARRAY_SIZE(type, function)\
int function (void *array, size_t size)\
    {\
      return read_array_size (array, size);\
    }
READ_ARRAY_SIZE(short,read_size_short)
READ_ARRAY_SIZE(char,read_size_char)
READ_ARRAY_SIZE(int,read_size_int)
READ_ARRAY_SIZE(float,read_size_float)
READ_ARRAY_SIZE(double,read_size_double)

int read_array_size(void *array, size_t size) {
    Memcontrol *p;
    if(array==NULL) return 0;
    p=(Memcontrol *)array;
    p-=2;
    if(p[0].size_element ==0 && size==0){
        fprintf(stderr, "ERROR - In read_array_size: trying to read the size of a malloced block\n");
    }
    else if( size ==0){
        return (int)p[0].size/p[0].size_element;
    }

    return (int)p[0].size/size;
}

int read_array_size_new(void *array) {
    return read_array_size( array, 0);
}

int is_dynamic_memory(void *array) {
    Memcontrol *p;
    if(array==NULL) return 0;
    p=(Memcontrol *)array;
    p-=2;
    if(strm(p[0].check, "dy")) {
        return 1;
    }
    return 0;
}

int verify_memory(int s) {

    alloc_mem+=s;
    tot_mem=(alloc_mem>tot_mem)?alloc_mem:tot_mem;

    if(max_mem && alloc_mem>max_mem) {
        fprintf(stderr, "ERROR - Requires Too Much Memory: %d Megabytes\n", (int)(alloc_mem/1024*1024));
        fprintf(stderr, "ERROR - Tip: Rerun your Job with a smaller dataset\n");
        exit(-1);
    }
    else {
        return 1;
    }
    return 0;
}

void vfree(void *p) {
   Memcontrol *M;
   size_t size;

   if(!p) return;
   else {
       M= (Memcontrol*) p;
       M-=2;
       size=M[0].size;

       p=M;
       free(p);

       verify_memory(-(size+2*sizeof(Memcontrol)));
    }
}

void *vmalloc(size_t size) {
    void *x;
    Memcontrol *M;

    verify_memory(size+2*sizeof(Memcontrol));

    if(size==0) {
        return NULL; /*crash ("\n0 bytes in vmalloc\n");*/
    }
    else {
        x= malloc(size + 2*sizeof(Memcontrol));
        if(x==NULL){
            fprintf(stderr, "ERROR - FAILED TO ALLOCATE REQUIRED MEMORY (vmalloc)\n");
            exit(-1);
        }
        else {
            M = (Memcontrol *) x;
            M[0].size=size;
            M[0].size_element=0;
            sprintf(M[0].check, "dy");
            M+=2;
            x=M;
            return x;
        }
    }
    return NULL;
}

void *sub_vcalloc(size_t nobj, size_t size, int MODE) {

    void *x;
    Memcontrol *M;

    if(nobj<=0 || size<=0) return NULL;/*crash ("\n0 bytes in vmalloc\n");*/
    else x=vmalloc(nobj*size);

    M = (Memcontrol *) x;
    M-=2;
    M[0].size_element=size;
    M+=2;
    x=M;
    if(x==NULL){
        fprintf(stderr, "\nFAILED TO ALLOCATE REQUIRED MEMORY (vcalloc)\n");
        return NULL;
    }
    else {
        if(MODE==MEMSET0) {
            x=memset(x,0, nobj*size);
        }
        else {
            if(nobj) x=memset (x, 0, size);
        }
        return x;
    }
}


void *vcalloc(size_t nobj, size_t size) {
    return sub_vcalloc(nobj,size, MEMSET0);
}

void *vcalloc_nomemset ( size_t nobj, size_t size){
  return sub_vcalloc (nobj, size, NO_MEMSET0);
}

void *vrealloc ( void *p, size_t size){

    void *x;
    Memcontrol *M;
    size_t i_size;
    int a;
    
    if(p==NULL) {
        x=vmalloc(size);
        memset(x, 0, size);
        return x;
    }
    else {
        M = (Memcontrol *) p;
        M-=2;
        i_size=M[0].size;
        p=M;
        if(size<=0){
            return NULL;
            vfree(p);
            return NULL;
        }
        else {
            verify_memory(size - i_size);
            x=realloc(p, size+2*sizeof(Memcontrol));
            if(x==NULL){
                fprintf(stderr, "\nFAILED TO ALLOCATE REQUIRED MEMORY (realloc)\n");
                return NULL;
            }
            M= (Memcontrol *) x;
            M[0].size=size;
            M+=2;
            x=M;
            for(a=i_size; a< size; a++) ((char*)x)[a]=0;
            return x;
        }
    }
    return NULL;
}

void *ckalloc(size_t bytes){
    register void *ret;
    extern void *vcalloc (size_t nelem, size_t elsize);

    if( (ret = vcalloc(bytes, sizeof(char))) == NULL){
        fprintf(stderr, "ERROR - Out of memory\n");
    }
    else {
        return ret;
    }
    return ret;
}


void *ckvrealloc(void *ptr, size_t bytes){
    register void *ret;
    extern void *vrealloc (void *ptr, size_t size);

    if( (ret = vrealloc(ptr, bytes)) == NULL){
            fprintf(stderr, "ERROR - Out of memory\n");
    }
    else{
        return ret;
    }
    return ret;
}


void ckfree(void *ptr){
    vfree(ptr);
}


/************************************************************************/
/*                                                                      */
/*             DECLARE 2d ARRAYS                                        */
/*                                                                      */
/*                                                                      */
/************************************************************************/

void * free_arrayN(void *p, int n)
{
  int a, s;
  void **i;


  if ( p==NULL) return NULL;
  else if ( n==1)vfree ((void *)p);
  else
    {
      i=(void**)p;
      s=read_array_size ( (void *)p, sizeof ( void *));
      for ( a=0; a< s; a++)free_arrayN ((void *)i[a], n-1);
      vfree (p);
    }
  return NULL;
}

void * declare_arrayNnomemset (int ndim, size_t size, ...)
{
  va_list ap;
  int *array;
  void **p;
  int a;

  va_start (ap, size);

  array=vcalloc (ndim, sizeof (int));
  for ( a=0; a< ndim; a++)
    {
      array[a]=va_arg (ap,int);
      if ( array[a]<0){va_end(ap);return NULL;}

    }
  va_end (ap);

  if ( ndim==2)
    {

      p=vcalloc_nomemset (array[0], sizeof ( void*));
      for (a=0; a< array[0]; a++)
	{
	p[a]=vcalloc_nomemset (array[1], size);
	}
    }
  else
    {
      p=declare_arrayN2nomemset (ndim, array, size);
    }
  vfree (array);
  return p;
}

void *declare_arrayN2nomemset ( int ndim, int *A, size_t size)
{
  int a;
  void **p;

  if ( ndim>1)
    {
      p=vcalloc_nomemset (A[0], sizeof (void*));
      for ( a=0; a<A[0]; a++)
    	p[a]=declare_arrayN2(ndim-1, A+1, size);
    }
  else
    {
      p=vcalloc_nomemset (A[0], size);
    }
  return p;
}

void * declare_arrayN (int ndim, size_t size, ...)
{
  va_list ap;
  int *array;
  void **p;
  int a;


  va_start (ap, size);

  array=vcalloc (ndim, sizeof (int));
  for ( a=0; a< ndim; a++)
    {
      array[a]=va_arg (ap,int);
      if ( array[a]<0){va_end(ap);return NULL;}

    }
  va_end (ap);

  if ( ndim==2)
    {

      p=vcalloc_nomemset (array[0], sizeof ( void*));
      for (a=0; a< array[0]; a++)
	p[a]=vcalloc (array[1], size);
    }
  else
    {
      p=declare_arrayN2 (ndim, array, size);
    }
  vfree (array);
  return p;
}

void *declare_arrayN2 ( int ndim, int *A, size_t size)
{
  int a;
  void **p;

  if ( ndim>1)
    {
      p=vcalloc_nomemset (A[0], sizeof (void*));
      for ( a=0; a<A[0]; a++)
    	p[a]=declare_arrayN2(ndim-1, A+1, size);
    }
  else
    {
      p=vcalloc (A[0], size);
    }
  return p;
}

void **declare_array (int first, int second, size_t size)
{
  return (void **)declare_arrayN (2,size,first, second);
}


#define DECLARE_ARRAY(type,wf,rf,function)\
type**  function (int first, int second)\
  {\
    return (type **)declare_arrayN (2,sizeof(type), first, second);\
   }

int **declare_int2 (int f, int *s, int d)
{
  int **r;
  int a;
  r=vcalloc ( f, sizeof (int*));
  for (a=0; a<f; a++)
    r[a]=vcalloc (s[a]+d, sizeof (int));
  return r;
}


DECLARE_ARRAY(short,write_size_short,read_size_short,declare_short)
DECLARE_ARRAY(char,write_size_char,read_size_char,declare_char)
DECLARE_ARRAY(int,write_size_int,read_size_int,declare_int)
DECLARE_ARRAY(float,write_size_float,read_size_float,declare_float)
DECLARE_ARRAY(double,write_size_double,read_size_double,declare_double)

void **declare_array_nomemset (int first, int second, size_t size)
{
  return (void **)declare_arrayNnomemset (2,size,first, second);
}
#define DECLARE_ARRAY_NMS(type,wf,rf,function)\
type**  function (int first, int second)\
  {\
    return (type **)declare_arrayNnomemset (2,sizeof(type), first, second);\
   }
DECLARE_ARRAY_NMS(short,write_size_short,read_size_short,declare_short_nomemset)
DECLARE_ARRAY_NMS(char,write_size_char,read_size_char,declare_char_nomemset)
DECLARE_ARRAY_NMS(int,write_size_int,read_size_int,declare_int_nomemset)
DECLARE_ARRAY_NMS(float,write_size_float,read_size_float,declare_float_nomemset)
DECLARE_ARRAY_NMS(double,write_size_double,read_size_double,declare_double_nomemset)


/************************************************************************/
/*                                                                      */
/*             Realloc 2d ARRAYS                                        */
/*                                                                      */
/*                                                                      */
/************************************************************************/
void ** realloc_arrayN(int ndim,void **main_array,size_t size, ...)
{
  va_list ap;
  int *array;
  void **p;
  int a;

  /*dim size==-1: keep current size*/
  /*Dim sizes are the absolute size (not the extension*/
  /*If array getting shorter, memory is Not Claimed back*/

  array=vcalloc (ndim, sizeof (int));
  va_start (ap, size);
  for ( a=0; a< ndim; a++)
    {
      array[a]=va_arg (ap,int);
      if ( array[a]<-1){va_end(ap);return NULL;}
    }
  va_end (ap);

  p=realloc_arrayN2 (ndim, main_array, array,size);
  vfree (array);
  return p;
}

void **realloc_arrayN2 ( int ndim, void ** p, int *A, size_t size)
{
  int a;
  int o, n;


  if (ndim==0)return NULL;

  if ( ndim>1)
    {
      o=read_array_size (p,sizeof (void*));
      if (A[0]>o)p=vrealloc (p, sizeof (void*)*A[0]);
      n=(A[0]==-1)?o:A[0];
      for ( a=0; a<n; a++)
    	p[a]=realloc_arrayN2(ndim-1,p[a], A+1, size);
    }
  else
    {
      o=read_array_size (p, size);
      if (A[0]>o)p=vrealloc (p, size*A[0]);
    }
  return p;
}

void ** realloc_array (void **array,size_t size, int first, int second, int ext1, int ext2)
{
  int a;
  int d1, d2;
  if ( array==NULL)return declare_array (((first==-1)?0:first)+ext1, ((second==-1)?0:second)+ext2, size);
  else if ( first==-1)
    {
      first=read_array_size (array, sizeof (void *));
    }
  if (second==-1)second=read_array_size(array[0], size);

  d1=first+ext1;
  d2=second+ext2;

  for ( a=d1; a<first; a++)vfree (array[a]);
  array=vrealloc (array, sizeof (void*)*d1);

  if ( d2!=second)
    {
      for (a=0; a<d1 && a<first; a++)
	array[a]=vrealloc ( array[a], size*d2);
    }
  for ( a=first; a< d1; a++)
    array[a]=vrealloc ( array[a], size*d2);
  return array;
}

#define REALLOC_ARRAY(type,wf,rf,function1,function2,function3)\
type ** function1 ( type **array, int first, int second, int ext1, int ext2)\
    {\
return (type **)realloc_array ((void **)array, sizeof (type),first, second, ext1, ext2);\
\
    }
REALLOC_ARRAY(short,write_size_short,read_size_short,realloc_short,declare_short,free_short)
REALLOC_ARRAY(char,write_size_char,read_size_char,realloc_char,declare_char,free_char)
REALLOC_ARRAY(int,write_size_int,read_size_int,realloc_int,declare_int,free_int)
REALLOC_ARRAY(float,write_size_float,read_size_float,realloc_float,declare_float,free_float)
REALLOC_ARRAY(double,write_size_double,read_size_double,realloc_double,declare_double,free_double)

#define NEW_REALLOC_ARRAY(type,wf,rf,function1,function2,function3)\
type ** function1 ( type **array, int ext1, int ext2)\
    {int a, b;\
     int first, l1;\
     int second, l2;\
     type **new_array;\
\
     first=rf(array,sizeof (type*));\
     second=rf(array[0],sizeof (type));\
     \
     if ( ext1==-1)ext1=first;\
     if ( ext2==-1)ext2=second;\
     l1=MIN(ext1, first);\
     l2=MIN(ext2, second);\
     new_array=declare_arrayN(2,sizeof(type),ext1, ext2);\
     for ( a=0; a<l1; a++)\
       for ( b=0; b<l2; b++)new_array[a][b]=array[a][b];\
     function3(array, -1);\
     return new_array;\
    }

NEW_REALLOC_ARRAY(short,write_size_short,read_size_short,new_realloc_short,declare_short,free_short)
NEW_REALLOC_ARRAY(char,write_size_char,read_size_char,new_realloc_char,declare_char,free_char)
NEW_REALLOC_ARRAY(int,write_size_int,read_size_int,new_realloc_int,declare_int,free_int)
NEW_REALLOC_ARRAY(float,write_size_float,read_size_float,new_realloc_float,declare_float,free_float)
NEW_REALLOC_ARRAY(double,write_size_double,read_size_double,new_realloc_double,declare_double,free_double)

/************************************************************************/
/*                                                                      */
/*            free 2d ARRAYS                                        */
/*                                                                      */
/*                                                                      */
/************************************************************************/
#define FREE_ARRAY(type,wf,rf,function) \
type ** function (type **array, int first)\
    {\
      return free_arrayN((void*)array, 2);\
    }
FREE_ARRAY(short,write_size_short,read_size_short,free_short)
FREE_ARRAY(char,write_size_char,read_size_char,free_char)
FREE_ARRAY(int,write_size_int,read_size_int,free_int)
FREE_ARRAY(float,write_size_float,read_size_float,free_float)
FREE_ARRAY(double,write_size_double,read_size_double,free_double)


int **duplicate_int(int **array, int len, int field){
    return copy_int (array, declare_int(len, field), len, field);
}

int **copy_int(int **array1, int **array2, int len, int number_field){

    int a;

    if (array1==NULL) return NULL;
    if (len==-1)len=read_size_int (array1, sizeof (int*));
    if (number_field==-1)number_field=read_size_int (array1[0],sizeof (int));


    if (array2)free_int (array2, -1);
    array2=declare_int ( len, number_field);

    for (a=0; a< len; a++){
	ga_memcpy_int(array1[a],array2[a],number_field);
    }

    return array2;
}

int *ga_memcpy_int(int *array1, int *array2, int n){
    int a;

    for ( a=0; a< n; a++){
        array2[a]=array1[a];
    }

    return array2;
}

int **sim_array2dist_array(int **p, int max){

    int s1, s2, a, b;

    s1=read_array_size ((void *)p, sizeof (void *));
    s2=read_array_size ((void*)p[0],sizeof (int));

    /*s2=read_array_size ((void*)p[0],sizeof (void *)); OLD before 64 Bits stuff*/
    for (a=0; a< s1; a++){
        for(b=0; b< s2; b++){
            p[a][b]=max-(int)p[a][b];
        }
    }

    return p;
}

int **normalize_array(int **p, int max, int norm){

    int s1, s2, a, b;

    s1=read_array_size ((void *)p, sizeof (void *));
    s2=read_array_size ((void*)p[0],sizeof (int));

    /*s2=read_array_size ((void*)p[0],sizeof (void *)); OLD before 64 Bits stuff*/
    for(a=0; a< s1; a++){
        for(b=0; b< s2; b++){
            p[a][b]=(p[a][b]*norm)/max;
        }
    }

    return p;
}

#ifdef MPI_FLAG

//Not USED
void write_lib(char *libname){
    
    FILE *fp;
    char *lib;
    char c;
    long int size=0, i=0;
    
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
   // printf("sizeo: %ld\n", size); 
    lib = vcalloc(size, sizeof(char));
    MPI_Bcast(lib, size, MPI_CHAR, 0, MPI_COMM_WORLD);
    //printf("%s\n", lib);
    
    fp = fopen(libname, "w");
/*
    for(i=0; i<size; i++){
        fprintf(fp, "%c", lib[i]);
        //printf("%c", lib[i]);
    }
*/
    fprintf(fp, "%s", lib);
    //printf("%s\n", lib);
    fclose(fp);
/*
    fp = fopen(libname, "r");
    while(c != EOF){
        c = fgetc(fp);
        printf("%d", c);
        //lib[i] = c;
        i++;
    }
    
    fclose(fp);
    
*/ 
    vfree(lib);   
}

void read_lib(char *libname){
    
    FILE *fp;
    char *lib;
    char c=' ';
    long int size=0, i=0;
    
    fp = fopen(libname, "r");
    
    fseek(fp,0,SEEK_END); //Nos vamos el final del archivo
    size = (ftell(fp) * sizeof(char));
    //printf("Size: %ld\n", size);
    lib = vcalloc(size, sizeof(char));
    fclose(fp);
    fp = fopen(libname, "r");  
   
    while(c != EOF){
        c = fgetc(fp);
        //printf("c: %d", c);
        lib[i] = c;
        i++;
    }

    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(lib, size, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    fclose(fp);
    
    vfree(lib);
    
}

#endif