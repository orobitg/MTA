
#include "guide_tree_utils.h"

static int g_distance_tree, g_rooted_tree;
static int g_tot_nseq;
static int first_seq, last_seq;
static int nseqs;
static char **names;       /* the seq. names */
static double **tmat;
static double *av;
static double *left_branch, *right_branch;
static int *tkill;
int rooted_tree;
int distance_tree;

FILE * create_tree(NT_node ptree, NT_node parent,int *nseq,int  *ntotal,int  *nnodes,NT_node **lu, FILE *fp){

  	int i, type, ret;
  	float dist=0;
	float bootstrap=0;
  	char *name;
  	int ch;


	name=vcalloc ( MAXNAMES+1, sizeof (char));
	sprintf ( name, "   ");
  	fp=skip_space(fp);
  	ch = (char)getc(fp);

  	if (ch == '(')
    		{
    		type = NODE;
      		name[0] = '\0';
      		lu[0][ntotal[0]] = lu[1][nnodes[0]] = ptree;
      		ntotal[0]++;
      		nnodes[0]++;
      		create_tree_node(ptree, parent);
      		fp=create_tree(ptree->left, ptree, nseq,ntotal,nnodes,lu,fp);
     		ch = (char)getc(fp);
                if ( ch == ',')
       			{
          		fp=create_tree(ptree->right, ptree,nseq,ntotal,nnodes,lu,fp);
			ch = (char)getc(fp);
          		if ( ch == ',')
            			{

               			ptree = insert_tree_node(ptree);
               			lu[0][ntotal[0]] = lu[1][nnodes[0]] = ptree;
               			ntotal[0]++;
               			nnodes[0]++;
               			fp=create_tree(ptree->right, ptree,nseq,ntotal,nnodes,lu,fp);
               			rooted_tree = FALSE;

            			}
       			}

      		fp=skip_space(fp);
      		ch = (char)getc(fp);
    		}
 	else
    		{
     	 	type=LEAF;
      		lu[0][ntotal[0]] = lu[2][nseq[0]] = ptree;
      		ntotal[0]++;
      		nseq[0]++;
      		name[0] = ch;
      		i=1;
      		ch = (char)getc(fp);
		if ( name[0]=='\'')
		  {
		    /*This protects names that are between single quotes*/
		    while ( ch!='\'')
		      {
			if (i < MAXNAMES) name[i++] = ch;
          		ch = (char)getc(fp);
		      }
		    if (i < MAXNAMES) name[i++] = ch;
		    while ((ch != ':') && (ch != ',') && (ch != ')'))ch = (char)getc(fp);
		  }
		else
		  {
		    while ((ch != ':') && (ch != ',') && (ch != ')'))
		      {
			if (i < MAXNAMES) name[i++] = ch;
			ch = (char)getc(fp);
		      }
		  }

      		name[i] = '\0';

      		if ( i>=(MAXNAMES+1)){fprintf (stderr, "\nName is too long"); exit(EXIT_FAILURE);}
      		if (ch != ':' && !isdigit(ch))
         		{
			  /*distance_tree = FALSE*/;
         		}
		}
	if (ch == ':')
     		{
       		fp=skip_space(fp);
       		ret=fscanf(fp,"%f",&dist);
       		fp=skip_space(fp);
		bootstrap=0;
       		}
	/*Tree with Bootstrap information*/
	else if (isdigit (ch))
	  {
	    ungetc(ch,fp);
	    ret=fscanf(fp,"%f",&bootstrap);
	    if(fscanf(fp,":%f",&dist)==1);
	    else dist=0;
	    fp=skip_space(fp);
	  }
	else
	  {
	    ungetc ( ch, fp);
	    skip_space(fp);
	  }

	set_info(ptree, parent, type, name, dist, bootstrap);


  	vfree (name);
  	return fp;
}

/*
 * Declare, create and fill information of a tree node.
 */

NT_node declare_tree_node (int nseq){

    NT_node p;

    p= (NT_node)vcalloc (1, sizeof ( Treenode));
    p->left = NULL;
    p->right = NULL;
    p->parent = NULL;
    p->dist = 0.0;
    p->leaf = 0;
    p->order = 0;
    p->maxnseq=nseq;
    p->name=(char*)vcalloc (MAXNAMES+1,sizeof (char));
    p->name[0]='\0';
    p->lseq=(int*)vcalloc ( nseq, sizeof (int));
    p->aligned = 0;
    p->iscluster = 0;
    p->isorcluster = 0;
    p->ant = -1;
    p->cl_id = -1;
    p->nnodes = -1;
    return p;
}

void set_info(NT_node p, NT_node parent, int pleaf, char *pname, float pdist, float bootstrap){
    p->parent = parent;
    p->leaf = pleaf;
    p->dist = pdist;
    p->bootstrap=bootstrap;
    p->order = 0;

    sprintf (p->name, "%s", pname);

    if (pleaf ==1){
        p->left = NULL;
        p->right = NULL;
    }
}

NT_node insert_tree_node(NT_node pptr){

    NT_node newnode;

    newnode = declare_tree_node( pptr->maxnseq);
    create_tree_node(newnode, pptr->parent);

    newnode->left = pptr;
    pptr->parent = newnode;

    set_info(newnode, pptr->parent, 0, "", 0.0, 0);

    return(newnode);
}

void create_tree_node(NT_node pptr, NT_node parent){

    pptr->parent = parent;
    pptr->left =declare_tree_node(pptr->maxnseq) ;
    (pptr->left)->parent=pptr;

    pptr->right =declare_tree_node(pptr->maxnseq) ;
    (pptr->right)->parent=pptr;
}


NT_node reroot(NT_node ptree, int nseq, int ntotal, int nnodes, NT_node **lu){

	NT_node p, rootnode, rootptr;
 	float   diff, mindiff=0, mindepth = 1.0, maxdist;
	int   i;
	int first = TRUE;



	rootptr = ptree;

   	for (i=0; i<ntotal; i++)
     		{
        	p = lu[0][i];
        	if (p->parent == NULL)
           		diff = calc_root_mean(p, &maxdist, nseq, lu);
        	else
           		diff = calc_mean(p, &maxdist, nseq, lu);

        	if ((diff == 0) || ((diff > 0) && (diff < 2 * p->dist)))
         		 {
              		if ((maxdist < mindepth) || (first == TRUE))
                 		{
                    		first = FALSE;
                    		rootptr = p;
                    		mindepth = maxdist;
                    		mindiff = diff;
                 		}
           		}

     		}
	if (rootptr == ptree)
     		{
        	mindiff = rootptr->left->dist + rootptr->right->dist;
        	rootptr = rootptr->right;
     		}

   	rootnode = insert_root(rootptr, mindiff);
	diff = calc_root_mean(rootnode, &maxdist, nseq, lu);

	return(rootnode);
}

float calc_root_mean(NT_node root, float *maxdist, int nseq, NT_node **lu){

	float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
	NT_node p;
	int i;
	int nl, nr;
	int direction;


   	dist = (*maxdist) = 0;
   	nl = nr = 0;
   	for (i=0; i< nseq; i++)
   	  {
          p = lu[2][i];
          dist = 0.0;
          while (p->parent != root)
         	{
               	dist += p->dist;
               	p = p->parent;
           	}
         if (p == root->left) direction = LEFT;
         else direction = RIGHT;
         dist += p->dist;

         if (direction == LEFT)
           	{
           	lsum += dist;
             	nl++;
           	}
         else
           	{
             	rsum += dist;
             	nr++;
           	}

        if (dist > (*maxdist)) *maxdist = dist;
     	}

      lmean = lsum / nl;
      rmean = rsum / nr;

      diff = lmean - rmean;
      return(diff);
}

float calc_mean(NT_node nptr, float *maxdist, int nseq,NT_node **lu){
    
   	float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
   	NT_node p, *path2root;
   	float *dist2node;
   	int depth = 0, i,j , n;
   	int nl , nr;
   	int direction, found;


	path2root = (NT_node *)vcalloc(nseq,sizeof(Treenode));
	dist2node = (float *)vcalloc(nseq,sizeof(float));

   	depth = (*maxdist) = dist = 0;
   	nl = nr = 0;
   	p = nptr;
   	while (p != NULL)
     		{
         	path2root[depth] = p;
         	dist += p->dist;
         	dist2node[depth] = dist;
         	p = p->parent;
         	depth++;
     		}

/***************************************************************************
   *nl = *nr = 0;
   for each leaf, determine whether the leaf is left or right of the node.
   (RIGHT = descendant, LEFT = not descendant)
****************************************************************************/
   	for (i=0; i< nseq; i++)
     		{
       		p = lu[2][i];
       		if (p == nptr)
         		{
            		direction = RIGHT;
            		dist = 0.0;
         		}
       		else
         		{
         		direction = LEFT;
         		dist = 0.0;

        		found = FALSE;
         		n = 0;
         		while ((found == FALSE) && (p->parent != NULL))
           			{
               			for (j=0; j< depth; j++)
                 			if (p->parent == path2root[j])
                    				{
                      				found = TRUE;
                      				n = j;
                    				}
               			dist += p->dist;
               			p = p->parent;
           			}

         		if (p == nptr) direction = RIGHT;

			}
         	if (direction == LEFT)
           		{
             		lsum += dist;
             		lsum += dist2node[n-1];
             		nl++;
           		}
         	else
           		{
             		rsum += dist;
             		nr++;
           		}

        	if (dist > (*maxdist)) *maxdist = dist;
     		}

   vfree(dist2node);
   vfree(path2root);



   if ( nl==0 || nr==0)
     {
       exit(EXIT_FAILURE);
     }
   lmean = lsum / nl;
   rmean = rsum / nr;

   diff = lmean - rmean;
   return(diff);
}

NT_node insert_root(NT_node p, float diff){

   NT_node newp, prev, q, t;
   float dist, prevdist,td;

   newp = declare_tree_node( p->maxnseq);
   t = p->parent;

   prevdist = t->dist;
   p->parent = newp;

   dist = p->dist;

   p->dist = diff / 2;
   if (p->dist < 0.0) p->dist = 0.0;
   if (p->dist > dist) p->dist = dist;

   t->dist = dist - p->dist;

   newp->left = t;
   newp->right = p;
   newp->parent = NULL;
   newp->dist = 0.0;
   newp->leaf = NODE;

   if (t->left == p) t->left = t->parent;
   else t->right = t->parent;

   prev = t;
   q = t->parent;

   t->parent = newp;

   while (q != NULL)
     {
        if (q->left == prev)
           {
              q->left = q->parent;
              q->parent = prev;
              td = q->dist;
              q->dist = prevdist;
              prevdist = td;
              prev = q;
              q = q->left;
           }
        else
           {
              q->right = q->parent;
              q->parent = prev;
              td = q->dist;
              q->dist = prevdist;
              prevdist = td;
              prev = q;
              q = q->right;
           }
    }
/*
   remove the old root node
*/
   q = prev;
   if (q->left == NULL)
      {
         dist = q->dist;
         q = q->right;
         q->dist += dist;
         q->parent = prev->parent;
         if (prev->parent->left == prev)
            prev->parent->left = q;
         else
            prev->parent->right = q;
         prev->right = NULL;
      }
   else
      {
         dist = q->dist;
         q = q->left;
         q->dist += dist;
         q->parent = prev->parent;
         if (prev->parent->left == prev)
            prev->parent->left = q;
         else
            prev->parent->right = q;
         prev->left = NULL;
      }

   return(newp);
}


int tree_file2nseq (char *fname){
    FILE *fp;
    char *string;
    int p, a, b, c, n;

    string=vcalloc (count_n_char_in_file(fname)+1, sizeof (char));

    fp=openfile(fname, "r");
    n=0;
    while((c=fgetc(fp))!=EOF){
        if (c=='(' || c==')' || c==',' || c==';'){
            string[n++]=c;
        }
    }
    fclose (fp);string[n]='\0';

    for (n=0, p=1; p<strlen (string); p++){
        a=string[p-1];
        b=string[p];

        if ( a=='(' && b==',')n++;
        else if ( a==',' && b==')')n++;
        else if ( a==',' && b==',')n++;
    }
    //if ( getenv4debug("DEBUG_TREE"))fprintf (stderr, "\n[DEBUG_TREE:tree_file2nseq]%s", string);
    vfree (string);
    return n;
}

/*Produce a guide tree*/

NT_node **make_tree(Distance_matrix *DM, Sequence *S, char *tree_file, char *mode, int tree, int nrand, char *method){

    NT_node **T=NULL;
    int **distances;
    int out_nseq;
    char **out_seq_name;
    char **out_seq;
    int *out_seq_id;
  
    if (S->nseq==2){
        fprintf(stderr, "Not enought sequences\n");
    }
    else if(strm(mode, "nj") || strm(mode, "nj_bgt")){
        out_nseq=S->nseq;
        out_seq_name=S->seq_names;
        out_seq=S->seq;
        out_seq_id=S->seq_id;
        distances=duplicate_int(DM->score_similarity_matrix, -1, -1);

        distances=sim_array2dist_array (distances, MAXID*SCORE_K);
        distances=normalize_array (distances, MAXID*SCORE_K, 100);

        T=make_nj_tree(distances, out_seq,out_seq_name, out_seq_id, out_nseq, tree_file, mode, tree, nrand, method);
    }
    else {
        fprintf(stderr, "Bad mode\n");
    }
    free_int(distances, out_nseq);

    return T;
}

/*NJ clustering mode functions*/

NT_node **make_nj_tree(int **distances, char **out_seq, char **out_seq_name, int *out_seq_id, int out_nseq, char *tree_file, char *mode, int tree, int nrand, char *method){

    if ( distances==NULL){
        fprintf ( stderr, "\nError: make_nj_tree, must provide an alignment or a distance matrix.");
        exit(EXIT_FAILURE);
    }

    return int_dist2nj_tree(distances,out_seq_name, out_seq_id, out_nseq, tree_file, mode, tree, nrand, method);
}

NT_node **int_dist2nj_tree(int **distances, char **out_seq_name, int *out_seq_id, int out_nseq, char *tree_file, char *mode, int tree, int nrand, char *method){

    int a, b;
    double **d;
    NT_node **T;

    d=declare_double(out_nseq, out_nseq);
    for ( a=0; a< out_nseq; a++){
        for ( b=0; b< out_nseq; b++){
            d[a][b]=distances[a][b];
        }
    }

    T = dist2nj_tree(d, out_seq_name, out_seq_id, out_nseq, tree_file, mode, tree, nrand, method);

    free_double (d, -1);
    return T;
}

NT_node **float_dist2nj_tree(float **distances, char **out_seq_name, int *out_seq_id, int out_nseq, char *tree_file, char *mode, int tree, int nrand, char *method){

    int a, b;
    double **d;
    NT_node **T;

    d=declare_double( out_nseq, out_nseq);
    for ( a=0; a< out_nseq; a++){
        for ( b=0; b< out_nseq; b++){
            d[a][b]=distances[a][b];
        }
    }

    T=dist2nj_tree (d, out_seq_name, out_seq_id, out_nseq, tree_file, mode, tree, nrand, method);

    free_double (d, -1);

    return T;
}

NT_node **dist2nj_tree (double **distances, char **out_seq_name, int *out_seq_id, int out_nseq, char *tree_file, char *mode, int tree, int nrand, char *method){

    int a, b;
    double **d_dist;
    int tot_node=0;
    NT_node **T;   

    if(!tree_file){
        fprintf(stderr, "No tree file\n");
    }
    d_dist=declare_double(out_nseq+1, out_nseq+1);

    for(a=0; a<out_nseq; a++){
        for(b=0; b< out_nseq; b++){
            if (a!=b){
                d_dist[a+1][b+1]=distances[a][b]/MAXID;
            }
            else{
              d_dist[a+1][b+1]=0;
            }
        }
    }
   
    guide_tree(tree_file, d_dist, out_seq_name, out_nseq, mode, tree, nrand, method);
    free_double(d_dist,-1);
    T=read_tree(tree_file,&tot_node,out_nseq, out_seq_name, out_seq_id, method);

    return T;
}

void guide_tree(char *fname, double **saga_tmat, char **saga_seq_name, int saga_nseq, char *mode, int tree, int nrand, char *method){

    int i;
    FILE *fp;
    char **standard_tree;

    nseqs=saga_nseq;
    first_seq=1;
    last_seq=nseqs;

    names=saga_seq_name;
    tmat=saga_tmat;

    standard_tree   = (char **) vmalloc( (nseqs+1) * sizeof (char *) );
    for(i=0; i<nseqs+1; i++){
        standard_tree[i]  = (char *) vmalloc( (nseqs+1) * sizeof(char));
    }

    nj_tree(standard_tree, nseqs, mode, tree, nrand);

    fp=openfile(fname, "w");

    print_phylip_tree(standard_tree, fp, FALSE);


    vfree(left_branch);
    vfree(right_branch);
    vfree(tkill);
    vfree(av);

    for (i=0;i<nseqs+1;i++){
        vfree(standard_tree[i]);
    }

    fclose(fp);
 
}

void nj_tree(char **tree_description, int nseq, char *mode, int tree, int nrand){  
    
    if(strm(mode, "nj")){   
        if (nseq<100){
            slow_nj_tree_mta(tree_description, tree, nrand);
        }
        else{
            fast_nj_tree_mta(tree_description, tree, nrand);
        }
    }
    else if(strm(mode, "nj_bgt")){
        if (nseq<100){
            slow_nj_tree_bgt(tree_description); //arreglar      
        }
        else{
            fast_nj_tree_bgt(tree_description); //arreglar
        }
    }
    else {
        fprintf(stderr, "Bad mode\n");
    }
}

void print_phylip_tree(char **tree_description, FILE *tree, int bootstrap){

    fprintf(tree,"(\n");

    two_way_split(tree_description, tree, nseqs-2,1,bootstrap);
    fprintf(tree,":%7.5f,\n",left_branch[nseqs-2]);
    two_way_split(tree_description, tree, nseqs-2,2,bootstrap);
    fprintf(tree,":%7.5f,\n",left_branch[nseqs-1]);
    two_way_split(tree_description, tree, nseqs-2,3,bootstrap);

    fprintf(tree,":%7.5f)",left_branch[nseqs]);
    if (bootstrap) fprintf(tree,"TRICHOTOMY");
    fprintf(tree,";\n");
}

/*
void print_phylip_tree_clo(char **tree_description, FILE *tree, int bootstrap){

    fprintf(tree,"(\n(\n");

    two_way_split(tree_description, tree, nseqs-2,1,bootstrap);
    fprintf(tree,":%7.5f,\n",left_branch[nseqs-2]);
    two_way_split(tree_description, tree, nseqs-2,2,bootstrap);
    fprintf(tree,":%7.5f)\n",left_branch[nseqs-1]);
    fprintf(tree,":%7.5f,\n",left_branch[nseqs-1]);
    two_way_split(tree_description, tree, nseqs-2,3,bootstrap);

    fprintf(tree,":%7.5f)",left_branch[nseqs]);
    if (bootstrap) fprintf(tree,"TRICHOTOMY");
    fprintf(tree,";\n");
}
*/

void two_way_split(char **tree_description, FILE *tree, int start_row, int flag, int bootstrap){

    int row, new_row, col, test_col=0;
    int single_seq;

    if(start_row != nseqs-2){
        fprintf(tree,"(\n");
    }

    for(col=1; col<=nseqs; col++){
        if(tree_description[start_row][col] == flag){
            test_col = col;
            break;
        }
    }

    single_seq = TRUE;
    for(row=start_row-1; row>=1; row--){
        if(tree_description[row][test_col] == 1){
            single_seq = FALSE;
            new_row = row;
            break;
        }
    }

    if(single_seq) {
        tree_description[start_row][test_col] = 0;
        fprintf(tree,"%s",names[test_col+0-1]);
    }
    else{
        for(col=1; col<=nseqs; col++) {
            if((tree_description[start_row][col]==1) && (tree_description[new_row][col]==1)){
                tree_description[start_row][col] = 0;
            }
        }
        two_way_split(tree_description, tree, new_row, (int)1, bootstrap);
    }

    if(start_row == nseqs-2) {
        /*		if (bootstrap && (flag==1)) fprintf(tree,"[TRICHOTOMY]");
        */
        return;
    }

    fprintf(tree,":%7.5f,\n",left_branch[start_row]);

    for(col=1; col<=nseqs; col++){
        if(tree_description[start_row][col] == flag) {
            test_col = col;
            break;
        }
    }

    single_seq = TRUE;
    for(row=start_row-1; row>=1; row--){
        if(tree_description[row][test_col] == 1) {
            single_seq = FALSE;
            new_row = row;
            break;
        }
    }

    if(single_seq){
        tree_description[start_row][test_col] = 0;
        fprintf(tree,"%s",names[test_col+0-1]);
    }
    else {
        for(col=1; col<=nseqs; col++){
            if((tree_description[start_row][col]==1) && (tree_description[new_row][col]==1)){
                tree_description[start_row][col] = 0;
            }
        }
        two_way_split(tree_description, tree, new_row, (int)1, bootstrap);
    }

    fprintf(tree,":%7.5f)\n",right_branch[start_row]);

}

/*Read a gruide tree from a file*/

NT_node** read_tree(char *treefile, int *tot_node,int nseq, char **seq_names, int *out_seq_id, char *method){

	/*The Tree Root is in the TREE[3][0]...*/
	/*TREE[0][ntot]--> pointer to each node and leave*/
    char ch;
    int  a,b;

    FILE *fp;
    int nseq_read = 0;
    int nnodes = 0;/*Number of Internal Nodes*/
    int ntotal = 0;/*Number of Internal Nodes + Number of Leaves*/
    int flag;
    int c_seq;
    NT_node **lu_ptr;
    NT_node seq_tree, root,p;

    rooted_tree=distance_tree=TRUE;

    fp = openfile(treefile, "r");
    fp=skip_space(fp);
    ch = (char)getc(fp);
    if(ch != '('){
        fprintf(stderr, "ERROR: Wrong format in tree file %s\n", treefile);
        exit(EXIT_FAILURE);
    }
    rewind(fp);

    lu_ptr=(NT_node **)vcalloc(4,sizeof(NT_node*));
    lu_ptr[0] = (NT_node *)vcalloc(10*nseq,sizeof(NT_node));
    lu_ptr[1] = (NT_node *)vcalloc(10*nseq,sizeof(NT_node));
    lu_ptr[2] = (NT_node *)vcalloc(10*nseq,sizeof(NT_node));
    lu_ptr[3] =(NT_node *) vcalloc(1,sizeof(NT_node));

    seq_tree =(NT_node) declare_tree_node(nseq);

    set_info(seq_tree, NULL, 0, "  ", 0.0, 0);
  
    fp=create_tree(seq_tree,NULL,&nseq_read, &ntotal, &nnodes, lu_ptr, fp);
    fclose (fp);
    
    if(distance_tree == FALSE){
        if(rooted_tree == FALSE){
            fprintf(stderr,"ERROR: input tree is unrooted and has no distances, cannot align sequences\n");
            exit(EXIT_FAILURE);
        }
    }

    if(rooted_tree == FALSE){
        root = reroot(seq_tree, nseq,ntotal,nnodes, lu_ptr);
        lu_ptr[1][nnodes++]=lu_ptr[0][ntotal++]=root;
    }
    else {
        root = seq_tree;
    }
    
    lu_ptr[3][0]=root;
    root->nnodes = ntotal;
    tot_node[0]=nnodes;

    for(a=0; a< ntotal; a++){
        (lu_ptr[0][a])->isseq=(lu_ptr[0][a])->leaf;
        (lu_ptr[0][a])->dp=(lu_ptr[0][a])->dist;
    }

    for(a=0; a< nseq; a++){
        if (!seq_names){
            flag=1;
            (lu_ptr[2][a])->order=(lu_ptr[2][a])->seq=a;
        }
        else {
            for(flag=0,b=0; b<nseq; b++){
                if(strncmp ( (lu_ptr[2][a])->name, seq_names[b], MAXNAMES)==0){
                    flag=1;
                    (lu_ptr[2][a])->order=(lu_ptr[2][a])->seq=b;
                        /*vfree ( (lu_ptr[2][a])->name);*/
                    sprintf ((lu_ptr[2][a])->name, "%s", seq_names[b]);
                    lu_ptr[2][a]->seq_id = out_seq_id[b];
                    //fprintf(stdout, "%d - %s\n", lu_ptr[2][a]->seq_id, lu_ptr[2][a]->name);
                }
            }
        }
        /*
        if ( flag==0  && (lu_ptr[0][a])->leaf==1)
              {
                fprintf ( stderr, "\n%s* not in tree",(lu_ptr[2][a])->name);
                for ( a=0; a< ntotal; a++)
                  {
                    fprintf ( stderr, "\n%d %s",(lu_ptr[2][a])->leaf, (lu_ptr[2][a])->name);
                  }
              }
        */
    }

    if (seq_names){
        int tnseq;
        char *s;
        char **tree_names;
        int fail_flag=0;
        tnseq=tree_file2nseq(treefile);
        tree_names=vcalloc ( tnseq, sizeof (char*));
        for (a=0; a<tnseq; a++){
            s=(lu_ptr[2][a])->name;
            tree_names[a]=s;
            if (name_is_in_list(s, seq_names, nseq, MAXNAMES+1)==-1){
                fprintf (stderr, "\nERROR: Sequence %s in the tree [%s] is not in the alignment[FATAL]\n", s, treefile);
                fail_flag=1;
            }
        }
        for (a=0; a<nseq; a++){
            s=seq_names[a];
            if (name_is_in_list(s, tree_names, nseq, MAXNAMES+1)==-1){
                fprintf (stderr, "\nERROR: Sequence %s in the sequences is not in the tree [%s][FATAL]\n", s, treefile);
                fail_flag=1;
            }
        }
        vfree (tree_names);
        if(fail_flag==1){
            exit(EXIT_FAILURE);
        }
    }

    for ( a=0; a< nseq; a++){
        p=lu_ptr[2][a];
        c_seq=p->seq;

        while ( p!=NULL){
            p->lseq[p->nseq]=c_seq;
            p->nseq++;
            p=p->parent;
        }
    }

    return lu_ptr;
}

void slow_nj_tree_mta(char **tree_description, int tree, int nrand){

    int i;
    int l[4],nude,k;
    int nc,mini,minj,j,ii,jj;
    double fnseqs,fnseqs2=0,sumd;
    double diq,djq,dij,d2r,dr,dio,djo,da;
    double tmin,total,dmin;
    double bi,bj,b1,b2,b3,branch[4];
    int typei,typej;             /* 0 = node; 1 = OTU */

    double average=0.0, rand_num=0.0;
    int cont=0, nseq_rand=0;

    fnseqs = (double)nseqs;
    nseq_rand = (nseqs * nrand)/100;

/*********************** First initialisation ***************************/

    mini = minj = 0;

    left_branch = (double *) vcalloc( (nseqs+2),sizeof (double)   );
    right_branch = (double *) vcalloc( (nseqs+2),sizeof (double)   );
    tkill = (int *)    vcalloc( (nseqs+1),sizeof (int) );
    av = (double *) vcalloc( (nseqs+1),sizeof (double)   );

    for(i=1;i<=nseqs;++i){
        tmat[i][i] = av[i] = 0.0;
        tkill[i] = 0;
    }

/*********************** Enter The Main Cycle ***************************/
       	/**start main cycle**/
    for(nc=1; nc<=(nseqs-3); ++nc){
        cont=0;
        sumd = 0.0;
        for(j=2; j<=nseqs; ++j){
            for(i=1; i<j; ++i){
                tmat[j][i] = tmat[i][j];
                sumd = sumd + tmat[i][j];
                if(nc==1){
                    cont++;
                }
            }
        }

        if(nc==1){
            average=sumd/cont;
            //intf("\nAVERAGE: %f\n", average);
        }
        tmin = 99999.0;
        

/*.................compute SMATij values and find the smallest one ........*/       
        for(jj=2; jj<=nseqs; ++jj){
            if(tkill[jj] != 1){
                for(ii=1; ii<jj; ++ii){
                    if(tkill[ii] != 1){
                        diq = djq = 0.0;

                        for(i=1; i<=nseqs; ++i){
			  diq = diq + tmat[i][ii];
			  djq = djq + tmat[i][jj];
			}
                        if(nc <= nseq_rand /*&& tree != 0*/){
                            rand_num = ((rand()/((double) RAND_MAX + 1)) * INTERVAL);
                            //printf("\tRand Num: %lf -", rand_num);
                            diq = diq + rand_num;
                            rand_num = ((rand()/((double) RAND_MAX + 1)) * INTERVAL);
                            //printf(" %lf\n", rand_num);
                            dij = djq + rand_num;
                        }
      
                        //printf("\tRAND: %f\n", rand_num);
                        dij = tmat[ii][jj];
                        d2r = diq + djq - (2.0*dij);
                        dr  = sumd - dij -d2r;
                        fnseqs2 = fnseqs - 2.0;
                        total= d2r+ fnseqs2*dij +dr*2.0;
                        total= total / (2.0*fnseqs2);
                      //total = total + rand_num;
		      /* commented out to have consistent results with CYGWIN: if(total < tmin)"*/
                        if(total < tmin)
                        {
                            if ( tmin-total>MY_EPS){
                                tmin = total;
                            }
                            mini = ii;
                            minj = jj;
                        }
                    }
                }
            }
        }
/*.................compute branch lengths and print the results ........*/
            //printf("mini: %d minj: %d - tmin: %f\n", mini, minj, tmin);

        dio = djo = 0.0;
        for(i=1; i<=nseqs; ++i) {
            dio = dio + tmat[i][mini];
            djo = djo + tmat[i][minj];
        }

        dmin = tmat[mini][minj];
        dio = (dio - dmin) / fnseqs2;
        djo = (djo - dmin) / fnseqs2;
        bi = (dmin + dio - djo) * 0.5;
        bj = dmin - bi;
        bi = bi - av[mini];
        bj = bj - av[minj];

        if(av[mini] > 0.0){
            typei = 0;
        }
        else{
            typei = 1;
        }
        if(av[minj] > 0.0){
            typej = 0;
        }
        else{
            typej = 1;
        }
/*
   set negative branch lengths to zero.  Also set any tiny positive
   branch lengths to zero.
*/
        if( fabs(bi) < 0.0001) bi = 0.0;
        if( fabs(bj) < 0.0001) bj = 0.0;


        left_branch[nc] = bi;
        right_branch[nc] = bj;

        for(i=1; i<=nseqs; i++){
            tree_description[nc][i] = 0;
        }

        if(typei == 0) {
            for(i=nc-1; i>=1; i--)
                if(tree_description[i][mini] == 1){
                    for(j=1; j<=nseqs; j++)
                        if(tree_description[i][j] == 1)
                            tree_description[nc][j] = 1;
                        break;
		}
        }
        else{
            tree_description[nc][mini] = 1;
        }

        if(typej == 0) {
            for(i=nc-1; i>=1; i--)
                if(tree_description[i][minj] == 1) {
                    for(j=1; j<=nseqs; j++)
                        if(tree_description[i][j] == 1)
                            tree_description[nc][j] = 1;
                        break;
            }
        }
        else{
            tree_description[nc][minj] = 1;
        }


/*
   Here is where the -0.00005 branch lengths come from for 3 or more
   identical seqs.
*/
/*		if(dmin <= 0.0) dmin = 0.0001; */
        if(dmin <= 0.0) dmin = 0.000001;
        av[mini] = dmin * 0.5;

 /*........................Re-initialisation................................*/

        fnseqs = fnseqs - 1.0;
        tkill[minj] = 1;

        for(j=1; j<=nseqs; ++j){
            if(tkill[j] != 1 ){
                da = ( tmat[mini][j] + tmat[minj][j] ) * 0.5;
            if((mini - j) < 0 )
                tmat[mini][j] = da;
            if((mini - j) > 0)
                tmat[j][mini] = da;
            }
        }

        for(j=1; j<=nseqs; ++j){
            tmat[minj][j] = tmat[j][minj] = 0.0;
        }

    }
	/*end main cycle**/

/******************************Last Cycle (3 Seqs. left)********************/

    nude = 1;

    for(i=1; i<=nseqs; ++i){
        if(tkill[i] != 1 ){
            l[nude] = i;
            nude = nude + 1;
        }
    }

    b1 = (tmat[l[1]][l[2]] + tmat[l[1]][l[3]] - tmat[l[2]][l[3]]) * 0.5;
    b2 =  tmat[l[1]][l[2]] - b1;
    b3 =  tmat[l[1]][l[3]] - b1;

    branch[1] = b1 - av[l[1]];
    branch[2] = b2 - av[l[2]];
    branch[3] = b3 - av[l[3]];

/* Reset tiny negative and positive branch lengths to zero */
    if( fabs(branch[1]) < 0.0001) branch[1] = 0.0;
    if( fabs(branch[2]) < 0.0001) branch[2] = 0.0;
    if( fabs(branch[3]) < 0.0001) branch[3] = 0.0;

    left_branch[nseqs-2] = branch[1];
    left_branch[nseqs-1] = branch[2];
    left_branch[nseqs]   = branch[3];

    for(i=1; i<=nseqs; i++){
        tree_description[nseqs-2][i] = 0;
    }



    for(i=1; i<=3; ++i) {
        if( av[l[i]] > 0.0) {
            for(k=nseqs-3; k>=1; k--)
                if(tree_description[k][l[i]] == 1) {
                    for(j=1; j<=nseqs; j++)
                        if(tree_description[k][j] == 1)
                            tree_description[nseqs-2][j] = i;
                        break;
                }
        }
        else{
            tree_description[nseqs-2][l[i]] = i;
        }
        if(i < 3){
        }
    }
}

void fast_nj_tree_mta(char **tree_description, int tree, int nrand){

    register int i;
    int l[4],nude,k;
    int nc,mini,minj,j,ii,jj;
    double fnseqs,fnseqs2=0,sumd;
    double diq,djq,dij,dio,djo,da;
    double tmin,total,dmin;
    double bi,bj,b1,b2,b3,branch[4];
    int typei,typej;             /* 0 = node; 1 = OTU */

    /* IMPROVEMENT 1, STEP 0 : declare  variables */
    double *sum_cols, *sum_rows, *join;

    /* IMPROVEMENT 2, STEP 0 : declare  variables */
    int loop_limit;
    typedef struct _ValidNodeID {
        int n;
        struct _ValidNodeID *prev;
        struct _ValidNodeID *next;
    } ValidNodeID;
    ValidNodeID *tvalid, *lpi, *lpj, *lpii, *lpjj, *lp_prev, *lp_next;

    /*
     * correspondence of the loop counter variables.
     *   i .. lpi->n,	ii .. lpii->n
     *   j .. lpj->n,	jj .. lpjj->n
     */

    int cont=0, nseq_rand=0;
    double average=0.0, rand_num=0.0;

    fnseqs = (double)last_seq-first_seq+1;

    nseq_rand = ((last_seq-first_seq+1) * nrand)/100;

/*********************** First initialisation ***************************/


    if (fnseqs == 2) {
        return;
    }

    mini = minj = 0;

    left_branch 	= (double *) ckalloc( (nseqs+2) * sizeof (double)   );
    right_branch    = (double *) ckalloc( (nseqs+2) * sizeof (double)   );
    tkill 		= (int *) ckalloc( (nseqs+1) * sizeof (int) );
    av   		= (double *) ckalloc( (nseqs+1) * sizeof (double)   );

    /* IMPROVEMENT 1, STEP 1 : Allocate memory */
    sum_cols	= (double *) ckalloc( (nseqs+1) * sizeof (double)   );
    sum_rows	= (double *) ckalloc( (nseqs+1) * sizeof (double)   );
    join		= (double *) ckalloc( (nseqs+1) * sizeof (double)   );

    /* IMPROVEMENT 2, STEP 1 : Allocate memory */
    tvalid	= (ValidNodeID *) ckalloc( (nseqs+1) * sizeof (ValidNodeID) );
    /* tvalid[0] is special entry in array. it points a header of valid entry list */
    tvalid[0].n = 0;
    tvalid[0].prev = NULL;
    tvalid[0].next = &tvalid[1];

    /* IMPROVEMENT 2, STEP 2 : Construct and initialize the entry chain list */
    for(i=1, loop_limit = last_seq-first_seq+1,
    lpi=&tvalid[1], lp_prev=&tvalid[0], lp_next=&tvalid[2] ;
    i<=loop_limit ;
    ++i, ++lpi, ++lp_prev, ++lp_next){
        tmat[i][i] = av[i] = 0.0;
        tkill[i] = 0;
        lpi->n = i;
        lpi->prev = lp_prev;
        lpi->next = lp_next;

        /* IMPROVEMENT 1, STEP 2 : Initialize arrays */
        sum_cols[i] = sum_rows[i] = join[i] = 0.0;
    }
    tvalid[loop_limit].next = NULL;

    /*
     * IMPROVEMENT 1, STEP 3 : Calculate the sum of score value that
     * is sequence[i] to others.
     */
    sumd = 0.0;
    for(lpj=tvalid[0].next ; lpj!=NULL ; lpj = lpj->next){
        double tmp_sum = 0.0;
        j = lpj->n;
        /* calculate sum_rows[j] */
        for (lpi=tvalid[0].next ; lpi->n < j ; lpi = lpi->next) {
            i = lpi->n;
            tmp_sum += tmat[i][j];
            /* tmat[j][i] = tmat[i][j]; */
        }
        sum_rows[j] = tmp_sum;

        tmp_sum = 0.0;
        /* Set lpi to that lpi->n is greater than j */
        if ((lpi != NULL) && (lpi->n == j)) {
            lpi = lpi->next;
        }
        /* calculate sum_cols[j] */
        for( ; lpi!=NULL ; lpi = lpi->next) {
            i = lpi->n;
            tmp_sum += tmat[j][i];
            cont++;
            /* tmat[i][j] = tmat[j][i]; */
        }
        sum_cols[j] = tmp_sum;
    }

    //srand(time(NULL));

/*********************** Enter The Main Cycle ***************************/

    for(nc=1, loop_limit = (last_seq-first_seq+1-3); nc<=loop_limit; ++nc){

        sumd = 0.0;
        /* IMPROVEMENT 1, STEP 4 : use sum value */
        for(lpj=tvalid[0].next ; lpj!=NULL ; lpj = lpj->next) {
            sumd += sum_cols[lpj->n];
        }
        if(nc==1){
            average=sumd/cont;
        }

        /* IMPROVEMENT 3, STEP 0 : multiply tmin and 2*fnseqs2 */
        fnseqs2 = fnseqs - 2.0;		/* Set fnseqs2 at this point. */
        tmin = 99999.0 * 2.0 * fnseqs2;

/*.................compute SMATij values and find the smallest one ........*/

        mini = minj = 0;

        /* jj must starts at least 2 */
        if ((tvalid[0].next != NULL) && (tvalid[0].next->n == 1)){
            lpjj = tvalid[0].next->next;
        } else {
            lpjj = tvalid[0].next;
        }

        for( ; lpjj != NULL; lpjj = lpjj->next){
            jj = lpjj->n;
            for(lpii=tvalid[0].next ; lpii->n < jj ; lpii = lpii->next){
                ii = lpii->n;
                diq = djq = 0.0;
                diq = sum_cols[ii] + sum_rows[ii];
                djq = sum_cols[jj] + sum_rows[jj];
                if(nc <= nseq_rand/* && tree != 0*/){
                /* IMPROVEMENT 1, STEP 4 : use sum value */

                    rand_num = (rand()/((double) RAND_MAX + 1)) * INTERVAL;
                    diq = diq +  rand_num;
                    rand_num = (rand()/((double) RAND_MAX + 1)) * INTERVAL;
                    djq = djq +  rand_num;
                }
                /*
                 * always ii < jj in this point. Use upper
                 * triangle of score matrix.
                 */
                dij = tmat[ii][jj];

                /*
                 * IMPROVEMENT 3, STEP 1 : fnseqs2 is
                 * already calculated.
                 */
                /* fnseqs2 = fnseqs - 2.0 */

                /* IMPROVEMENT 4 : transform the equation */
  /*-------------------------------------------------------------------*
   * OPTIMIZE of expression 'total = d2r + fnseqs2*dij + dr*2.0'       *
   * total = d2r + fnseq2*dij + 2.0*dr                                 *
   *       = d2r + fnseq2*dij + 2(sumd - dij - d2r)                    *
   *       = d2r + fnseq2*dij + 2*sumd - 2*dij - 2*d2r                 *
   *       =       fnseq2*dij + 2*sumd - 2*dij - 2*d2r + d2r           *
   *       = fnseq2*dij + 2*sumd - 2*dij - d2r                         *
   *       = fnseq2*dij + 2*sumd - 2*dij - (diq + djq - 2*dij)         *
   *       = fnseq2*dij + 2*sumd - 2*dij - diq - djq + 2*dij           *
   *       = fnseq2*dij + 2*sumd - 2*dij + 2*dij - diq - djq           *
   *       = fnseq2*dij + 2*sumd  - diq - djq                          *
   *-------------------------------------------------------------------*/
                total = fnseqs2*dij + 2.0*sumd  - diq - djq;

                /*
                 * IMPROVEMENT 3, STEP 2 : abbrevlate
                 * the division on comparison between
                 * total and tmin.
                 */
                /* total = total / (2.0*fnseqs2); */

                if(total < tmin) {
                    tmin = total;
                    mini = ii;
                    minj = jj;
                }
            }
        }

		/* MEMO: always ii < jj in avobe loop, so mini < minj */

/*.................compute branch lengths and print the results ........*/
                //printf("mini: %d, minj: %d, tmin: %f\n", mini, minj, tmin);

        dio = djo = 0.0;

        /* IMPROVEMENT 1, STEP 4 : use sum value */
        dio = sum_cols[mini] + sum_rows[mini];
        djo = sum_cols[minj] + sum_rows[minj];

        dmin = tmat[mini][minj];
        dio = (dio - dmin) / fnseqs2;
        djo = (djo - dmin) / fnseqs2;
        bi = (dmin + dio - djo) * 0.5;
        bj = dmin - bi;
        bi = bi - av[mini];
        bj = bj - av[minj];

        if(av[mini] > 0.0 ){
                typei = 0;
        }
        else{
            typei = 1;
        }
        if(av[minj] > 0.0 ){
            typej = 0;
        }
        else{
            typej = 1;
        }

/*
   set negative branch lengths to zero.  Also set any tiny positive
   branch lengths to zero.
*/
        if( fabs(bi) < 0.0001) bi = 0.0;
        if( fabs(bj) < 0.0001) bj = 0.0;


        left_branch[nc] = bi;
        right_branch[nc] = bj;

        for(i=1; i<=last_seq-first_seq+1; i++){
            tree_description[nc][i] = 0;
        }

        if(typei == 0){
            for(i=nc-1; i>=1; i--)
                if(tree_description[i][mini] == 1) {
                    for(j=1; j<=last_seq-first_seq+1; j++)
                         if(tree_description[i][j] == 1)
                            tree_description[nc][j] = 1;
                        break;
                    }
        }
        else{
            tree_description[nc][mini] = 1;
        }

        if(typej == 0){
            for(i=nc-1; i>=1; i--)
                if(tree_description[i][minj] == 1) {
                    for(j=1; j<=last_seq-first_seq+1; j++)
                        if(tree_description[i][j] == 1)
                            tree_description[nc][j] = 1;
                        break;
                    }
        }
        else{
            tree_description[nc][minj] = 1;
        }

/*
   Here is where the -0.00005 branch lengths come from for 3 or more
   identical seqs.
*/
/*		if(dmin <= 0.0) dmin = 0.0001; */
        if(dmin <= 0.0) dmin = 0.000001;
        av[mini] = dmin * 0.5;

/*........................Re-initialisation................................*/

        fnseqs = fnseqs - 1.0;
        tkill[minj] = 1;

        /* IMPROVEMENT 2, STEP 3 : Remove tvalid[minj] from chain list. */
        /* [ Before ]
         *  +---------+        +---------+        +---------+
         *  |prev     |<-------|prev     |<-------|prev     |<---
         *  |    n    |        | n(=minj)|        |    n    |
         *  |     next|------->|     next|------->|     next|----
         *  +---------+        +---------+        +---------+
         *
         * [ After ]
         *  +---------+                           +---------+
         *  |prev     |<--------------------------|prev     |<---
         *  |    n    |                           |    n    |
         *  |     next|-------------------------->|     next|----
         *  +---------+                           +---------+
         *                     +---------+
         *              NULL---|prev     |
         *                     | n(=minj)|
         *                     |     next|---NULL
         *                     +---------+
         */
        (tvalid[minj].prev)->next = tvalid[minj].next;
        if (tvalid[minj].next != NULL) {
                (tvalid[minj].next)->prev = tvalid[minj].prev;
        }
        tvalid[minj].prev = tvalid[minj].next = NULL;

		/* IMPROVEMENT 1, STEP 5 : re-calculate sum values. */
        for(lpj=tvalid[0].next ; lpj != NULL ; lpj = lpj->next) {
            double tmp_di = 0.0;
            double tmp_dj = 0.0;
            j = lpj->n;

            /*
             * subtrace a score value related with 'minj' from
             * sum arrays .
             */
            if(j < minj){
                tmp_dj = tmat[j][minj];
                sum_cols[j] -= tmp_dj;
            }
            else if(j > minj){
                tmp_dj = tmat[minj][j];
                sum_rows[j] -= tmp_dj;
            } /* nothing to do when j is equal to minj. */


            /*
             * subtrace a score value related with 'mini' from
             * sum arrays .
             */
            if(j < mini){
                tmp_di = tmat[j][mini];
                sum_cols[j] -= tmp_di;
            }
            else if (j > mini){
                tmp_di = tmat[mini][j];
                sum_rows[j] -= tmp_di;
            } /* nothing to do when j is equal to mini. */

                /*
                 * calculate a score value of the new inner node.
                 * then, store it temporary to join[] array.
                 */
            join[j] = (tmp_dj + tmp_di) * 0.5;
        }

        /*
         * 1)
         * Set the score values (stored in join[]) into the matrix,
         * row/column position is 'mini'.
         * 2)
         * Add a score value of the new inner node to sum arrays.
         */
        for(lpj=tvalid[0].next ; lpj != NULL; lpj = lpj->next) {
            j = lpj->n;
            if(j < mini){
                tmat[j][mini] = join[j];
                sum_cols[j] += join[j];
            }
            else if(j > mini){
                tmat[mini][j] = join[j];
                sum_rows[j] += join[j];
            } /* nothing to do when j is equal to mini. */
        }

        /* Re-calculate sum_rows[mini],sum_cols[mini]. */
        sum_cols[mini] = sum_rows[mini] = 0.0;

        /* calculate sum_rows[mini] */
        da = 0.0;
        for(lpj=tvalid[0].next ; lpj->n < mini ; lpj = lpj->next) {
            da += join[lpj->n];
        }
        sum_rows[mini] = da;

        /* skip if 'lpj->n' is equal to 'mini' */
        if ((lpj != NULL) && (lpj->n == mini)) {
            lpj = lpj->next;
        }

        /* calculate sum_cols[mini] */
        da = 0.0;
        for( ; lpj != NULL; lpj = lpj->next) {
            da += join[lpj->n];
        }
        sum_cols[mini] = da;

        /*
         * Clean up sum_rows[minj], sum_cols[minj] and score matrix
         * related with 'minj'.
         */
        sum_cols[minj] = sum_rows[minj] = 0.0;
        for(j=1; j<=last_seq-first_seq+1; ++j){
            tmat[minj][j] = tmat[j][minj] = join[j] = 0.0;
        }

/****/
    }						/*end main cycle**/

/******************************Last Cycle (3 Seqs. left)********************/

    nude = 1;

    for(lpi=tvalid[0].next; lpi != NULL; lpi = lpi->next){
        l[nude] = lpi->n;
        ++nude;
    }

    b1 = (tmat[l[1]][l[2]] + tmat[l[1]][l[3]] - tmat[l[2]][l[3]]) * 0.5;
    b2 =  tmat[l[1]][l[2]] - b1;
    b3 =  tmat[l[1]][l[3]] - b1;

    branch[1] = b1 - av[l[1]];
    branch[2] = b2 - av[l[2]];
    branch[3] = b3 - av[l[3]];

/* Reset tiny negative and positive branch lengths to zero */
    if( fabs(branch[1]) < 0.0001) branch[1] = 0.0;
    if( fabs(branch[2]) < 0.0001) branch[2] = 0.0;
    if( fabs(branch[3]) < 0.0001) branch[3] = 0.0;

    left_branch[last_seq-first_seq+1-2] = branch[1];
    left_branch[last_seq-first_seq+1-1] = branch[2];
    left_branch[last_seq-first_seq+1]   = branch[3];

    for(i=1; i<=last_seq-first_seq+1; i++){
        tree_description[last_seq-first_seq+1-2][i] = 0;
    }


    for(i=1; i<=3; ++i) {
       if( av[l[i]] > 0.0){
            for(k=last_seq-first_seq+1-3; k>=1; k--)
                if(tree_description[k][l[i]] == 1){
                    for(j=1; j<=last_seq-first_seq+1; j++)
                        if(tree_description[k][j] == 1)
                            tree_description[last_seq-first_seq+1-2][j] = i;
                        break;
                    }
        }
        else{
            tree_description[last_seq-first_seq+1-2][l[i]] = i;
        }
        if(i < 3) {;
        }
    }
    ckfree(sum_cols);
    ckfree(sum_rows);
    ckfree(join);
    ckfree(tvalid);
}


/*
 * NJ Balanced Guide Tree slow generation method
 */

void slow_nj_tree_bgt(char **tree_description){
    
    int l[4], nude, k;
    int nc, mini, minj, entra, i, j, ii, jj;
    double fnseqs, fnseqs2, sumd;
    double diq, djq, dij, d2r, dr, dio, djo, da;
    double tmin, total, dmin;
    double bi, bj, b1, b2, b3, branch[4];
    int typei, typej;		/* 0 = node; 1 = OTU */

    int cont=0;
    double average=0;
    int *eliminate;
    int *tkill_or;
    double fnseqs_or, fnseqs2_or;
    double original[nseqs+1][nseqs+1];

    fnseqs = fnseqs_or = (double)nseqs;


/* First Initialisation */
    mini = minj = entra = 0;
    left_branch = (double *) vcalloc((nseqs+2),sizeof(double));
    right_branch = (double *) vcalloc((nseqs+2),sizeof(double));
    tkill = (int *) vcalloc((nseqs+1), sizeof(int));
    av = (double *) vcalloc((nseqs+1),sizeof(double));
    eliminate = (int *) vcalloc((nseqs+1), sizeof(int));
    tkill_or = (int *) vcalloc((nseqs+1), sizeof(int));

    for (i = 1; i <= nseqs; ++i) {
        tmat[i][i] = av[i] = 0.0;
        tkill[i] = 0;
        eliminate[i]=0;
        for(j=1; j<=nseqs; ++j){
            original[i][j]=tmat[i][j];
        }
    }
/* Enter The Main Cycle */
        /* start main cycle */     

    for (nc = 1; nc <= (nseqs - 3); ++nc){

        sumd = 0.0;
        cont = 0;
        for (j = 2; j <= nseqs; ++j){
            for (i = 1; i < j; ++i) {
                tmat[j][i] = tmat[i][j];
                sumd = sumd + tmat[i][j];
                if(nc==1){
                    cont++;
                }
            }
        }
        if(nc==1){
            average=sumd/cont;
        }
        tmin = 99999.0;
        /* compute SMATij values and find the smallest one */
        for (jj = 2; jj <= nseqs; ++jj){
            if (tkill[jj] != 1){
                for (ii = 1; ii < jj; ++ii){
                    if (tkill[ii] != 1){
                        diq = djq = 0.0;
                        for (i = 1; i <= nseqs; ++i){
                            diq = diq + tmat[i][ii];
                            djq = djq + tmat[i][jj];
                        }
                        dij = tmat[ii][jj];
                        d2r = diq + djq - (2.0 * dij);
                        dr = sumd - dij - d2r;  
                        fnseqs2 = fnseqs - 2.0;
                        fnseqs2_or = fnseqs_or - 2.0;
                        total = d2r + fnseqs2_or * dij + dr * 2.0;
                        total = total / (2.0 * fnseqs2_or);
                        //total=double(int(total*1000.0+.5))/1000.0;
                        if (total < tmin) {
                            if (tmin - total > MY_EPS){
                                    tmin = total;
                            }
                            mini = ii;
                            minj = jj;
                        }
                    }
                }
            }
        }
    //std::cout<<"\n\nminj = "<<minj<<" --- "<<"mini = "<< mini <<" --- "<<"tmin = "<<tmin<<" --- "<<"dmin = "<<dmin<<"\n\n";
    //Fins aqui he fet la comparació de cada columna amb les altres a partir del calcul d'un escor i hem quedo amb
    //la comparació mínima i hem guardo els index de les columnes.

    /* compute branch lengths and print the results */

        dmin = tmat[mini][minj];
        if(fabs(sumd)<0.000001){
            sumd = 0.0;
        }
        if(sumd==0.0  || (mini == 0 && minj == 0)){
            dmin=999999.99999;
        }

        if(dmin<=average || entra !=0){
            
            if(dmin<=average){
                entra=0;
            }
            dio = djo = 0.0;
            eliminate[minj]=1;
            for (i = 1; i <= nseqs; ++i) {
                    dio = dio + tmat[i][mini]; //Calculo el sumatori de les columnes minimes i i j.
                    djo = djo + tmat[i][minj];
            }
            dio = (dio - dmin) / fnseqs2;
            djo = (djo - dmin) / fnseqs2;
            bi = (dmin + dio - djo) * 0.5;
            bj = dmin - bi;
            bi = bi - av[mini];
            bj = bj - av[minj];

            if (av[mini] > 0.0)
                typei = 0;
            else
                    typei = 1;
            if (av[minj] > 0.0)
                typej = 0;
            else
                typej = 1;
            /* Set negative branch lengths to zero;
             * also set any tiny positive branch lengths to zero. */
            if (fabs(bi) < 0.0001)
                bi = 0.0;
            if (fabs(bj) < 0.0001)
                bj = 0.0;
            left_branch[nc] = bi;
            right_branch[nc] = bj;

            for (i = 1; i <= nseqs; i++)
                    tree_description[nc][i] = 0;
            if (typei == 0) {
                for (i = nc - 1; i >= 1; i--){
                    if (tree_description[i][mini] == 1) {
                        for (j = 1; j <= nseqs; j++){
                            if (tree_description[i][j] == 1){
                                tree_description[nc][j] = 1;
                            }
                        }
                        break;
                    }
                }
            }
            else{
                    tree_description[nc][mini] = 1;
            }
            if (typej == 0) {
                for (i = nc - 1; i >= 1; i--){
                    if (tree_description[i][minj] == 1) {
                        for (j = 1; j <= nseqs; j++){
                            if (tree_description[i][j] == 1)
                                tree_description[nc][j] = 1;
                        }
                        break;
                    }
                }
            }
            else{
                tree_description[nc][minj] = 1;
            }

                /* Here is where the -0.00005 branch lengths come from for 3 or more identical seqs. */
            if (dmin <= 0.0)
                dmin = 0.000001;
            av[mini] = dmin * 0.5;
            /* re-initialisation */
            if(entra==0){
                fnseqs = fnseqs - 2.0;
            }
            else{
                fnseqs = fnseqs - 1.0;
            }
            fnseqs_or = fnseqs_or - 1.0;
            tkill[minj] = 1;
            tkill_or[minj] =1;
            for (j = 1; j <= nseqs; ++j){
                if (tkill_or[j] != 1) {
                    da = (original[mini][j] + original[minj][j]) * 0.5;
                    if ((mini - j) < 0){
                        tmat[mini][j] = da;
                        original[mini][j] = da;
                    }
                    if ((mini - j) > 0){
                        tmat[j][mini] = da;
                        original[j][mini] = da;
                    }
                }
            }
            for (j = 1; j <= nseqs; ++j){
                tmat[minj][j] = tmat[j][minj] = original[minj][j] = original[j][minj] = 0.0;
                if(nc<(nseqs - 3) && entra==0){
                    tmat[mini][j] = tmat[j][mini] =  0.0;
                    tkill[mini]= 1;
                }
                else if(nc==(nseqs-3)){
                    for(i=1; i<=nseqs; i++){
                        tmat[j][i]=original[j][i];
                    }
                    if(eliminate[j]==0){
                        tkill[j]=0;
                    }
                }
            }
        }
        else {
            entra++;
            fnseqs = fnseqs_or;
            for(i=1; i<=nseqs; i++){
                for(j=1; j<=nseqs; j++){
                    tmat[i][j]=original[i][j];
                }
                if(eliminate[i]==0){
                    tkill[i]=0;
                }
            }
            nc=nc-1;
        }
    }
    /* end main cycle */
/* Last Cycle (3 seqs. left) Brankes a esquerres*/
    nude = 1;
    for (i = 1; i <= nseqs; ++i){
        if (tkill[i] != 1) {
            l[nude] = i;
            nude = nude + 1;
        }
    }
    b1 = (tmat[l[1]][l[2]] + tmat[l[1]][l[3]] - tmat[l[2]][l[3]]) * 0.5;
    b2 = tmat[l[1]][l[2]] - b1;
    b3 = tmat[l[1]][l[3]] - b1;
    branch[1] = b1 - av[l[1]];
    branch[2] = b2 - av[l[2]];
    branch[3] = b3 - av[l[3]];


    /* Reset tiny negative and positive branch lengths to zero */
    if (fabs(branch[1]) < 0.0001)
        branch[1] = 0.0;
    if (fabs(branch[2]) < 0.0001)
        branch[2] = 0.0;
    if (fabs(branch[3]) < 0.0001)
        branch[3] = 0.0;

    left_branch[nseqs - 2] = branch[1];
    left_branch[nseqs - 1] = branch[2];
    left_branch[nseqs] = branch[3];

    for (i = 1; i <= nseqs; i++)
        tree_description[nseqs - 2][i] = 0;
    for (i = 1; i <= 3; ++i) {
        if (av[l[i]] > 0.0){
            for (k = nseqs - 3; k >= 1; k--){
                if (tree_description[k][l[i]] == 1) {
                    for (j = 1; j <= nseqs; j++){
                        if (tree_description[k][j] == 1)
                            tree_description[nseqs - 2][j] = i;
                    }
                        break;
                }
            }
        }
        else{
            tree_description[nseqs - 2][l[i]] = i;
        }
    }
    vfree(tkill_or);
    vfree(eliminate);
}

/*
 * NJ Balanced Guide Tree fast generation method
 */

void fast_nj_tree_bgt(char **tree_description){

    int l[4], nude, k;
        int nc, mini, minj, i, j, ii, jj;
        double fnseqs, fnseqs2 = 0, sumd;
        double diq, djq, dij, dio, djo, da;
        double tmin, total, dmin;
        double bi, bj, b1, b2, b3, branch[4];
        int typei, typej;		/* 0 = node; 1 = OTU */

/* IMPROVEMENT 1, STEP 0: declare  variables */
        double *sum_cols, *sum_rows, *join;

/* IMPROVEMENT 2, STEP 0: declare  variables */
        int loop_limit;
        typedef struct _ValidNodeID {
            int n;
            struct _ValidNodeID *prev;
            struct _ValidNodeID *next;
        } ValidNodeID;
        ValidNodeID *tvalid, *lpi, *lpj, *lpii, *lpjj, *lp_prev, *lp_next;

        double cont=0.0;
        int entra=0, nseq=0;
        double average=0.0;
        double original[nseqs+1][nseqs+1];
        int *eliminate;
        int *tkill_or;
        double fnseqs_or, fnseqs2_or;
        double *sum_cols_or;
        double *sum_rows_or;

        /*
         * correspondence of the loop counter variables.
         *   i .. lpi->n,  ii .. lpii->n
         *   j .. lpj->n,  jj .. lpjj->n
         */
        fnseqs = fnseqs_or = (double)(last_seq - first_seq + 1);

/* First Initialisation */
        if (fnseqs == 2)
                return;

        mini = minj = 0;
        left_branch = (double *) ckalloc((nseqs + 2) * sizeof(double));
        right_branch = (double *) ckalloc((nseqs + 2) * sizeof(double));
        tkill = (int *) ckalloc((nseqs + 1) * sizeof(int));
        av = (double *) ckalloc((nseqs + 1) * sizeof(double));

/* IMPROVEMENT 1, STEP 1: Allocate memory */
        sum_cols = (double *) ckalloc((nseqs + 1) * sizeof(double));
        sum_rows = (double *) ckalloc((nseqs + 1) * sizeof(double));
        join = (double *) ckalloc((nseqs + 1) * sizeof(double));
/* IMPROVEMENT 2, STEP 1: Allocate memory */
        tvalid = (ValidNodeID *) ckalloc((nseqs + 1) * sizeof(ValidNodeID));
        /* tvalid[0] is special entry in array, it points a header of valid entry list */
        tvalid[0].n = 0;
        tvalid[0].prev = 0;
        tvalid[0].next = &tvalid[1];

        tkill_or = (int *) ckalloc((nseqs + 1) * sizeof(int));
        eliminate = (int *) ckalloc((nseqs + 1) * sizeof(int));
        av = (double *) ckalloc((nseqs + 1) * sizeof(double));

/* IMPROVEMENT 1, STEP 1: Allocate memory */
        sum_cols_or = (double *) ckalloc((nseqs + 1) * sizeof(double));
        sum_rows_or = (double *) ckalloc((nseqs + 1) * sizeof(double));

/* IMPROVEMENT 2, STEP 2: Construct and initialize the entry chain list */
        for (i = 1, loop_limit = last_seq - first_seq + 1, lpi = &tvalid[1], lp_prev = &tvalid[0],
                lp_next = &tvalid[2]; i <= loop_limit; ++i, ++lpi, ++lp_prev, ++lp_next) {
                tmat[i][i] = av[i] = 0.0;
                tkill[i] = 0;
                tkill_or[i] = 0;
                eliminate[i]=0;
                lpi->n = i;
                lpi->prev = lp_prev;
                lpi->next = lp_next;
/* IMPROVEMENT 1, STEP 2: Initialize arrays */
                sum_cols[i] = sum_rows[i] = sum_cols_or[i] = sum_rows_or[i] = join[i] = 0.0;
                for(j=1; j<=loop_limit; j++){
                    original[i][j] = tmat[i][j];
                }
        }

        tvalid[loop_limit].next = 0;

/* IMPROVEMENT 1, STEP 3: Calculate the sum of score value that is sequence[i] to others. */
        sumd = 0.0;
        for (lpj = tvalid[0].next; lpj != 0; lpj = lpj->next) {
                double tmp_sum = 0.0;
                j = lpj->n;
                /* calculate sum_rows[j] */
                for (lpi = tvalid[0].next; lpi->n < j; lpi = lpi->next) {
                        i = lpi->n;
                        tmp_sum += tmat[i][j];
                }
                sum_rows[j] = tmp_sum;
                sum_rows_or[j] = tmp_sum;
                tmp_sum = 0.0;
                /* set lpi to that lpi->n is greater than j */
                if (lpi != 0 && lpi->n == j)
                        lpi = lpi->next;
                /* calculate sum_cols[j] */
                for (; lpi != 0; lpi = lpi->next) {
                        i = lpi->n;
                        tmp_sum += tmat[j][i];
                        cont++;
                }
                sum_cols[j] = tmp_sum;
                sum_cols_or[j] = tmp_sum;
        }

/* Enter The Main Cycle */
        for (nc = 1, loop_limit = (last_seq - first_seq + 1 - 3); nc <= loop_limit; ++nc) {
                sumd = 0.0;
/* IMPROVEMENT 1, STEP 4: use sum value */

                for (lpj = tvalid[0].next; lpj != 0; lpj = lpj->next){
                    if(tkill[lpj->n] != 1){
                        sumd += sum_cols[lpj->n];
                    }
                }
                if(nc==1){
                    average=sumd/cont;
                }

/* IMPROVEMENT 3, STEP 0: multiply tmin and 2*fnseqs2 */
                fnseqs2 = fnseqs - 2.0;	/* Set fnseqs2 at this point. */
                fnseqs2_or = fnseqs_or - 2.0;
                tmin = 99999.0 * 2.0 * fnseqs2;

                /* compute SMATij values and find the smallest one */
                mini = minj = 0;

                /* jj must starts at least 2 */
                if (tvalid[0].next != 0 && tvalid[0].next->n == 1){
                        lpjj = tvalid[0].next->next;
                }
                else{
                        lpjj = tvalid[0].next;
                }
                for (; lpjj != 0; lpjj = lpjj->next) {
                        jj = lpjj->n;
                        if(tkill[jj]!=1){
                            for (lpii = tvalid[0].next; lpii->n < jj; lpii = lpii->next) {
                                    ii = lpii->n;
                                    if(tkill[ii]!=1){
                                        diq = djq = 0.0;
        /* IMPROVEMENT 1, STEP 4: use sum value */
                                        diq = sum_cols[ii] + sum_rows[ii];
                                        djq = sum_cols[jj] + sum_rows[jj];
                                        /* always ii < jj in this point.
                                         * Use upper triangle of score matrix. */
                                        dij = tmat[ii][jj];
        /* IMPROVEMENT 3, STEP 1: fnseqs2 is already calculated. */
                                        //fnseqs2 = fnseqs - 2.0

        /* IMPROVEMENT 4: transform the equation */
        /* --------------------------------------------------------------- *
         * OPTIMIZE of expression 'total = d2r + fnseqs2 * dij + dr * 2.0' *
         *   total = d2r + fnseq2*dij + 2.0*dr                             *
         *         = d2r + fnseq2*dij + 2(sumd - dij - d2r)                *
         *         = d2r + fnseq2*dij + 2*sumd - 2*dij - 2*d2r             *
         *         =       fnseq2*dij + 2*sumd - 2*dij - 2*d2r + d2r       *
         *         = fnseq2*dij + 2*sumd - 2*dij - d2r                     *
         *         = fnseq2*dij + 2*sumd - 2*dij - (diq + djq - 2*dij)     *
         *         = fnseq2*dij + 2*sumd - 2*dij - diq - djq + 2*dij       *
         *         = fnseq2*dij + 2*sumd - 2*dij + 2*dij - diq - djq       *
         *         = fnseq2*dij + 2*sumd  - diq - djq                      *
         * --------------------------------------------------------------- */
                                        total = fnseqs2_or * dij + 2.0 * sumd - diq - djq;
                                        //total=double(int(total*1000.0+.5))/1000.0;

        /* IMPROVEMENT 3, STEP 2: abbrevlate the division on comparison between total and tmin. */
                                        //total = total / (2.0 * fnseqs2);
                                        if (total < tmin) {
                                                tmin = total;
                                                mini = ii;
                                                minj = jj;
                                        }
                                    }
                            }
                        }
                }
                 
                /* MEMO: always ii < jj in avobe loop, so mini < minj */

                /* compute branch lengths and print the results */
                dio = djo = 0.0;
                dmin = tmat[mini][minj];          
                if(fabs(sumd)<0.000001){
                    sumd = 0.0;
                }
                if(sumd == 0.0 || (mini == 0 && minj == 0)){
                    dmin=999999.99999;
                }

                if(dmin<=average || entra!=0){
                    if(dmin<=average){
                        entra=0;
                    }

/* IMPROVEMENT 1, STEP 4: use sum value */
                    dio = sum_cols[mini] + sum_rows[mini];
                    djo = sum_cols[minj] + sum_rows[minj];

                    dio = (dio - dmin) / fnseqs2;
                    djo = (djo - dmin) / fnseqs2;
                    bi = (dmin + dio - djo) * 0.5;
                    bj = dmin - bi;
                    bi = bi - av[mini];
                    bj = bj - av[minj];
                    if (av[mini] > 0.0){
                            typei = 0;
                    }
                    else {
                            typei = 1;
                    }
                    if (av[minj] > 0.0){
                            typej = 0;
                    }
                    else {
                            typej = 1;

                    }
                    /* Set negative branch lengths to zero;
                     * also set any tiny positive branch lengths to zero. */
                    if (fabs(bi) < 0.0001){
                            bi = 0.0;
                    }
                    if (fabs(bj) < 0.0001){
                            bj = 0.0;
                    }
                    left_branch[nc] = bi;
                    right_branch[nc] = bj;
                    for (i = 1; i <= last_seq - first_seq + 1; i++){
                            tree_description[nc][i] = 0;
                    }
                    if (typei == 0) {
                            for (i = nc - 1; i >= 1; i--){
                                    if (tree_description[i][mini] == 1) {
                                            for (j = 1; j <= last_seq - first_seq + 1; j++){
                                                    if (tree_description[i][j] == 1){
                                                            tree_description[nc][j] = 1;
                                                    }
                                            }
                                            break;

                                    }
                            }
                    }
                    else {
                            tree_description[nc][mini] = 1;
                    }
                    if (typej == 0) {
                            for (i = nc - 1; i >= 1; i--){
                                    if (tree_description[i][minj] == 1) {
                                            for (j = 1; j <= last_seq - first_seq + 1; j++){
                                                    if (tree_description[i][j] == 1){
                                                            tree_description[nc][j] = 1;

                                                    }
                                            }
                                            break;
                                    }

                            }
                    }
                    else {
                            tree_description[nc][minj] = 1;
                    }
                    /* Here is where the -0.00005 branch lengths come from for 3 or more identical seqs. */
                    if (dmin <= 0.0){
                            dmin = 0.000001;
                    }
                    av[mini] = dmin * 0.5;

                    /* re-initialisation */
                    if(entra==0){
                        fnseqs = fnseqs - 2.0;
                    }
                    else{
                        fnseqs = fnseqs - 1.0;
                    }
                    fnseqs_or = fnseqs_or - 1.0;
                    eliminate[minj]=1;
                    tkill[minj] = 1;
                    tkill_or[minj] = 1;

    /* IMPROVEMENT 2, STEP 3: Remove tvalid[minj] from chain list. */
                    (tvalid[minj].prev)->next = tvalid[minj].next;
                    if (tvalid[minj].next != 0){
                        (tvalid[minj].next)->prev = tvalid[minj].prev;
                    }
                    tvalid[minj].prev = tvalid[minj].next = 0;
    /* IMPROVEMENT 1, STEP 5: re-calculate sum values. */
                    for (lpj = tvalid[0].next; lpj != 0; lpj = lpj->next) {
                            double tmp_di = 0.0;
                            double tmp_dj = 0.0;

                            j = lpj->n;
                            /* subtrace a score value related with 'minj' from sum arrays. */
                            if (j < minj) {
                                    tmp_dj = original[j][minj];
                                    if(tkill[j]!=1){
                                        sum_cols[j] -= tmp_dj;
                                    }
                                    sum_cols_or[j] -= tmp_dj;                                 
                            }
                            else if (j > minj) {
                                    tmp_dj = original[minj][j];
                                    if(tkill[j]!=1){
                                       sum_rows[j] -= tmp_dj;
                                    }
                                    sum_rows_or[j] -= tmp_dj;                             
                            }
                            /* nothing to do when j is equal to minj. */
                            /* subtrace a score value related with 'mini' from sum arrays. */
                            if (j < mini) {
                                    tmp_di = original[j][mini];
                                    if(tkill[j]!=1){
                                        sum_cols[j] -= tmp_di;
                                    }
                                    sum_cols_or[j] -= tmp_di;
                            }
                            else if (j > mini) {
                                    tmp_di = original[mini][j];
                                    if(tkill[j]!=1){
                                        sum_rows[j] -= tmp_di;
                                    }
                                    sum_rows_or[j] -= tmp_di;
                            }

                            /* nothing to do when j is equal to mini. */
                            /* calculate a score value of the new inner node;
                             * then, store it temporary to join[] array. */
                            join[j] = (tmp_dj + tmp_di) * 0.5;
                    }
                    /*
                     * 1) Set the score values (stored in join[]) into the matrix, row/column position is 'mini'.
                     * 2) Add a score value of the new inner node to sum arrays.
                     */
                    for (lpj = tvalid[0].next; lpj != 0; lpj = lpj->next) {
                            j = lpj->n;
                            if (j < mini) {
                                    original[j][mini] = join[j];
                                    tmat[j][mini] = join[j];
                                    sum_cols_or[j] += join[j];
                            }
                            else if (j > mini) {
                                    original[mini][j] = join[j];
                                    tmat[mini][j] = join[j];
                                    sum_rows_or[j] += join[j];
                            }

                            /* nothing to do when j is equal to mini. */
                    }
                    /* re-calculate sum_rows[mini], sum_cols[mini]. */
                    sum_cols[mini] = sum_rows[mini] = sum_cols_or[mini] = sum_rows_or[mini] = 0.0;

                    /* calculate sum_rows[mini] */
                    da = 0.0;
                    for (lpj = tvalid[0].next; lpj->n < mini; lpj = lpj->next){
                            da += join[lpj->n];
                    }
                    sum_rows_or[mini] = da;

                    /* skip if 'lpj->n' is equal to 'mini' */
                    if (lpj != 0 && lpj->n == mini){
                        lpj = lpj->next;
                    }

                    /* calculate sum_cols[mini] */
                    da = 0.0;
                    for (; lpj != 0; lpj = lpj->next){
                        da += join[lpj->n];
                    }
                    sum_cols_or[mini] = da;

                    /* Clean up sum_rows[minj], sum_cols[minj] and score matrix related with 'minj'. */
                    sum_cols[minj] = sum_rows[minj] = sum_rows_or[minj] = sum_cols_or[minj] = 0.0;
                    for (j = 1; j <= last_seq - first_seq + 1; ++j){
                        tmat[minj][j] = tmat[j][minj] = original[minj][j] = original[j][minj] = join[j] = 0.0;
                        if(nc<(nseqs - 3) && entra==0){
                            tmat[mini][j] = tmat[j][mini] = 0.0;
                            sum_cols[mini] = sum_rows[mini] = 0.0;
                            tkill[mini]=1;
                        }
                    }
                }
                else {
                    
                    entra++;
                    fnseqs = fnseqs_or;
                    for(i=1; i<=last_seq - first_seq + 1; i++){
                        sum_cols[i] = sum_cols_or[i];
                        sum_rows[i] = sum_rows_or[i];
//                        if(eliminate[i]==1){
//                            std::cout<<"\nEntra a i: "<<i<<"\n";
                        for(j=1; j<=last_seq - first_seq + 1; j++){
                            tmat[i][j]=original[i][j];

                        }

                        if(eliminate[i]==0){                  
                            tkill[i]=0;
                        }
                    }

                    nc=nc-1;

                }

        }
        /* end main cycle */

/* Last Cycle (3 seqs. left) */
        nude = 1;
        for (lpi = tvalid[0].next; lpi != 0; lpi = lpi->next){
                l[nude] = lpi->n;
                ++nude;
        }

        b1 = (tmat[l[1]][l[2]] + tmat[l[1]][l[3]] - tmat[l[2]][l[3]]) * 0.5;
        b2 = tmat[l[1]][l[2]] - b1;
        b3 = tmat[l[1]][l[3]] - b1;
        branch[1] = b1 - av[l[1]];
        branch[2] = b2 - av[l[2]];
        branch[3] = b3 - av[l[3]];

        /* Reset tiny negative and positive branch lengths to zero */
        if (fabs(branch[1]) < 0.0001)
                branch[1] = 0.0;
        if (fabs(branch[2]) < 0.0001)
                branch[2] = 0.0;
        if (fabs(branch[3]) < 0.0001)
                branch[3] = 0.0;

        left_branch[last_seq - first_seq + 1 - 2] = branch[1];
        left_branch[last_seq - first_seq + 1 - 1] = branch[2];
        left_branch[last_seq - first_seq + 1] = branch[3];
        for (i = 1; i <= last_seq - first_seq + 1; i++)
                tree_description[last_seq - first_seq + 1 - 2][i] = 0;
        for (i = 1; i <= 3; ++i) {
                if (av[l[i]] > 0.0) {
                        for (k = last_seq - first_seq + 1 - 3; k >= 1; k--)
                                if (tree_description[k][l[i]] == 1) {
                                        for (j = 1; j <= last_seq - first_seq + 1; j++)
                                                if (tree_description[k][j] == 1)
                                                        tree_description[last_seq - first_seq + 1 - 2][j] = i;
                                        break;
                                }
                }
                else
                        tree_description[last_seq - first_seq + 1 - 2][l[i]] = i;
        }
//
        ckfree(sum_cols);
        ckfree(sum_rows);
        ckfree(sum_cols_or);
        ckfree(sum_rows_or);
        ckfree(tkill_or);
        ckfree(eliminate);
        ckfree(join);
        ckfree(tvalid);

}

char* identify_tree(NT_node root, char *tree_id, int nseqs){
    
    NT_node LL, RL;
    
    LL = root->left;
    RL = root->right;

    if(!LL->isseq && !RL->isseq){
        //printf("Node Left\n");
        LL->level = (LL->parent)->level + 1;
        tree_id = identify_tree(LL, tree_id, nseqs);
         //printf("Node Right\n");
        RL->level = (RL->parent)->level + 1;
        tree_id = identify_tree(RL,  tree_id, nseqs);
    }
    else if(LL->isseq && !RL->isseq){  
        //printf("Leaf left: %s\n", LL->name);
        LL->level = (LL->parent)->level + 1;
        sprintf(tree_id, "%s%d%dl", tree_id, LL->seq_id, LL->level);
        //printf("Node Right\n");
        RL->level = (RL->parent)->level + 1;
        tree_id = identify_tree(RL, tree_id, nseqs);
    }
    else if(!LL->isseq && RL->isseq){
        //printf("Node Left\n");
        LL->level = (LL->parent)->level + 1;     
        tree_id = identify_tree(LL, tree_id, nseqs);
        //printf("Leaf right: %s\n", RL->name);
        RL->level = (RL->parent)->level + 1;
        sprintf(tree_id, "%s%d%dl", tree_id, RL->seq_id, RL->level);
    }
    else if(LL->isseq && RL->isseq){ 
        //printf("Leaf left: %s\n", LL->name);
        LL->level = (LL->parent)->level + 1;
        sprintf(tree_id, "%s%d", tree_id, LL->seq_id);
        //printf("Leaf right: %s\n", RL->name);
        RL->level = (RL->parent)->level + 1;
        sprintf(tree_id, "%s%d%dl", tree_id, RL->seq_id, RL->level);      
    }
    
    return tree_id;
}

int search_tree_id(char** tree_id_list, char* tree_id, int ntreeid){
    int i=0;

    for(i=0; i<ntreeid; i++){
        if(strm(tree_id_list[i], tree_id)){
            return 1;
        }
    }
    return 0;
}

int count_entries_file(char *fname, int nseq){

    FILE *fp;
    int entries=0;
    char *tree_id;

    tree_id = (char *) vcalloc(1000, sizeof(char));
    fp = openfile(fname, "r");

    while(fscanf(fp, "%s", tree_id) != EOF){
        entries++;
    }
    vfree(tree_id);
    return entries;
}

char** file2tree_id_list(char *fname, int ntrees, int nseqs, int entries){

    FILE *fp;
    int i=0, tree_id_lenght=0, new_ntrees=0;
    char *tree_id;
    char **tree_id_list;

    tree_id_lenght = (4 * nseqs);
    new_ntrees = ntrees + entries;
    tree_id = (char *) vcalloc(1000, sizeof(char));

    tree_id_list = declare_char(new_ntrees, tree_id_lenght);
    fp = openfile(fname, "r");

    while(fscanf(fp, "%s", tree_id) != EOF){
        sprintf(tree_id_list[i], "%s", tree_id);
        i++;
    }

    if(entries != i){
        fprintf(stderr, "ERROR - Bad number of tree id entries");
    }

    fclose(fp);
    vfree(tree_id);
    return tree_id_list;
}

char* root_unrooted_tree(char *treename, int ntree, char *retree_bin){

    FILE *fp;
    char *intmpname, *outtmpname, *treename2;
    char *cmd;
    int ret = 0; 
    
    intmpname = (char *) vcalloc(FILENAMELEN, sizeof(char));
    outtmpname = (char *) vcalloc(FILENAMELEN, sizeof(char));
    treename2 = (char *) vcalloc(FILENAMELEN, sizeof(char));
    cmd = (char *) vcalloc(ALLPATH, sizeof(char));
    
    sprintf(intmpname, "./input_%d.run", ntree);
    sprintf(outtmpname, "./output_%d.txt", ntree);
    sprintf(treename2, "%s.rooted", treename);
    fp = fopen("outtree", "w");
    fclose(fp);
    fp = fopen(intmpname, "w");
    fprintf(fp, "Y\n");
    fprintf(fp, "%s\n", treename);
    fprintf(fp, "W\n");
    fprintf(fp, "F\n");
    fprintf(fp, "%s\n", treename2);
    fprintf(fp, "R\n");
    fprintf(fp, "Q\n");
    
    fclose(fp);

    sprintf(cmd, "%s < %s > %s", retree_bin, intmpname, outtmpname);
    ret=system(cmd);

    
    remove(intmpname);
    remove(outtmpname);
    //remove("outtree");
    //vfree(treename2);
    vfree(intmpname);
    vfree(outtmpname);
    vfree(cmd);
  
    return treename2;
    
}

