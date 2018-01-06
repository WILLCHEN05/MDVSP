#include<iostream>
#include<fstream>
#include<ostream>
#include<istream>
#include<cstring>
#include<stack>
using namespace std;

const double eps=0.0000001;

bool eq(const double& x, const double& y){
   return (x>=y-eps && x<=y+eps);
}

bool gt(const double& x, const double& y){
   	return (x-eps>y);
}

bool low(const double& x, const double& y){
    return (x+eps<y);	
}

bool loweq(const double& x, const double& y){
    return (x+eps<y) || eq(x,y);	
}

const int NR_LINES=500;
const int NR_VAR=50000;

struct simplex {
  int n,m;
  int n1, m1, semn[NR_VAR];
  double b[500][NR_VAR];
  double a[500][NR_VAR];
  double bv[500];
  double z[NR_VAR];	
  int has_solution=1;
  
  //transform to standard  form (add slack variables)
  void canonical_to_standard() {
  	
  	int maux=0;
  	for (int i=1; i<=m; ++i) if (semn[i]!=0) ++maux;
  	int iaux=0;
  	
  	for (int i=1; i<=m; ++i)
	  if (semn[i]==0) { a[i][n+maux+1]=bv[i]; }
	  else {
	   ++iaux;
	   if (bv[i]>=0) {
	     a[i][n+iaux]=1.0;
	     a[i][n+maux+1]=bv[i];   
	   }
	  else {
	 	a[i][n+iaux]=-1.0;
	    a[i][n+maux+1]=-bv[i]; 
	  }
     }
   
    n=n+maux; 
  }
  
  //read the problem in canonical form Ax<=b
  void read(istream &in) {
  	
  	in>>m>>n;
  	
  	//coeficientii functiei obiectiv
	for (int i=1; i<=n; ++i) in>>z[i];
	for (int i=1; i<=m; ++i) in>>semn[i];
	
	for (int i=1; i<=m; ++i) {
	 for (int j=1; j<=n; ++j) in>>a[i][j];
	 cin>>bv[i];
    }
	 
	canonical_to_standard();
	 
	for (int i=1; i<=n; ++i) a[m+1][i]=z[i];
  }
  
  void build_from_memory(int n2, int m2, double z1[], double a1[][NR_VAR], double bv1[], int semn1[]) {
    n=n2; m=m2;
    memset(z,0,sizeof(z));
    for (int i=1; i<=n; ++i) z[i]=z1[i];
    
    memset(a,0,sizeof(a));
    for (int i=1; i<=m; ++i)
     for (int j=1; j<=n; ++j) a[i][j]=a1[i][j];
    
    memset(bv,0,sizeof(bv));
    for (int i=1; i<=m; ++i) { bv[i]=bv1[i]; semn[i]=semn1[i]; }
    
    canonical_to_standard();
 	 
	for (int i=1; i<=n; ++i) a[m+1][i]=z[i];
  }
  
  void printState(ostream &out, double (*a)[NR_VAR], int n, int m) {
	out<<"\n";
  	for (int i=1; i<=m+1; ++i) {
  		out<<a[i][0]<<" | ";
  		for (int j=1; j<=n; ++j) out<<a[i][j]<<" ";
  		out<<"| "<<a[i][n+1]<<"\n";
  	}
	out<<"\n";
  }
  
  int solve_simplex_table(ostream &out, double (*a)[NR_VAR], int n, int m) {
  	
  	out<<"Enter solve simplex table:\n";
  	out<<m<<" "<<n<<"\n";
  	int it=0;
	while (1) {
	  
	  //check if the solution is optim
	  bool optim=1;
	  int col;
	  double maxaux=1000000000;
	  
	  if (it%2==0) {
	  for (int i=1; i<=n; ++i)
	   if (low(a[m+1][i],0)) { 
	       optim=0; 
		   double min=1000000000;
	       int lin=-1;
	  
	       for (int j=1; j<=m; ++j)
	        if ( gt(a[j][i],0) && loweq(a[j][n+1]/a[j][i],min) ){
	   	    if ( low(a[j][n+1]/a[j][i],min) ) { lin=j; min=a[j][n+1]/a[j][i]; }
	   	    else if ( (int)a[j][0]<(int)a[lin][0] ) lin=j;
	      }
	      
	      double inc=a[lin][n+1]*a[m+1][i]/a[lin][i];
	      
	      if (inc<=maxaux) { maxaux=inc; col=i; }
	    }
      }
      else {
      	for (int i=1; i<=n; ++i)
         if (low(a[m+1][i],0))  { optim=0; col=i; break; }	
      }
	  	
	  if (optim)  { out<<"Reached final state:\n"; break; }
	  
	  ++it;
	  if (it%100==0) { out<<"Iteration #"<<it<<"\n"; out<<"current z="<<a[m+1][n+1]<<"\n"; }
	  //printState(out, a, n, m);
	  
	  double min=1000000000;
	  int lin=-1;
	  
	  for (int i=1; i<=m; ++i)
	   if ( gt(a[i][col],0) && loweq(a[i][n+1]/a[i][col],min) ){
	   	  if ( low(a[i][n+1]/a[i][col],min) ) { lin=i; min=a[i][n+1]/a[i][col]; }
	   	  else if ( (int)a[i][0]<(int)a[lin][0] ) lin=i;
	   }
	   
	  if (lin==-1) { 
	     out<<"Infinite Optimum:\n"; 
		 return -1; //unbounded
	  }
	  
	  //out<<"Pivot: "<<lin<<" "<<col<<"\n";
	  
	  //pivot location is (lin,col)
	  for (int i=1; i<=m+1; ++i)
	   for (int j=1; j<=n+1; ++j)
	    if (i!=lin && j!=col) a[i][j]=a[i][j]-a[lin][j]*a[i][col]/a[lin][col];
	  
	  for (int i=1; i<=n+1; ++i) if (i!=col) a[lin][i]/=a[lin][col];
	  for (int i=1; i<=m+1; ++i) a[i][col]=0;
	  
	  a[lin][col]=1;
	  
	  a[lin][0]=col;
		
	}
	
	//printState(out, a, n, m);	
  	return 1; //optimal solution exists
  }
  
  void normalize_base(double (*a)[NR_VAR], int n, int m) {
     //replace base variables with non basic so the reduced cost of basic variables will become 0	
	 for (int i=1; i<=m; ++i)
	  if (!eq(a[m+1][(int)a[i][0]],0)) {
	    double mult = -a[m+1][(int)a[i][0]];
	    for (int j=1; j<=n+1; ++j)
	      a[m+1][j]+=a[i][j]*mult;
      }
  }
  
  bool buildBase() {
  	 int k=0;
     for (int col=1; col<=n; ++col)
	  if (eq(a[m+1][col],0)) {
	     int lin, nr=0;
	     for (int j=1; j<=m; ++j)
	      if (!eq(a[j][col],0)) { lin=j; ++nr; }
		  
		 if (nr==1 && eq(a[lin][col],1) && (int)a[lin][0]==0) { a[lin][0]=col; ++k; } 
      }
  	
  	 if (k==m) return 1;
  	 else { 
  	 	//add artificial variables and build table for phase 1
  	 	memset(b,0,sizeof(b));
  	 	for (int i=1; i<=m; ++i)
  	 	 for (int j=0; j<=n; ++j) b[i][j]=a[i][j];
  	 	 
  	 	n1=n; m1=m;
		   
		for (int i=1; i<=m; ++i)
		 if ((int)b[i][0]==0) {
		    ++n1;
		    b[i][n1]=1;
		    b[i][0]=n1;
		    b[m+1][n1]=1;
	     }
	     
	    for (int i=1; i<=m; ++i) b[i][n1+1]=a[i][n+1];
	    
	    normalize_base(b, n1, m1);
	    return 0;
  	 }
  }
  
  bool solve_phase1(ostream &out) {	
  	  solve_simplex_table(out,b,n1,m1);
  	  
  	  if (!eq(b[m1+1][n1+1],0)) { has_solution=0; return 0; }
  	  
  	  //remove artificial variables from base
  	  for (int i=1; i<=m; ++i)
  	   if ((int)b[i][0]>n) {
  	   	 int idx=-1;
  	   	 for (int j=1; j<=n; ++j)
  	   	  if (!eq(b[i][j],0)) { idx=j; break; }
  	   	  
  	   	 if (idx==-1) b[i][0]=-1;
  	   	 else {
  	   	   	  //pivot location is (i,idx)
	          for (int s=1; s<=m+1; ++s)
	           for (int t=1; t<=n1+1; ++t)
	            if (s!=i && t!=idx) b[s][t]=b[s][t]-b[i][t]*b[s][idx]/b[i][idx];
	  
	          for (int s=1; s<=n1+1; ++s) if (s!=idx) b[i][s]/=b[i][idx];
	          for (int s=1; s<=m1+1; ++s) b[s][idx]=0;
	  
	          b[i][idx]=1;
	          b[i][0]=idx;
  	   	 }
  	   }
  	   
  	  // rebuild a for second phase
  	  m=0;
  	  for (int i=1; i<=m1; ++i)
  	   if ((int)b[i][0]!=-1) {
  	   	++m;
  	    for (int j=0; j<=n; ++j) a[m][j]=b[i][j];
  	    a[m][n+1]=b[i][n1+1];
  	   }
  	   
  	  for (int i=1; i<=n+1; ++i) a[m+1][i]=z[i];
  	  
  	  normalize_base(a,n,m);
  	  
  	  return 1;
  }
  
  int solve(ostream &out) {
  	
    bool ok = buildBase();
	
	if (ok) return solve_simplex_table(out, a, n, m);
	else {
	  out<<"phase 1 required\n";
	  ok = solve_phase1(out);
	  if (ok) solve_simplex_table(out, a, n, m);
	  else { out<<"No feasable solutions exists\n"; return 2; }
    }
  }
};


struct ILP {
  int n,m,semn[NR_VAR];
  double a[500][NR_VAR];
  double bv[500];
  double z[NR_VAR];
  bool is_integer[NR_VAR];
  static double x[NR_VAR];
  static int xsz;
  static double zz;
  
  //read the problem in canonical form Ax<=b
  void read(istream &in) {
  	
  	in>>m>>n;
  	
  	//coeficientii functiei obiectiv
	for (int i=1; i<=n; ++i) in>>z[i];
	for (int i=1; i<=n; ++i) in>>is_integer[i];
	for (int i=1; i<=m; ++i) in>>semn[i];
	
	for (int i=1; i<=m; ++i) {
	 for (int j=1; j<=n; ++j) in>>a[i][j];
	 in>>bv[i];
    }
  }
  
  static int solve_milp(ILP instance, ostream &cout){	
  	
  	for (int i=1; i<=instance.n; ++i) ILP::x[i]=-1;
  	ILP::xsz=instance.n;
  	ILP::zz=-1000000000;
    
    simplex t;
  	t.build_from_memory(instance.n,instance.m,instance.z,instance.a,instance.bv,instance.semn);
  		
  	int code = t.solve(cout);
    
    while (1) {
    	
      if (code!=1) break; //infeasable or unbounded   	
    	
      int idx=-1;
  	  double max_diff=0;
      for (int i=1; i<=t.m; ++i) 
  	   if (t.a[i][t.n+1]-(int)(t.a[i][t.n+1]+eps)>max_diff) {
  			       max_diff=t.a[i][t.n+1]-(int)t.a[i][t.n+1];
				   idx=i;
  		}
    	
      if (idx==-1 || max_diff<eps) {//solution found
      	       ILP::zz=t.a[t.m+1][t.n+1];
			   for (int i=1; i<=t.n; ++i) ILP::x[i]=0;
			   for (int i=1; i<=t.m; ++i) ILP::x[(int)t.a[i][0]]=t.a[i][t.n+1];
			   
			   break;
      }
      else {//add constraint
      	
        ++t.n;
        ++t.m;
        
        for (int i=1; i<=t.m; ++i) { t.a[i][t.n+1]=t.a[i][t.n]; t.a[i][t.n]=0; }
        for (int i=1; i<=t.n+1; ++i) { t.a[t.m+1][i]=t.a[t.m][i]; t.a[t.m][i]=0; }
      	
      	for (int i=1; i<=t.n+1; ++i) 
      	  if (i!=t.a[idx][0]) {
      	     if (t.a[idx][i]>=0) t.a[t.m][i]=t.a[idx][i]-(int)t.a[idx][i];
			 else t.a[t.m][i]=t.a[idx][i]-(int)t.a[idx][i]+1.0;	
      	  }
      	  
      	t.a[t.m][t.n]=-1;
      	
      	code=t.solve(cout);
      }
    	
    }
  	
  	if (ILP::x[1]!=-1) {
  	cout<<"Integer solution found:\n";
  	for (int i=1; i<=ILP::xsz; ++i) cout<<ILP::x[i]<<" ";
  	cout<<"\nZ="<<ILP::zz<<"\n";
  	return 1;
    }
    else {
	 cout<<"No feasable solution found\n";
	 return -1;
    }
  	
  	
  }	
};

int ILP::xsz=-1;
double ILP::x[NR_VAR]={};
double ILP::zz=-1000000000;

ILP pr;

int main(void) {
	
	//ofstream fout("file.out");
	ifstream cin("output.txt");
	
	
    pr.read(cin);
    ILP::solve_milp(pr,cout);
	
	return 0;
}
