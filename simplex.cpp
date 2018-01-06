#include<iostream>
#include<fstream>
#include<ostream>
#include<istream>
#include<cstring>
using namespace std;

const double eps=0.0000000001;

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

struct simplex {
  int n,m;
  int n1, m1;
  double b[200][200];
  double a[200][200];
  double z[200];	
  bool has_solution=1;
  
  //read the problem in standard form (Ax=b)
  void read(istream &in) {
  	
  	in>>m>>n;
  	
  	//coeficientii functiei obiectiv
	for (int i=1; i<=n; ++i) in>>z[i];
	
	for (int i=1; i<=m; ++i)
	 for (int j=1; j<=n+1; ++j) in>>a[i][j];
	 
	for (int i=1; i<=n; ++i) a[m+1][i]=z[i];
  }
  
  void printState(ostream &out, double (*a)[200], int n, int m) {
	out<<"\n";
  	for (int i=1; i<=m+1; ++i) {
  		out<<a[i][0]<<" | ";
  		for (int j=1; j<=n; ++j) out<<a[i][j]<<" ";
  		out<<"| "<<a[i][n+1]<<"\n";
  	}
	out<<"\n";
  }
  
  void solve_simplex_table(ostream &out, double (*a)[200], int n, int m) {
  	
  	out<<"Enter solve simplex table:\n";
  	int it=0;
	while (1) {
	  
	  //check if the solution is optim
	  bool optim=1;
	  int col;
	  for (int i=1; i<=n; ++i)
	   if (low(a[m+1][i],0)) { optim=0; col=i; break; }	
	  	
	  if (optim)  { out<<"Reached final state:\n"; break; }
	  
	  out<<"Iteration #"<<it++<<"\n";
	  printState(out, a, n, m);
	  
	  double min=1000000000;
	  int lin=-1;
	  
	  for (int i=1; i<=m; ++i)
	   if ( gt(a[i][col],0) && loweq(a[i][n+1]/a[i][col],min) ){
	   	  if ( low(a[i][n+1]/a[i][col],min) ) { lin=i; min=a[i][n+1]/a[i][col]; }
	   	  else if ( (int)a[i][0]<(int)a[lin][0] ) lin=i;
	   }
	   
	  if (lin==-1) { out<<"Infinite Optimum:\n"; break; }
	  
	  out<<"Pivot: "<<lin<<" "<<col<<"\n";
	  
	  //pivot location is (lin,col)
	  for (int i=1; i<=m+1; ++i)
	   for (int j=1; j<=n+1; ++j)
	    if (i!=lin && j!=col) a[i][j]=a[i][j]-a[lin][j]*a[i][col]/a[lin][col];
	  
	  for (int i=1; i<=n+1; ++i) if (i!=col) a[lin][i]/=a[lin][col];
	  for (int i=1; i<=m+1; ++i) a[i][col]=0;
	  
	  a[lin][col]=1;
	  
	  a[lin][0]=col;
		
	}
	
	printState(out, a, n, m);	
  	
  }
  
  void normalize_base(double (*a)[200], int n, int m) {
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
  
  void solve(ostream &out) {
  	
    bool ok = buildBase();
	
	if (ok) solve_simplex_table(out, a, n, m);
	else {
	  ok = solve_phase1(out);
	  if (ok) solve_simplex_table(out, a, n, m);
	  else out<<"No feasable solutions exists\n";
    }
    
    if (has_solution) {
    //check if this is the only solution	
     bool only=1;
	 for (int i=1; i<=m; ++i)
	  if (a[i][n+1]==0) { only=0; break; }	
    	
     for (int i=1; i<=n; ++i) 
	  if (a[m+1][i]==0) {
	  	 int cnt=0;
	     for (int j=1; j<=m; ++j)
	       if (a[j][0]==i) { cnt=1; break; }
	       
	     if (cnt==0) { only=0; break; }  
      }
    	
      if (only) out<<"This is the unique optimal solution\n";
      else out<<"Infinite number of optimal solutions\n";
    }
  	
  }
};

simplex t;

int main(void) {
	
    t.read(cin);
    t.solve(cout);
	
	return 0;
}
