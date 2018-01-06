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
const int NR_VAR=25000;

struct simplex {
  int n,m;
  int n1, m1, semn[NR_VAR];
  double b[NR_LINES][NR_VAR];
  double a[NR_LINES][NR_VAR];
  double bv[NR_LINES];
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
	   if ( i==lin || a[i][col]==0 ) continue;
	   else for (int j=1; j<=n+1; ++j)
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


struct MILP {
  int n,m,semn[NR_VAR];
  double a[NR_LINES][NR_VAR];
  double bv[NR_LINES];
  double z[NR_VAR];	
  bool is_integer[NR_VAR];
  bool has_solution=1;
  static stack<MILP> st;
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
  
  static int solve_milp(MILP instance, ostream &cout){
  	
    while (!MILP::st.empty()) MILP::st.pop();	
    MILP::st.push(instance); 	
  	
  	for (int i=1; i<=instance.n; ++i) MILP::x[i]=-1;
  	MILP::xsz=instance.n;
  	MILP::zz=-1000000000;
  	
  	while (!MILP::st.empty()) {
  		
  		MILP current = MILP::st.top();
  		MILP::st.pop();
  		
  		cout<<"Solving problem:\n";
  		
  	/*	for (int i=1; i<=current.m; ++i) {
  		   for (int j=1; j<=current.n; ++j) cout<<current.a[i][j]<<" ";	
  		   if (current.bv[i]>0) cout<<" = "<<current.bv[i]<<"\n";
  		   else cout<<" >= "<<-current.bv[i]<<"\n";
  		}
  	*/	
  		//solve using two phase method
  		simplex t;
  		t.build_from_memory(current.n,current.m,current.z,current.a,current.bv,current.semn);
  		
  		int code = t.solve(cout);
  		
        if (code!=1) continue; //infeasable or unbounded		
  		
  		//get current optimal solution
  		double zc=t.a[t.m+1][t.n+1];
  		cout<<"found z="<<zc<<"\n";
  		
  		if (zc>MILP::zz) { //possibility to branch
  			
  			//check integrality
  			bool integrality=1;
  			for (int i=1; i<=t.m; ++i) 
  			 if (current.is_integer[(int)t.a[i][0]] && t.a[i][t.n+1]-(int)(t.a[i][t.n+1]+eps)>eps) { integrality=0; break; }
  			 
  			if (integrality) {
  			   MILP::zz=zc;
  			   cout<<"best found:"<<zz<<"\n";
			   for (int i=1; i<=t.n; ++i) MILP::x[i]=0;
			   for (int i=1; i<=t.m; ++i) MILP::x[(int)t.a[i][0]]=t.a[i][t.n+1];
  			}
  			else { //split by non-integrality
  			
  			    int idx;
  			    int val;
  			    double max_diff=0;
  			    for (int i=1; i<=t.m; ++i) 
  			     if (current.is_integer[(int)t.a[i][0]] && t.a[i][t.n+1]-(int)(t.a[i][t.n+1]+eps)>max_diff) {
  			       max_diff=t.a[i][t.n+1]-(int)t.a[i][t.n+1];
				   idx=t.a[i][0];	
				   val=(int)t.a[i][t.n+1];
  			     }
  			     
  			     cout<<idx<<" "<<max_diff<<"\n";
				   
				   MILP lower_bound = current;
				   ++lower_bound.m;
				   for (int i=1; i<=lower_bound.n; ++i) lower_bound.a[lower_bound.m][i]=0;
				   lower_bound.a[lower_bound.m][idx]=1;
				   lower_bound.bv[lower_bound.m]=val;
				   MILP::st.push(lower_bound);
				   
				   MILP upper_bound = current;
				   ++upper_bound.m;
				   for (int i=1; i<=upper_bound.n; ++i) upper_bound.a[upper_bound.m][i]=0;
				   upper_bound.a[upper_bound.m][idx]=1;
				   upper_bound.bv[upper_bound.m]=val+1;
				   MILP::st.push(upper_bound);
				   /*MILP upper_bound = current;
				   ++upper_bound.m;
				   for (int i=1; i<=upper_bound.n; ++i) upper_bound.a[upper_bound.m][i]=0;
				   upper_bound.a[upper_bound.m][idx]=1;
				   upper_bound.bv[upper_bound.m]=-(val+1);
				   MILP::st.push(upper_bound);
				   */
  			}
  			
  		}
  	}
  	
  	if (MILP::x[1]!=-1) {
  	cout<<"Integer solution found:\n";
  	for (int i=1; i<=MILP::xsz; ++i) cout<<MILP::x[i]<<" ";
  	cout<<"\nZ="<<MILP::zz<<"\n";
  	return 1;
    }
    else {
	 cout<<"No feasable solution found\n";
	 return -1;
    }
  	
  }	
};

int MILP::xsz=-1;
double MILP::x[NR_VAR]={};
stack<MILP> MILP::st;
double MILP::zz=-1000000000;

MILP pr;

int main(void) {
	
//	ofstream fout("file.out");
    ifstream cin("output.txt");

    pr.read(cin);
    cout<<"read complete\n";
    MILP::solve_milp(pr,cout);
	
	return 0;
}
