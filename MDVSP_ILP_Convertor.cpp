#include<fstream>
#include<cstring>
using namespace std;

/*
2 3 1 1

1 3 6
2 5 7

0 100000000 2 2 2
100000000 0 2 2 2
2 2 100000000 1 1
2 2 100000000 100000000 3
2 2 100000000 100000000 100000000
*/

int n,k,m,st[1000], en[1000], v[10], cost[1000][1000];
int inf=100000000;
int coef[1000];

int main(void) {
   
    ifstream cin("input.txt");
	ofstream cout("output.txt");
	
	cin>>k>>n;
	for (int i=1; i<=k; ++i) cin>>v[i];	
	for (int i=1; i<=n; ++i) cin>>st[i];
	for (int i=1; i<=n; ++i) cin>>en[i];
	
	int nrx=0;
	
	for (int i=1; i<=n+k; ++i)
	 for (int j=1; j<=n+k; ++j) cin>>cost[i][j];
    
    for (int i=1; i<=n; ++i)
     for (int j=1; j<=n; ++j) 
      if (cost[i][j]==inf) ++nrx;
    
    ++nrx;
    cout<<nrx*k<<"\n";
	 
	 //write m and n for ILP solver
	cout<<n+k+(n+1)*k<<" "<<k*(n+1)*(n+1)+k<<"\n";
	
	//write coeficients for objective function
	for (int i=1; i<=k; ++i) {
	   cout<<"100000000 ";
	   for (int j=1; j<=n; ++j) cout<<cost[i][j+k]<<" ";
	   	
	   for (int j=1; j<=n; ++j) {
	   	 cout<<cost[j+k][i]<<" ";
	   	 for (int t=1; t<=n; ++t) cout<<cost[j+k][t+k]<<" ";
	   }
	}
	
	for (int i=1; i<=k; ++i) cout<<"0 ";
	
	cout<<"\n";
	
	//write vector marking integern and non integer variables
	for (int i=1; i<=k*(n+1)*(n+1)+k; ++i) cout<<"1 "; //all variables are integers
	cout<<"\n";
	
	for (int i=1; i<=n+k+(n+1)*k; ++i) cout<<"0 ";
	cout<<"\n";
	//write equations for system Ax=b
	
	//first n equalites ensures that only one edge leaves a vertex in some of those k graphs associated with depots
	for (int i=1; i<=n; ++i) { 
	 for (int j=1; j<=k; ++j) 
	  for (int t=0; t<=n; ++t)
	   for (int s=0; s<=n; ++s) {
	   	
	   	if (t==i && (s==0 || cost[t+k][s+k]!=inf)) cout<<"1 ";
		else cout<<"0 "; 
	   	
	   }
	  
	 for (int j=1; j<=k; ++j) cout<<"0 ";
	 cout<<"1\n";   
    }
	
	//next k equalities ensures that maximum number of vehicles in each depot is respected
	for (int i=1; i<=k; ++i) {
	 
	 for (int j=1; j<=k; ++j) 
	  for (int t=0; t<=n; ++t)
	   for (int s=0; s<=n; ++s) {
	   	
	   	if (j==i && t==0 && s>0) cout<<"1 ";
		else cout<<"0 "; 
	   	
	   }
	 
	 for (int j=1; j<=k; ++j)
	  if (j==i) cout<<"1 "; else cout<<"0 ";
	  
	 cout<<v[i]<<"\n";   
    }
    
    //last (n+1)*k equalities ensures that flow is conserved
    for (int i=0; i<=n; ++i) 
	 for (int g=1; g<=k; ++g) {
	 
	 for (int j=1; j<=k; ++j) 
	  for (int t=0; t<=n; ++t)
	   for (int s=0; s<=n; ++s) {
	   	
	     if (j!=g || t==s) cout<<"0 ";
	     else {//i am in the right graph
	     	
	        if (t==i) {//current edge goes from i to s
			   if (i==0 && cost[j][s+k]!=inf) cout<<"1 ";
			   else if (i>0 && s>0 && cost[i+k][s+k]!=inf) cout<<"1 ";
			   else if (i>0 && s==0 && cost[i+k][j]!=inf) cout<<"1 ";
			   else cout<<"0 ";
		    }
		    else if (s==i) {//current edge goes from t to i
		    	if (i==0 && cost[t+k][j]!=inf) cout<<"-1 ";
		    	else if (i>0 && t>0 && cost[t+k][i+k]!=inf) cout<<"-1 ";
		    	else if (i>0 && t==0 && cost[j][i+k]!=inf) cout<<"-1 ";
		    	else cout<<"0 ";
		    }
		    else cout<<"0 ";
	     	
	     }
	   	
	   }
	 
	 for (int j=1; j<=k; ++j) cout<<"0 ";
	  
	 cout<<"0\n";   
    }
    
	return 0;
}
