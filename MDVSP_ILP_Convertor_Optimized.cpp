#include<fstream>
#include<cstring>
using namespace std;

int n,k,m,st[1000], en[1000], v[10], cost[1000][1000];
int inf=100000000, lim=1399;
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
      if (cost[i+k][j+k]>lim) ++nrx;
	 
	 //write m and n for ILP solver
	cout<<n+k+(n+1)*k<<" "<<k*(n+1)*(n+1)-nrx*k<<"\n";
	
	int ver=k*(n+1)*(n+1)-nrx*k;
	
	int auxver=0;
	//write coeficients for objective function
	for (int i=1; i<=k; ++i) {

	   for (int j=1; j<=n; ++j) cout<<cost[i][j+k]<<" ";
	   auxver+=n;
	   	
	   for (int j=1; j<=n; ++j) {
	   	 cout<<cost[j+k][i]<<" ";
	   	 ++auxver;
	   	 
	   	 for (int t=1; t<=n; ++t) 
			if (cost[j+k][t+k]<=lim) { cout<<cost[j+k][t+k]<<" "; ++auxver; }
	   }
	}
	
	for (int i=1; i<=k; ++i) cout<<"0 ";
	auxver+=k;
	cout<<"\n";
	
//	cout<<ver<<" "<<auxver<<"\n";
	auxver=0;
	
	//write vector marking integer and non integer variables
	for (int i=1; i<=k*(n+1)*(n+1)-nrx*k; ++i) cout<<"1 "; //all variables are integers
	cout<<"\n";
	
	//write inequalities signs (0 means =)
	for (int i=1; i<=n+k+(n+1)*k; ++i) cout<<"0 ";
	cout<<"\n";
	//write equations for system Ax=b
	
	//first n equalites ensures that only one edge leaves a vertex in some of those k graphs associated with depots
	for (int i=1; i<=n; ++i) { 
	 auxver=0;
	 for (int j=1; j<=k; ++j) 
	  for (int t=0; t<=n; ++t)
	   for (int s=0; s<=n; ++s) 
	    if (s==0 && t==0 ) continue;
	    else if (s==0 || t==0 || cost[t+k][s+k]<=lim) {
	   	
	   	if (t==i) cout<<"1 ";
		else cout<<"0 "; 
	   	
	   	++auxver;
	   }
	  
	 for (int j=1; j<=k; ++j) cout<<"0 ";
	 auxver+=k;
	 cout<<"1\n";
	 
	// cout<<auxver<<"\n";   
    }
	
	//next k equalities ensures that maximum number of vehicles in each depot is respected
	for (int i=1; i<=k; ++i) {
	 auxver=0;
	 for (int j=1; j<=k; ++j) 
	  for (int t=0; t<=n; ++t)
	   for (int s=0; s<=n; ++s)
	    if ( s==0 && t==0 ) continue;
	    else if (s==0 || t==0 || cost[t+k][s+k]<=lim) {
	   	
	   	if (j==i && t==0 && s>0) cout<<"1 ";
		else cout<<"0 "; 
	   	++auxver;
	   }
	 
	 for (int j=1; j<=k; ++j)
	  if (j==i) cout<<"1 "; else cout<<"0 ";
	  
	 cout<<v[i]<<"\n";   
	 auxver+=k;
	 //cout<<auxver<<"\n";
    }
    
    //last (n+1)*k equalities ensures that flow is conserved
    for (int i=0; i<=n; ++i) //for each node
	 for (int g=1; g<=k; ++g) { //for each graph
	 auxver=0; 
	 for (int j=1; j<=k; ++j) 
	  for (int t=0; t<=n; ++t)
	   for (int s=0; s<=n; ++s) 
	    if ( s==0 && t==0 ) continue;
	    else if (s==0 || t==0 || cost[t+k][s+k]<=lim) {
	   	
	     if (j!=g || t==s) cout<<"0 ";
	     else {//i am in the right graph
	     	
	        if (t==i) {//current edge goes from i to s, mark with 1
			   if (i==0) cout<<"1 ";
			   else if (i>0 && s>=0) cout<<"1 ";
			   else cout<<"0 ";
		    }
		    else if (s==i) {//current edge goes from t to i, mark with -1
		    	if (i==0) cout<<"-1 ";
		    	else if (i>0 && t>=0) cout<<"-1 ";
		    	else cout<<"0 ";
		    }
		    else cout<<"0 ";
	     	
	     }
	   	++auxver;
	   }
	 
	 for (int j=1; j<=k; ++j) cout<<"0 ";
	  
	 cout<<"0\n";   
	 auxver+=k;
	 //cout<<auxver<<"\n";
    }
    
	return 0;
}
