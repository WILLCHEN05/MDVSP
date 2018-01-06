#include<iostream>
#include<fstream>
#include<cstring>
#include<vector>
using namespace std;

struct path {
    int s, e, c, ci;	
};

vector<path> paths;
int n,k,m,st[1000], en[1000], vv[10], cost[1000][1000];
int inf=100000000;

//max flow helper
typedef struct celula{
              int nod;
              celula *next;
              }*lista;

typedef struct {int nod,dist; } tip;

lista v;

struct maxFlowMinCost {
	
  tip heap[100000];
  lista graf[505];
  int n,cap[505][505],cost[505][505],s,d,dist[505],coada[1300000],aux[505],tata[505],sol,hp,newdist[505],st,en,flow;
  bool viz[505];
  
  maxFlowMinCost() {
    zero();
  }
  
  void zero() {
  	 memset(heap,0,sizeof(heap));
	 memset(graf,0,sizeof(graf));
	 n=m=0;
	 memset(cap,0,sizeof(cap));
	 memset(cost,0,sizeof(cost));
	 memset(dist,0,sizeof(dist));
	 memset(aux,0,sizeof(aux));
	 memset(tata,0,sizeof(tata));
	 sol=0;
	 hp=0;
	 memset(newdist,0,sizeof(newdist));
	 st=en=0;
	 memset(viz,0,sizeof(viz));	
  }
  
  void swap(int x, int y) { tip w=heap[x]; heap[x]=heap[y]; heap[y]=w; }
  
  void upheap(int nod) { 
     while (nod>1&&heap[nod].dist<=heap[nod/2].dist) { swap(nod,nod/2); nod/=2; }
  }
  
  void downheap(int nod) {
     bool ok=1;
     while (2*nod+1<=hp&&ok){ 
           ok=0;
           int mmin=min(heap[nod].dist,min(heap[2*nod].dist,heap[2*nod+1].dist));
           if (mmin==heap[2*nod].dist) {ok=1; swap(nod,2*nod); nod*=2; }
            else if (mmin==heap[2*nod+1].dist) { ok=1; swap(nod,2*nod+1); nod*=2; ++nod; }
            }
     if (2*nod==hp&&heap[2*nod].dist<heap[nod].dist) swap(nod,2*nod); 
  }
  
  void sterge(){
     heap[1]=heap[hp]; --hp; downheap(1);
  }
  
  void baga(int nod, int dist) {
     ++hp; heap[hp].nod=nod; heap[hp].dist=dist; upheap(hp);
  }
  
  bool exista_drum(){
     for (int i=1; i<=n; ++i) viz[i]=0, newdist[i]=inf;
       newdist[s]=aux[s]=0;
      hp=1; heap[hp].nod=s; heap[hp].dist=0;
     while (hp) {
           int nod=heap[1].nod; sterge();
           if (viz[nod]==0) {
                            viz[nod]=1;
                             for (lista p=graf[nod]; p; p=p->next)
                              if (cap[nod][p->nod]>0) { 
                              int v=cost[nod][p->nod]+dist[nod]-dist[p->nod];
                                
                              if (newdist[p->nod]>v+newdist[nod]) {
                                    newdist[p->nod]=v+newdist[nod];
                                    tata[p->nod]=nod; baga(p->nod,newdist[p->nod]);
                                     aux[p->nod]=aux[nod]+cost[nod][p->nod];
                                     }
                                }
                              }
            }
     return(viz[d]);
  }
  
  void solve() {
  
   //aflu distanta minima cu Ford si modific costurile muchiilor astfel incit sa nu am arce de cost negativ
   st=1; en=1; coada[en]=s;
   for (int i=1; i<=n; ++i) dist[i]=inf, viz[i]=0; dist[s]=0; viz[s]=1;
   
   while (st<=en) {
          int nod=coada[st];
          for (lista p=graf[nod]; p; p=p->next) 
           if (dist[p->nod]>dist[nod]+cost[nod][p->nod]&&cap[nod][p->nod]>0) {
                                                         dist[p->nod]=dist[nod]+cost[nod][p->nod];
                                                         if (viz[p->nod]==0) { ++en; coada[en]=p->nod; viz[p->nod]=1; }
                                                         }
          viz[nod]=0; ++st;
          }
          
   //caut drumuri de ameliorare de cost minim cu Dijkstra
   while (exista_drum()){
         int q=d,mmin=inf;
         for (int i=1; i<=n; ++i) dist[i]=aux[i];
         while (q!=s) { mmin=min(mmin,cap[tata[q]][q]); q=tata[q]; }
         sol+=mmin*dist[d];
         flow+=mmin;
         q=d;
         while (q!=s) { cap[tata[q]][q]-=mmin; cap[q][tata[q]]+=mmin; q=tata[q]; }
         }
  }
	
};

maxFlowMinCost mfmc;
maxFlowMinCost secondFlow;
maxFlowMinCost mcmc;

int main(void) {
	ifstream cin("input.txt");
	
	cin>>k>>n;
	for (int i=1; i<=k; ++i) cin>>vv[i];	
	for (int i=1; i<=n; ++i) cin>>st[i];
	for (int i=1; i<=n; ++i) cin>>en[i];
	
	for (int i=1; i<=n+k; ++i)
	 for (int j=1; j<=n+k; ++j) cin>>cost[i][j];
	 
    //prepare maxflow graph
	//add edges from in to each node
	
	for (int i=1; i<=n; ++i){
	    int mindist=inf;
	    
	    for (int j=1; j<=k; ++j) mindist=min(mindist,cost[j][i+k]);
	    
	    v=new celula(); v->nod=i+1; v->next=mfmc.graf[1]; mfmc.graf[1]=v;
	    v=new celula(); v->nod=1; v->next=mfmc.graf[i+1]; mfmc.graf[i+1]=v;
	    
	    mfmc.cap[1][i+1]=1; mfmc.cap[i+1][1]=0;
	    //mfmc.cost[1][i+1]=mindist; mfmc.cost[i+1][1]=-mindist;
	    mfmc.cost[1][i+1]=0; mfmc.cost[i+1][1]=0;
    }
    
    //add edges from each node to out
	for (int i=1; i<=n; ++i){
	    int mindist=inf;
	    
	    for (int j=1; j<=k; ++j) mindist=min(mindist,cost[i+k][j]);
	    
	    v=new celula(); v->nod=i+1+n; v->next=mfmc.graf[2*n+2]; mfmc.graf[2*n+2]=v;
	    v=new celula(); v->nod=2*n+2; v->next=mfmc.graf[n+i+1]; mfmc.graf[n+i+1]=v;
	    
	    mfmc.cap[n+i+1][2*n+2]=1; mfmc.cap[2*n+2][n+i+1]=0;
	    //mfmc.cost[n+i+1][2*n+2]=mindist; mfmc.cost[2*n+2][n+i+1]=-mindist;
	    mfmc.cost[n+i+1][2*n+2]=0; mfmc.cost[2*n+2][n+i+1]=0;
    }    
    
    //add interior edges
    for (int i=1; i<=n; ++i)
     for (int j=1; j<=n; ++j)
      if (cost[i+k][j+k]!=inf) {
	    v=new celula(); v->nod=i+1; v->next=mfmc.graf[n+j+1]; mfmc.graf[n+j+1]=v;
	    v=new celula(); v->nod=n+j+1; v->next=mfmc.graf[i+1]; mfmc.graf[i+1]=v;         	
      	
	    mfmc.cap[i+1][n+j+1]=1; mfmc.cap[n+j+1][i+1]=0;
	    mfmc.cost[i+1][n+j+1]=cost[i+k][j+k]; mfmc.cost[n+j+1][i+1]=-cost[i+k][j+k];		      	
      }
	
	mfmc.d=2*n+2;
	mfmc.s=1;
	mfmc.n=2*n+2;
	
	mfmc.solve();
	
	//cout<<mfmc.sol<<"\n";
	
	//get found internal paths
	//get starting nodes
	for (int i=2; i<=n+1; ++i) 
	 if (mfmc.cap[i+n][2*n+2]==1) {
	 	path auxp;
	    int mindist=inf;
	    
	    for (int j=1; j<=k; ++j) mindist=min(mindist,cost[j][i+k-1]);
	 	
	 	auxp.c=mindist;
	 	auxp.ci=0;
	 	auxp.s=i;
	 	auxp.e=i;
	 	paths.push_back(auxp);
	 }
	 
	for (int i=0; i<paths.size(); ++i) {
		int next=paths[i].e;
		
		while (1) {
		  for (lista p=mfmc.graf[paths[i].e]; p; p=p->next)
		   if (mfmc.cap[paths[i].e][p->nod]==0) { next=p->nod; break; }
			
		   if (next==1) break;
		   
		   paths[i].c+=mfmc.cost[paths[i].e][next];
		   paths[i].ci+=mfmc.cost[paths[i].e][next];
		   paths[i].e=next-n;
		}
		
	    int mindist=inf;
	    
	    for (int j=1; j<=k; ++j) mindist=min(mindist,cost[paths[i].e+k-1][j]);
		
		paths[i].e+=n;
	    paths[i].c+=mindist;
	}
	
	cout<<"Paths: "<<paths.size()<<"\n";
	int total_cost=0;
	
	for (int i=0; i<paths.size(); ++i) {
	   cout<<paths[i].s<<" "<<paths[i].e<<" "<<paths[i].c<<"\n";
	   total_cost+=paths[i].c;
    }
	
	cout<<"Total Cost (Lower Bound): "<<total_cost<<"\n";
	//after the number of paths has been determined we need to solve the network flow model
	//building the second network
	
	int nrv=paths.size();
	
    for (int i=1; i<=n; ++i){
	    int mindist=inf;
	    
	    for (int j=1; j<=k; ++j) mindist=min(mindist,cost[j][i+k]);
	    
	    v=new celula(); v->nod=i+2; v->next=secondFlow.graf[2]; secondFlow.graf[2]=v;
	    v=new celula(); v->nod=2; v->next=secondFlow.graf[i+2]; secondFlow.graf[i+2]=v;
	    
	    secondFlow.cap[2][i+2]=1; secondFlow.cap[i+2][2]=0;
	    secondFlow.cost[2][i+2]=mindist; secondFlow.cost[i+2][2]=-mindist;
	    //secondFlow.cost[2][i+2]=0; secondFlow.cost[i+2][2]=0;
    }
    
    //add edges from each node to out
	for (int i=1; i<=n; ++i){
	    int mindist=inf;
	    
	    for (int j=1; j<=k; ++j) mindist=min(mindist,cost[i+k][j]);
	    
	    v=new celula(); v->nod=i+2+n; v->next=secondFlow.graf[2*n+3]; secondFlow.graf[2*n+3]=v;
	    v=new celula(); v->nod=2*n+3; v->next=secondFlow.graf[n+i+2]; secondFlow.graf[n+i+2]=v;
	    
	    secondFlow.cap[n+i+2][2*n+3]=1; secondFlow.cap[2*n+3][n+i+2]=0;
	    secondFlow.cost[n+i+2][2*n+3]=mindist; secondFlow.cost[2*n+3][n+i+2]=-mindist;
	    //secondFlow.cost[n+i+2][2*n+3]=0; secondFlow.cost[2*n+3][n+i+2]=0;
    }
    
    //add interior edges
    for (int i=1; i<=n; ++i)
     for (int j=1; j<=n; ++j)
      if (cost[i+k][j+k]!=inf) {
	    v=new celula(); v->nod=i+n+2; v->next=secondFlow.graf[j+2]; secondFlow.graf[j+2]=v;
	    v=new celula(); v->nod=j+2; v->next=secondFlow.graf[i+n+2]; secondFlow.graf[i+n+2]=v;         	
      	
	    secondFlow.cap[i+n+2][j+2]=1; secondFlow.cap[j+2][i+n+2]=0;
	    secondFlow.cost[i+n+2][j+2]=cost[i+k][j+k]; secondFlow.cost[j+2][n+i+2]=-cost[i+k][j+k];		      	
      }    
	
	//add new start and end nodes
	v=new celula(); v->nod=1; v->next=secondFlow.graf[2]; secondFlow.graf[2]=v;
	v=new celula(); v->nod=2; v->next=secondFlow.graf[1]; secondFlow.graf[1]=v;
	
	secondFlow.cap[1][2]=nrv; secondFlow.cap[2][1]=0;
	secondFlow.cost[1][2]=0; secondFlow.cost[2][1]=0; //can be replaced with vehicle cost
	
	v=new celula(); v->nod=2*n+3; v->next=secondFlow.graf[2*n+4]; secondFlow.graf[2*n+4]=v;
	v=new celula(); v->nod=2*n+4; v->next=secondFlow.graf[2*n+3]; secondFlow.graf[2*n+3]=v;
	
	secondFlow.cap[2*n+3][2*n+4]=nrv; secondFlow.cap[2*n+4][2*n+3]=0;
	secondFlow.cost[2*n+3][2*n+4]=0; secondFlow.cost[2*n+4][2*n+3]=0;
	
	for (int i=1; i<=n; ++i) {
	    v=new celula(); v->nod=1; v->next=secondFlow.graf[n+i+2]; secondFlow.graf[n+i+2]=v;
	    v=new celula(); v->nod=n+i+2; v->next=secondFlow.graf[1]; secondFlow.graf[1]=v;         	
      	
	    secondFlow.cap[1][n+i+2]=1; secondFlow.cap[n+i+2][1]=0;
	    secondFlow.cost[1][n+i+2]=0; secondFlow.cost[n+i+2][1]=0;

	    v=new celula(); v->nod=2*n+4; v->next=secondFlow.graf[i+2]; secondFlow.graf[i+2]=v;
	    v=new celula(); v->nod=i+2; v->next=secondFlow.graf[2*n+4]; secondFlow.graf[2*n+4]=v;         	
      	
	    secondFlow.cap[i+2][2*n+4]=1; secondFlow.cap[2*n+4][i+2]=0;
	    secondFlow.cost[i+2][2*n+4]=0; secondFlow.cost[2*n+4][i+2]=0;		 
    }
    
    secondFlow.s=1;
    secondFlow.d=2*n+4;
    secondFlow.n=2*n+4;
    
    secondFlow.solve();
    
    cout<<"Second Flow Value: "<<secondFlow.sol<<"\n";
    
    //get paths from secondFlow
    paths.clear();
    
	for (int i=1; i<=n; ++i) 
	 if (secondFlow.cap[2][i+2]==0) {
	 	path auxp;
	    int mindist=inf;
	    
	    for (int j=1; j<=k; ++j) mindist=min(mindist,cost[j][i+k]);
	 	
	 	auxp.c=mindist;
	 	auxp.ci=0;
	 	auxp.s=i+2;
	 	auxp.e=i+2;
	 	paths.push_back(auxp);
	 }
	 
	for (int i=0; i<paths.size(); ++i) {
		int next=paths[i].e;
		
		while (1) {
		  for (lista p=secondFlow.graf[paths[i].e+n]; p; p=p->next)
		   if (secondFlow.cap[paths[i].e+n][p->nod]==0) { next=p->nod; break; }
			
		   if (next==2*n+3) break;
		   
		   paths[i].c+=secondFlow.cost[paths[i].e+n][next];
		   paths[i].ci+=secondFlow.cost[paths[i].e+n][next];
		   paths[i].e=next;
		}
		
	    int mindist=inf;
	    
	    for (int j=1; j<=k; ++j) mindist=min(mindist,cost[paths[i].e+k-2][j]);
		
		paths[i].e+=n;
	    paths[i].c+=mindist;
	}
	
	cout<<"Paths2: "<<paths.size()<<"\n";
    total_cost=0;
	
	for (int i=0; i<paths.size(); ++i) {
	   cout<<paths[i].s<<" "<<paths[i].e<<" "<<paths[i].c<<"\n";
	   total_cost+=paths[i].c;
    }
	
	cout<<"Total Cost (Lower Bound)2: "<<total_cost<<"\n";
	
	//now each path must be transformed in a single node and solve the problem of assigning nodes to depots according to depot capacity
	//add edges from source to depots
	for (int i=1; i<=k; ++i) {
	    v=new celula(); v->nod=i+1; v->next=mcmc.graf[1]; mcmc.graf[1]=v;
	    v=new celula(); v->nod=1; v->next=mcmc.graf[i+1]; mcmc.graf[i+1]=v;
		
	    mcmc.cap[1][i+1]=vv[i]; mcmc.cap[i+1][1]=0;
	    mcmc.cost[1][i+1]=0; mcmc.cost[i+1][1]=0;		
	}
	
	//add edges from paths to sink
	for (int i=k+1; i<=k+1+paths.size(); ++i) {
	    v=new celula(); v->nod=i; v->next=mcmc.graf[k+paths.size()+2]; mcmc.graf[k+paths.size()+2]=v;
	    v=new celula(); v->nod=k+paths.size()+2; v->next=mcmc.graf[i]; mcmc.graf[i]=v;
		
	    mcmc.cap[i][k+paths.size()+2]=1; mcmc.cap[k+paths.size()+2][i]=0;
	    mcmc.cost[k+paths.size()+2][i]=0; mcmc.cost[i][k+paths.size()+2]=0;		
	}	
	
	//add interior edges
	for (int i=1; i<=k; ++i)
	 for (int j=1; j<=paths.size(); ++j) {
	    v=new celula(); v->nod=i+1; v->next=mcmc.graf[k+j+1]; mcmc.graf[k+j+1]=v;
	    v=new celula(); v->nod=k+j+1; v->next=mcmc.graf[i+1]; mcmc.graf[i+1]=v;	 	
	 	
	    mcmc.cap[i+1][k+j+1]=1; mcmc.cap[k+j+1][i+1]=0;
	   // mcmc.cost[i+1][k+j+1]=cost[i][paths[j-1].s-1+k]+cost[paths[j-1].e-n-1+k][i]; mcmc.cost[k+j+1][i+1]=-mcmc.cost[i+1][k+j+1]; //without first flow part
		mcmc.cost[i+1][k+j+1]=cost[i][paths[j-1].s-2+k]+cost[paths[j-1].e-n-2+k][i]; mcmc.cost[k+j+1][i+1]=-mcmc.cost[i+1][k+j+1];		 	
	 }

	mcmc.d=k+paths.size()+2;
	mcmc.s=1;
	mcmc.n=k+paths.size()+2;	 
	
	mcmc.solve();
	
	int found_solution=mcmc.sol;
	for (int i=0; i<paths.size(); ++i)
	 found_solution+=paths[i].ci;
	 
	cout<<"Found Integer Solution (Upper Bound): "<<found_solution<<"\n";
	cout<<"Optimality Gap: "<<found_solution-total_cost<<"\n";
	
	return 0;
}
