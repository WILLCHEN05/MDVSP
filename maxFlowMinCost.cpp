#include<fstream>
using namespace std;
const int inf=100000000;
typedef struct celula{
              int nod;
              celula *next;
              }*lista;
typedef struct {int nod,dist; } tip;
tip heap[10000];
lista graf[355],v;
int i,j,n,m,cap[355][355],cost[355][355],s,d,dist[355],coada[130000],aux[355],tata[355],sol,hp,newdist[355],st,en;
bool viz[355];
  
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
  
int main(void){
    ifstream fin("fmcm.in");
    ofstream fout("fmcm.out");
    fin>>n>>m>>s>>d; int x,y;
    for (i=1; i<=m; ++i) {
          fin>>x>>y; fin>>cap[x][y]>>cost[x][y]; cost[y][x]=-cost[x][y];
          v=new celula; v->nod=x; v->next=graf[y]; graf[y]=v;
          v=new celula; v->nod=y; v->next=graf[x]; graf[x]=v;
          }
   //aflu distanta minima cu Ford si modific costurile muchiilor astfel incit sa nu am arce de cost negativ
    st=1; en=1; coada[en]=s;
   for (i=1; i<=n; ++i) dist[i]=inf, viz[i]=0; dist[s]=0; viz[s]=1;
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
         for (i=1; i<=n; ++i) dist[i]=aux[i];
         while (q!=s) { mmin=min(mmin,cap[tata[q]][q]); q=tata[q]; }
         sol+=mmin*dist[d];
         q=d;
         while (q!=s) { cap[tata[q]][q]-=mmin; cap[q][tata[q]]+=mmin; q=tata[q]; }
         }
   fout<<sol;
  return(0);
}
