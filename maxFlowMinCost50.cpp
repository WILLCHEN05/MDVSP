#include<fstream>
#define inf 1000000000
using namespace std;
typedef struct celula{
        int nod;
        celula *next;
        }*lista;
lista graf[355],v,f;
int i,j,n,m,s,d,tata[355],dist[355],cost[355][355],cap[355][355],sol;
bool ok=true;
 
inline void update(){
       int nod=d,min=1000000;
        while (tata[nod]!=nod){
          if (cap[tata[nod]][nod]<min) min=cap[tata[nod]][nod];
              nod=tata[nod];
           }
        nod=d; sol+=min*dist[d];
        while (tata[nod]!=nod){
            cap[tata[nod]][nod]-=min;
            cap[nod][tata[nod]]+=min;
            nod=tata[nod];
           }
}
 
inline void ford(){
       for (i=1; i<=n; ++i) dist[i]=inf; dist[s]=0; tata[s]=s;
       v=new celula; v->nod=s; v->next=0; f=v; 
       while (v!=0) {
             for (lista p=graf[v->nod];p; p=p->next)
              if ( cap[v->nod][p->nod]>0&&dist[v->nod]+cost[v->nod][p->nod]<dist[p->nod] ){
                   dist[p->nod]=dist[v->nod]+cost[v->nod][p->nod];
                   lista aux=new celula; aux->nod=p->nod; aux->next=0; f->next=aux; f=aux;
                   tata[p->nod]=v->nod;
                  }
         v=v->next;
        }
       if (dist[d]==inf) ok=false;
}
 
int main(void){
    ifstream fin("fmcm.in");
    ofstream fout("fmcm.out");
    fin>>n>>m>>s>>d; int x,y,c,z;
    for (i=1; i<=m; ++i){
        fin>>x>>y>>c>>z; cap[x][y]=c; cost[x][y]=z; cost[y][x]=-z;
        v=new celula; v->nod=y; v->next=graf[x]; graf[x]=v;
        v=new celula; v->nod=x; v->next=graf[y]; graf[y]=v;
    }
    while (ok==true){ford(); if (ok) update(); }
     fout<<sol;
  return(0);
} 
