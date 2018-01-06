#include<fstream>
#include<cstring>
using namespace std;
typedef struct celula {
        int nod;
        celula *next;
}*lista;
lista graf[1005],v;
int n,m,i,flow,cap[1005][1005];
int st,en,tata[1005],coada[1005];
bool viz[1005];
 
void update(){
 
    int minim=1000000000, aux=n;
 
    while (aux!=1) {
        minim=min(minim,cap[tata[aux]][aux]);
        aux=tata[aux];
    }
 
    flow+=minim;
 
    aux=n;
 
    while (aux!=1) {
 
        cap[tata[aux]][aux]-=minim;
        cap[aux][tata[aux]]+=minim;
        aux=tata[aux];
 
    }
 
}
 
bool bfs(){
 
    st=en=1;
    tata[1]=1;
    memset(viz,0,sizeof(viz));
    viz[1]=1;
    coada[st]=1;
    bool ok=0;
 
    while (st<=en) {
        int nodc=coada[st++];
 
        for (lista p=graf[nodc]; p; p=p->next) 
            if (cap[nodc][p->nod]>0&&viz[p->nod]==0) {
 
                viz[p->nod]=1;
                tata[p->nod]=nodc;
                coada[++en]=p->nod;
 
                if (p->nod==n) { update(); ok=1; viz[n]=0; --en; }
 
            }
 
    }
 
    return ok;
}
 
int main(void) {
 
    ifstream cin("maxflow.in");
    ofstream cout("maxflow.out");
 
    cin>>n>>m;
 
    for (i=1; i<=m; ++i) {
        int x,y,z;
        cin>>x>>y>>z;
        v=new celula; v->nod=y; v->next=graf[x]; graf[x]=v; cap[x][y]=z;
        v=new celula; v->nod=x; v->next=graf[y]; graf[y]=v; 
    }
 
    while (bfs()) {};
 
    cout<<flow;
 
    return 0;
}
