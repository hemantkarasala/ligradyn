// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include "ligra.h"

struct CC_F {
  uintE* IDs, *prevIDs;
  CC_F(uintE* _IDs, uintE* _prevIDs) : 
    IDs(_IDs), prevIDs(_prevIDs) {}
  inline bool update(uintE s, uintE d){ //Update function writes min ID
    uintE origID = IDs[d];
    if(IDs[s] < origID) {
      IDs[d] = min(origID,IDs[s]);
      if(origID == prevIDs[d]) return 1;
    } return 0; }
  inline bool updateAtomic (uintE s, uintE d) { //atomic Update
    uintE origID = IDs[d];
    return (writeMin(&IDs[d],IDs[s]) && origID == prevIDs[d]);
  }
  inline bool cond (uintE d) { return cond_true(d); } //does nothing
};

//function used by vertex map to sync prevIDs with IDs
struct CC_Vertex_F {
  uintE* IDs, *prevIDs;
  CC_Vertex_F(uintE* _IDs, uintE* _prevIDs) :
    IDs(_IDs), prevIDs(_prevIDs) {}
  inline bool operator () (uintE i) {
    prevIDs[i] = IDs[i];
    return 1; }};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long n = GA.n;vertex *V=GA.V;
  uintE* IDs = newA(uintE,n), *prevIDs = newA(uintE,n);
  {parallel_for(long i=0;i<n;i++) IDs[i] = i;} //initialize unique IDs

  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;} 
  vertexSubset Frontier(n,n,frontier); //initial frontier contains all vertices

  while(!Frontier.isEmpty()){ //iterate until IDS converge
    vertexMap(Frontier,CC_Vertex_F(IDs,prevIDs));
    vertexSubset output = edgeMap(GA, Frontier, CC_F(IDs,prevIDs));
    Frontier.del();
    Frontier = output;
  }
  Frontier.del(); free(prevIDs);
  for(long j=0;j<n;j++)cout <<IDs[j]<<endl;

  char* data=NULL;
  char *temp;
  size_t len = 0;
  int edge=0;//edge is 1 when edge is removed.
  int a,b;
  FILE *fp;
  fp =fopen("graphedits","r");



  while(1){
    edge = getline(&data,&len,fp);
    edge = 0;
    temp = strtok(data," ");
    cout << "hi";
    if(temp[0]=='x')break;
    else if(temp[0]=='e')edge=1;
    a=atoi(strtok(NULL," "));
    b=atoi(strtok(NULL,"\n"));
    //else if(temp[0]=='a')a=atoi(strtok(NULL,"\n"));
    startTime();
    if(edge){


      uintE* t = (uintE*)V[a].getOutNeighbors();
      uintE* t2 = (uintE*)malloc(sizeof(uintE)*(V[a].getOutDegree()-1));
      int counter = 0;
      for(long i =0;i<V[a].getOutDegree();i++)if(t[i]!=b)t2[counter++]=t[i];
      V[a].setOutNeighbors(t2);
      V[a].setOutDegree((V[a].getOutDegree()-1));

      t = (uintE*)V[b].getInNeighbors();
      t2 = (uintE*)malloc(sizeof(uintE)*V[b].getInDegree());
      counter = 0;
      for(long i =0;i<V[b].getInDegree();i++)if(t[i]!=a)t2[counter++]=t[i];
      V[b].setInNeighbors(t2);
      V[b].setInDegree((V[b].getInDegree()-1));






      long curr=0;
      long culp = IDs[a];
      for(long l=0;l<n;l++)if(IDs[l]>curr)curr=IDs[l];
      
      long cct = 0;
      bool* frontier2 = newA(bool,n);uintE* prevIDs2 = newA(uintE,n);
      parallel_for(long i=0;i<n;i++){if(IDs[i]==culp){frontier2[i] = 1;cct++;}else frontier2[i]=0;} 
      for(long l=0;l<n;l++)if(IDs[l]==culp)IDs[l]=curr++;
      vertexSubset Frontier2(n,cct,frontier2); //initial frontier contains all vertices
      
      while(!Frontier2.isEmpty()){ //iterate until IDS converge
        vertexMap(Frontier2,CC_Vertex_F(IDs,prevIDs2));
        vertexSubset output = edgeMap(GA, Frontier2, CC_F(IDs,prevIDs2));
        Frontier2.del();
        Frontier2 = output;
      }
      Frontier2.del(); free(prevIDs2);



    }
    else{


      uintE* t = (uintE*)V[a].getOutNeighbors();
      uintE* t2 = (uintE*)malloc(sizeof(uintE)*(V[a].getOutDegree()+1));
      int counter = 0;
      for(long i =0;i<V[a].getOutDegree();i++)t2[counter++]=t[i];
      t2[counter]=b;
      V[a].setOutNeighbors(t2);
      V[a].setOutDegree((V[a].getOutDegree()+1));

      t = (uintE*)V[b].getInNeighbors();
      t2 = (uintE*)malloc(sizeof(uintE)*V[b].getInDegree());
      counter = 0;
      for(long i =0;i<V[b].getInDegree();i++)t2[counter++]=t[i];
      t2[counter]=a;
      V[b].setInNeighbors(t2);
      V[b].setInDegree((V[b].getInDegree()+1));





      long sm;
      long big;
      if(IDs[a]<IDs[b]){sm=IDs[a];big=IDs[b];}
      else {sm=b;big=a;}
      
      if(IDs[a]!=IDs[b])for(long l=0;l<n;l++)if(IDs[l]==big)IDs[l]=sm;



    }
    for(long j=0;j<n;j++)cout <<IDs[j]<<endl;


  }



}
