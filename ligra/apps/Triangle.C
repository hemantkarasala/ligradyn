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

// Triangle counting code (assumes a symmetric graph, so pass the "-s"
// flag). This is not optimized (no ordering heuristic is used)--for
// optimized code, see "Multicore Triangle Computations Without
// Tuning", ICDE 2015. Currently only works with uncompressed graphs,
// and not with compressed graphs.
#include "ligra.h"
#include "quickSort.h"



//assumes sorted neighbor lists
template <class vertex>
long countCommon(vertex& A, vertex& B, uintE a, uintE b) { 
  uintT i=0,j=0,nA = A.getInDegree(), nB = B.getOutDegree();
  uintE* nghA = (uintE*) A.getInNeighbors(), *nghB = (uintE*) B.getOutNeighbors();
  
  long ans=0;

  while (i < nA && j < nB && nghA[i] < a && nghB[j] < b) { //count "directed" triangles
    if (nghA[i]==nghB[j]) i++, j++, ans++;
    else if (nghA[i] < nghB[j]) i++;
    else j++;
  }
  //for(;i<nA;i++)for(;j<nB;j++)if(nghA[i]==nghB[j] && a>nghA[i])ans++;
  
  return ans;
}

template <class vertex>
struct countF { //for edgeMap
  vertex* V;
  long* counts; 
  countF(vertex* _V, long* _counts) : V(_V), counts(_counts) {}
  inline bool update (uintE s, uintE d) {
    if(s > d) //only count "directed" triangles
      writeAdd(&counts[s], countCommon<vertex>(V[s],V[d],s,d));
    return 1;
  }
  inline bool updateAtomic (uintE s, uintE d) {
    if (s > d) //only count "directed" triangles
      writeAdd(&counts[s], countCommon<vertex>(V[s],V[d],s,d));
    return 1;
  }
  inline bool cond (uintE d) { return cond_true(d); } //does nothing
};

struct intLT { bool operator () (uintT a, uintT b) { return a < b; }; };

template <class vertex>
struct initF { //for vertexMap to initial counts and sort neighbors for merging
  vertex* V;
  long* counts;
  initF(vertex* _V, long* _counts) : V(_V), counts(_counts) {}
  inline bool operator () (uintE i) {
    counts[i] = 0;
    quickSort(V[i].getOutNeighbors(),V[i].getOutDegree(),intLT());
    quickSort(V[i].getInNeighbors(),V[i].getInDegree(),intLT());
    return 1;
  }
};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  uintT n = GA.n;
  vertex *V=GA.V;
  char* data=NULL;
  long count2;
  long count3;
  long* counts = newA(long,n);
  bool* frontier = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) frontier[i] = 1;} 
  vertexSubset Frontier(n,n,frontier); //frontier contains all vertices
  startTime();
  vertexMap(Frontier,initF<vertex>(GA.V,counts));
  edgeMap(GA,Frontier,countF<vertex>(GA.V,counts));
  long count = sequence::plusReduce(counts,n);
  for(long i=0;i<n;i++)cout<<counts[i]<<endl;
  cout << "triangle count = " << count << endl;
  Frontier.del(); free(counts);
  nextTime("Running time");

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
    
    if(temp[0]=='x')break;
    else if(temp[0]=='e')edge=1;
    a=atoi(strtok(NULL," "));
    b=atoi(strtok(NULL,"\n"));
    //else if(temp[0]=='a')a=atoi(strtok(NULL,"\n"));
    startTime();
    if(edge){
      


      long* counts2 = newA(long,n);
      for(long g=0;g<n;g++)counts2[g]=0;
      bool* frontier2 = newA(bool,n);
      {parallel_for(long i=0;i<n;i++)frontier2[i]=0;}
      
      uintE* tt1=(uintE*)V[a].getOutNeighbors();
      uintE* tt2=(uintE*)V[a].getInNeighbors();
      uintE* tt3=(uintE*)V[b].getOutNeighbors();
      uintE* tt4=(uintE*)V[b].getInNeighbors();
      long cnt = 2;
      frontier2[a]=1;frontier2[b]=1;
      for(long it=0;it<V[a].getOutDegree();it++)if(frontier2[tt1[it]]!=1){frontier2[tt1[it]]=1;cnt++;}
      for(long it=0;it<V[a].getInDegree();it++)if(frontier2[tt2[it]]!=1){frontier2[tt2[it]]=1;cnt++;}
      for(long it=0;it<V[b].getOutDegree();it++)if(frontier2[tt3[it]]!=1){frontier2[tt3[it]]=1;cnt++;}
      for(long it=0;it<V[b].getInDegree();it++)if(frontier2[tt4[it]]!=1){frontier2[tt4[it]]=1;cnt++;}

      vertexSubset Frontier2(n,cnt,frontier2); //frontier contains all vertices

      vertexMap(Frontier2,initF<vertex>(GA.V,counts2));
      edgeMap(GA,Frontier2,countF<vertex>(GA.V,counts2));
      count2 = 0;
      for(long ii=0;ii<n;ii++){count2+=counts2[ii];cout<<counts2[ii]<<endl;}
      Frontier2.del(); free(counts2);
      


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

      long* counts3 = newA(long,n);
      for(long g=0;g<n;g++)counts3[g]=0;
      bool* frontier3 = newA(bool,n);

      {parallel_for(long i=0;i<n;i++)frontier3[i]=0;}    
      tt1=(uintE*)V[a].getOutNeighbors();
      tt2=(uintE*)V[a].getInNeighbors();
      tt3=(uintE*)V[b].getOutNeighbors();
      tt4=(uintE*)V[b].getInNeighbors();
      cnt = 2;
      frontier3[a]=1;frontier3[b]=1;
      for(long it=0;it<V[a].getOutDegree();it++)if(frontier3[tt1[it]]!=1){frontier3[tt1[it]]=1;cnt++;}
      for(long it=0;it<V[a].getInDegree();it++)if(frontier3[tt2[it]]!=1){frontier3[tt2[it]]=1;cnt++;}
      for(long it=0;it<V[b].getOutDegree();it++)if(frontier3[tt3[it]]!=1){frontier3[tt3[it]]=1;cnt++;}
      for(long it=0;it<V[b].getInDegree();it++)if(frontier3[tt4[it]]!=1){frontier3[tt4[it]]=1;cnt++;}

      vertexSubset Frontier3(n,cnt,frontier3); //frontier contains all vertices

      vertexMap(Frontier3,initF<vertex>(GA.V,counts3));
      edgeMap(GA,Frontier3,countF<vertex>(GA.V,counts3));
      count3 = 0;
      for(long ii=0;ii<n;ii++)count3+=counts3[ii];
      Frontier3.del(); free(counts3);
      //cout << "triangle count2 = " << count2 << endl;
      

      count += count3 - count2;
      cout << "triangle count2 = " << count <<endl;
      cout << count2<<endl;
      cout << count3<<endl;
      
      
    }

    else{

      long* counts3 = newA(long,n);
      for(long g=0;g<n;g++)counts3[g]=0;
      bool* frontier3 = newA(bool,n);
      {parallel_for(long i=0;i<n;i++)frontier3[i]=0;}    
      uintE* tt1=(uintE*)V[a].getOutNeighbors();
      uintE* tt2=(uintE*)V[a].getInNeighbors();
      uintE* tt3=(uintE*)V[b].getOutNeighbors();
      uintE* tt4=(uintE*)V[b].getInNeighbors();
      long cnt = 2;
      frontier3[a]=1;frontier3[b]=1;
      for(long it=0;it<V[a].getOutDegree();it++)if(frontier3[tt1[it]]!=1){frontier3[tt1[it]]=1;cnt++;}
      for(long it=0;it<V[a].getInDegree();it++)if(frontier3[tt2[it]]!=1){frontier3[tt2[it]]=1;cnt++;}
      for(long it=0;it<V[b].getOutDegree();it++)if(frontier3[tt3[it]]!=1){frontier3[tt3[it]]=1;cnt++;}
      for(long it=0;it<V[b].getInDegree();it++)if(frontier3[tt4[it]]!=1){frontier3[tt4[it]]=1;cnt++;}

      vertexSubset Frontier3(n,cnt,frontier3);

      vertexMap(Frontier3,initF<vertex>(GA.V,counts3));
      edgeMap(GA,Frontier3,countF<vertex>(GA.V,counts3));
      count3 = 0;
      for(long ii=0;ii<n;ii++)count3+=counts3[ii];
      
      //cout << "triangle count2 = " << count2 << endl;
      Frontier3.del(); free(counts3);


   
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

      long* counts2 = newA(long,n);
      for(long g=0;g<n;g++)counts2[g]=0;
      bool* frontier2 = newA(bool,n);
      {parallel_for(long i=0;i<n;i++)frontier2[i]=0;}
      
      tt1=(uintE*)V[a].getOutNeighbors();
      tt2=(uintE*)V[a].getInNeighbors();
      tt3=(uintE*)V[b].getOutNeighbors();
      tt4=(uintE*)V[b].getInNeighbors();
      cnt = 2;
      frontier2[a]=1;frontier2[b]=1;
      for(long it=0;it<V[a].getOutDegree();it++)if(frontier2[tt1[it]]!=1){frontier2[tt1[it]]=1;cnt++;}
      for(long it=0;it<V[a].getInDegree();it++)if(frontier2[tt2[it]]!=1){frontier2[tt2[it]]=1;cnt++;}
      for(long it=0;it<V[b].getOutDegree();it++)if(frontier2[tt3[it]]!=1){frontier2[tt3[it]]=1;cnt++;}
      for(long it=0;it<V[b].getInDegree();it++)if(frontier2[tt4[it]]!=1){frontier2[tt4[it]]=1;cnt++;}

      vertexSubset Frontier2(n,cnt,frontier2); //frontier contains all vertices

      vertexMap(Frontier2,initF<vertex>(GA.V,counts2));
      edgeMap(GA,Frontier2,countF<vertex>(GA.V,counts2));
      count2 = 0;
      for(long ii=0;ii<n;ii++)count2+=counts2[ii];
      Frontier2.del(); free(counts2);
      count +=count2-count3;
      cout << "triangle count2 = " << count << endl;

      


    }
    nextTime("Running time");


  }


}
