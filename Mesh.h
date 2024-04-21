// manifold and watertight mesh

#ifndef MESH_H
#define MESH_H

#include "VertexEdgeFace.h"
#include <vector>
#include <string>

using std::ofstream;
using std::vector;
using std::string;

class way
{
public:
    vector<int> tree;
    unsigned int min;
	
    way()
    {
		min=0;
    }
};


class Mesh
{

public:

  vector<Vertex> vertex;
  vector<Edge> edge;
  vector<Face> face;

  vector<Face> faceT;
  vector<Edge> edgeT;
  vector<Vertex> vertexT;
  
   vector<int> vertexR;
	
   bool target;
	
  bool reverse;
  vector<int> comb;
  int gain;
  int moy;
	
  int interval;

  int nbflip;

  Mesh();
  ~Mesh();
  Mesh(const string &filename);

  int vertexAdd(float x,float y,float z);
  int vertexAddT(float x,float y,float z);
	
  int faceAdd(int v0,int v1,int v2,int e0,int e1,int e2);
  int faceAddT(int v0,int v1,int v2,int e0,int e1,int e2);
	
  int edgeAdd(int v0,int v1,int f0,int f1);
  int edgeAddT(int v0,int v1,int f0,int f1);
	
  short findEdgeFace(int e,int f);
  short findEdgeVertex(int e, int v);
  short findFaceEdge(int f, int e);
  short findFaceVertex(int f, int v);
	
  int findVertexNeighbours(int i, vector<int> &V, vector<int> &E, vector<int> &F);
	
  void saveOFF();
  void saveOBJ();
  void saveOFF_test();
  void saveOFF_testT();
	
  int OrienterFaces(int a);
  void PermuteOrient();

  void Flip(int i);
  void FlipR(const vector<int>& hist);
	
  bool ombrelle(const int& i);

  int over(class way A, int v, int v2, int f);
  void addBranch(class way& A, vector<bool>& facette, const vector<bool>& block);
  bool composante_connexe(int e, vector<bool>& block);
  int edgeExist(int v1, int v2);
  int edgeExist(int v1, int v2, vector<int>& result);
  bool edgeExistNoBlock(int v1, int v2, vector<int>& result, const vector<bool>& block);
  bool problem(const int& f0, const int& f1);
  void addBranchStar(way& A, const int& v0, const vector<bool>& block, vector<bool>& facette, const bool& noStar);
  void initWayStar(class way& A, int v1, vector<bool>& facette);
  int BuildEdgeStar(int v0, int v1, vector<bool>& block, const bool& c, const bool& noStar);
  int BuildEdgeStarFromEdgeWithClose(const int& v0, const int& v1, const int& e1, vector<bool>& block, int t);
  int BuildEdgeStarFromEdge(const int& v0, const int& v1, const int& e, vector<bool>& block, const int& tr, bool c);
  int DifferentFaceEdge(int f, int e1, int e2);
  bool CheckOrientation(int f, const int& t);
  int EdgeCorresponding(int f, int center, int u);	
  int troisieme(const int& e1, const int& e2, const vector<bool>& block);
  int CloseTriangle(const int& ed, const int& u, const vector<bool>& block);
  bool buildTriangle(const int& t, vector<bool>& block, vector<int>& FaceOk, vector<int>& vertexR);
  int CloseTriangleWithOrientation(const int& ed, const int& u, const int& t, const vector<bool>& block);
  int init(const int& t, vector<bool>& block, const int& color);
  bool NoConstructible(const int& t, vector<int>& FaceOk, vector<bool>& block);
  bool buildTriangleWithModifyGenus(const int& t, vector<bool>& block, vector<int>& FaceOk, vector<int>& vertexR);
  void buildEdgeWithCondition(bool simply);
  void generateFile();
	
  int simplifysequence(); 
  bool commute(const int& a, const int& b);
  bool simplificationNoCommutativity(const int& a, const int& b);
  void orderFlip(vector<int>& comb, int i);
  int Simplify(int i, int& first);
  void simplifyLastFlip(int& gain, int i);
	
  void SimplifySeq();
  void SimplifySubSeq(vector<vector<int> >& seq, int& pos, const int& start, const int& end, const int& col1, const int& col2);
  void EmptyColumn(vector<vector<int> >& seq, const int& i);
  bool simplificationWithNoCommutativity(vector<vector<int> >& seq, int& pos, const int& col, int& i, int& other);
  bool FillOneColumn(vector<vector<int> >& seq, int& pos, const int& col);
  void MovePos(const vector<vector<int> >& seq, int& pos, const int& souhaite);
	
  int chooseVertex(vector<int>& vertexR, const int& v, const int & e, const vector<bool>& block);
	
  void FindSequence();
  void CodageSimultaneousFlip2();
  void ReadCodageSimultaneousFlip2();
};

#endif // MESH_H
