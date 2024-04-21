//smooth when join neighbours2
//remove inactive when IDLE
//put symmetry here inside elementary/simple deformations
//if boxes - test for region>0

//why modified
#include "Mesh.h"
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <deque>
#include "vtkArithmeticCoderBase.h"
#include "QuasiStaticModel.h"

#include <sys/times.h>
#include <time.h>

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::ostringstream;

Mesh::Mesh()
{
    int count=300000;
    vertex.reserve(count);
    edge.reserve(count);
    face.reserve(count);
}

//________________________________________________________________________________________________________________________________

Mesh::~Mesh()
{
}

//________________________________________________________________________________________________________________________________

int Mesh::vertexAdd(float x,float y,float z)
{
    vertex.push_back(Vertex(x,y,z));
    return vertex.size()-1;
}

//________________________________________________________________________________________________________________________________

int Mesh::vertexAddT(float x,float y,float z)
{
    vertexT.push_back(Vertex(x,y,z));
    return vertexT.size()-1;
}

//________________________________________________________________________________________________________________________________

int Mesh::faceAdd(int v0,int v1,int v2,int e0,int e1,int e2)
{
    face.push_back(Face(v0,v1,v2,e0,e1,e2));
    return face.size()-1;
}

//________________________________________________________________________________________________________________________________


int Mesh::faceAddT(int v0,int v1,int v2,int e0,int e1,int e2)
{
    faceT.push_back(Face(v0,v1,v2,e0,e1,e2));
    return faceT.size()-1;
}

//________________________________________________________________________________________________________________________________

int Mesh::edgeAdd(int v0,int v1,int f0,int f1)
{
    edge.push_back(Edge(v0,v1,f0,f1));
    return edge.size()-1;
}

//________________________________________________________________________________________________________________________________

int Mesh::edgeAddT(int v0,int v1,int f0,int f1)
{
    edgeT.push_back(Edge(v0,v1,f0,f1));
    return edgeT.size()-1;
}

//________________________________________________________________________________________________________________________________

short Mesh::findEdgeFace(int e, int f)
{
    if(edge[e].f[0]==f)return 0;
    if(edge[e].f[1]==f)return 1;
    return -1;
}

//________________________________________________________________________________________________________________________________

short Mesh::findEdgeVertex(int e, int v)
{
    if(edge[e].v[0]==v)return 0;
    if(edge[e].v[1]==v)return 1;
    return -1;
}

//________________________________________________________________________________________________________________________________

short Mesh::findFaceEdge(int f, int e)
{
    if(face[f].e[0]==e)return 0;
    if(face[f].e[1]==e)return 1;
    if(face[f].e[2]==e)return 2;
    return -1;
}

//________________________________________________________________________________________________________________________________

short Mesh::findFaceVertex(int f, int v)
{
    if(face[f].v[0]==v)return 0;
    if(face[f].v[1]==v)return 1;
    if(face[f].v[2]==v)return 2;
    return -1;
}

//________________________________________________________________________________________________________________________________

int Mesh::findVertexNeighbours(int i, vector<int> &V, vector<int> &E, vector<int> &F)
{
    V.clear();
    E.clear();
    F.clear();
    int f=vertex[i].f;
    int f0=f;
    do
    {
        F.push_back(f);
        short k=findFaceVertex(f,i);
        int e=face[f].e[(k+1)%3];
        E.push_back(e);
        int v=face[f].v[(k+2)%3];
        V.push_back(v);
        f=edge[e].f[1-findEdgeFace(e,f)];
    } while(f!=f0);
    return V.size();
}

//________________________________________________________________________________________________________________________________

void Mesh::saveOBJ()
//Generer le maillge en format .obj
{
    ofstream file("Final_Mesh.obj");

    for(unsigned int i=0;i<vertex.size();i++)
        file<<"v "<<vertex[i].V.v[0]<<" "<<vertex[i].V.v[1]<<" "<<vertex[i].V.v[2]<<std::endl;

    for(unsigned int i=0;i<face.size();i++)
        file<<"f "<<face[i].v[0]+1<<" "<<face[i].v[1]+1<<" "<<face[i].v[2]+1<<std::endl;

    file.close();
}

//________________________________________________________________________________________________________________________________

void Mesh::saveOFF()
//Generer le maillge en format .off
{
    ofstream file("Final_Mesh.off");

    file<<"OFF"<<std::endl;
    file<<vertex.size()<<" "<<face.size()<<" 0"<<std::endl;

    for(int i=0;i<vertex.size();i++)
        file<<vertex[i].V.v[0]<<" "<<vertex[i].V.v[1]<<" "<<vertex[i].V.v[2]<<std::endl;

    for(int i=0;i<face.size();i++)
        file<<"3 "<<face[i].v[0]<<" "<<face[i].v[1]<<" "<<face[i].v[2]<<std::endl;

    file.close();
}

//________________________________________________________________________________________________________________________________

//void Mesh::loadOBJ(const string &filename)
////Charger un maillage à partir d'un fichier .obj
//{
//    vertex.clear();
//    edge.clear();
//    face.clear();
//	
//    std::string line;
//    ifstream file(filename.c_str());
//    while(getline(file,line))
//    {
//        char first=toupper(line[0]);
//        char second=toupper(line[1]);
//        if((first=='V')&&(second==' '))
//        {
//            line=line.substr(1,line.length()-1);
//            istringstream content(line);
//            float x,y,z;
//            content>>x;
//            content>>y;
//            content>>z;
//            vertexAdd(x,y,z);
//        }
//        else
//            if(first=='F')
//            {
//				QString s=QString::fromStdString(line);
//				QStringList groups=s.split(" ",QString::SkipEmptyParts);
//				int v[3];
//				for(int i=1;i<groups.size();i++)
//				{
//					QStringList numbers=groups.at(i).split("/",QString::SkipEmptyParts);
//					switch(i)
//					{
//						case 1:v[0]=numbers.at(0).toInt();break;
//						case 2:v[1]=numbers.at(0).toInt();break;
//						default:
//						{
//							v[2]=numbers.at(0).toInt();
//							faceAdd(v[0]-1,v[1]-1,v[2]-1);
//							v[1]=v[2];
//						}
//					}
//				}
//			}
//    }
//    file.close();
//}

//________________________________________________________________________________________________________________________________

//void Mesh::loadOFF(const string &filename)
////Charger un maillage à partir d'un fichier .off
//{
//    vertex.clear();
//    edge.clear();
//    face.clear();
//	
//    std::string line;
//    ifstream file(filename.c_str());
//	
//    getline(file,line);
//    getline(file,line);
//    istringstream content(line);
//    int nv,nf;
//    content>>nv>>nf;
//	
//    vector<vector<int> > V;
//    V.resize(nv);
//    for(int i=0;i<nv;i++)
//    {
//        getline(file,line);
//        istringstream content(line);
//        float x,y,z;
//        content>>x>>y>>z;
//        vertexAdd(x,y,z);
//		//        V.push_back(vector<int>());
//    } 
//	
//    for(int i=0;i<nf;i++)
//    {
//        getline(file,line);
//        istringstream content(line);
//        int n,a,b,c;
//        content>>n>>a>>b>>c;
//        faceAdd(a,b,c,-1,-1,-1);
//        V[a].push_back(i);
//        V[b].push_back(i);
//        V[c].push_back(i);
//    }
//    file.close();
//    center();
//}


//______________________________________________________________________________________________________________________________________

void Mesh::PermuteOrient()
//Permuter l'orientation des faces
{
    int temp;
    for(int a=0; a<face.size(); a++)
    {
        temp=face[a].v[0];
        face[a].v[0]=face[a].v[1];
        face[a].v[1]=temp;

        temp=face[a].e[0];
        face[a].e[0]=face[a].e[1];
        face[a].e[1]=temp;
    }
}

//______________________________________________________________________________________________________________________________________


/*********************/
// Fonction for Flip */
/*********************/

void Mesh::Flip(int i)
//Bascule de l'arête i
{
    int f0=edge[i].f[0];
    int f1=edge[i].f[1];
    int e0[3],v0[3],e1[3],v1[3];

    int k0=findFaceEdge(f0,i);
    for(int j=0;j<3;j++)
    {
        e0[j]=face[f0].e[(k0+j)%3];
        v0[j]=face[f0].v[(k0+j)%3];
    }

    int k1=findFaceEdge(f1,i);
    for(int j=0;j<3;j++)
    {
        e1[j]=face[f1].e[(k1+j)%3];
        v1[j]=face[f1].v[(k1+j)%3];
    }


    face[f0].setE(i,e0[2],e1[1]);
    face[f1].setE(i,e1[2],e0[1]);
    edge[e1[1]].f[findEdgeFace(e1[1],f1)]=f0;
    edge[e0[1]].f[findEdgeFace(e0[1],f0)]=f1;

    face[f0].setV(v0[1],v1[0],v0[0]);
    face[f1].setV(v1[1],v0[0],v1[0]);


    for(int j=0;j<3;j++)vertex[face[f0].v[j]].f=f0;
    vertex[v1[1]].f=f1;

    edge[i].setV(v0[0],v1[0]);
}

//________________________________________________________________________________________________________________________________


void Mesh::FlipR(const vector<int>& hist)
//Revenir au maillage avant bascule à partir du tableau hist
{
    for(unsigned int i=0; i<hist.size(); i=i+22)
    {
        for(unsigned int w=0; w<3; w++)
        {
            face[hist[i]].v[w]=hist[i+7*w+1];
            vertex[hist[i+7*w+1]].f=hist[i+7*w+2];
            face[hist[i]].e[w]=hist[i+7*w+3];
            for(unsigned int d=0; d<2; d++)
            {
                edge[hist[i+7*w+3]].v[d]=hist[i+7*w+4+d*2];
                edge[hist[i+7*w+3]].f[d]=hist[i+7*w+5+d*2];
            }
        }
    }
}

//___________________________________________________________________________________________________________________________________________


bool Mesh::ombrelle(const int& i)
{
    int som=0;
    for(int e=0; e<edge.size(); e++)
    {
        if(edge[e].v[0]==i || edge[e].v[1]==i)
            som++;
    }
    int f=vertex[i].f;
    int f0=f;
    int count=0;
    do
    {
        short k=findFaceVertex(f,i);
        int e=face[f].e[(k+1)%3];
        if(e==-1)
        {
            return true;
        }
        f=edge[e].f[1-findEdgeFace(e,f)];
        count++;
    } while(f!=f0 && count<edge.size());
    if(count==edge.size())
    {
        return true;
    }
    if(count!=som)
    {
        return false;
    }
    return true;
}

/****************************/
// Find Edge Flips Sequence //
/****************************/

//__________________________________________________________________________________________________________________________________


int Mesh::over(class way A, int v, int v2, int f)
//Fonction pour trouver un chemin de faces entre les sommets v et v2
{
	//std::cout << " Over : " << "(" << v << " , " << v2 << ") - A.size()=" << A.tree.size() << std::endl;
    if(f<A.tree.size())
    {
        for(unsigned int i=f; i<A.tree.size(); i=i+5)
        {
            if(A.tree[i+2]==v)
            {
                if(i>0)
                {
				//	std::cout << " Over : Trouve " << "(" << v << " , " << v2 << ") - A.size()=" << A.tree.size() << std::endl;
                    return i;
                }
            }
        }
    }
    return -1;
}

//__________________________________________________________________________________________________________________________________

void Mesh::addBranch(class way& A, vector<bool>& facette, const vector<bool>& block)
//Fonction pour trouver un chemin de faces
{
    unsigned int minf=A.tree.size();
    for(unsigned int i=A.min; i<minf; i=i+5)
    {
        int m=A.tree[i+3];//j est le dernier triangle traité

        int temp1=findFaceVertex(m,A.tree[i+1]);
        if(temp1>-1)
        {
            temp1=face[m].e[temp1];//arete
            if(temp1>-1 && block[temp1])
            {
                int temp2=findEdgeFace(temp1, m);
                if(temp2>-1)
                {
                    int temp3=edge[temp1].f[(temp2+1)%2];// facette en question
                    if(temp3>-1 && !facette[temp3])
                    {
                        A.tree.push_back(A.tree[i]);
                        A.tree.push_back(A.tree[i+2]);
                            A.tree.push_back(face[temp3].v[findFaceEdge(temp3,temp1)]);
                        A.tree.push_back(temp3);
                        A.tree.push_back(i);
                        facette[temp3]=true;
                    }
                }
            }
        }

        temp1=findFaceVertex(m,A.tree[i]);
        if(temp1>-1)
        {
            temp1=face[m].e[temp1];//arete
            if(temp1>-1 && block[temp1])
            {
                int temp2=findEdgeFace(temp1, m);
                if(temp2>-1)
                {
                    int temp3=edge[temp1].f[(temp2+1)%2];// facette en question
                    if(temp3>-1 && !facette[temp3])
                    {
                        A.tree.push_back(A.tree[i+1]);
                        A.tree.push_back(A.tree[i+2]);
                        A.tree.push_back(face[temp3].v[findFaceEdge(temp3,temp1)]);
                        A.tree.push_back(temp3);
                        A.tree.push_back(i);
                        facette[temp3]=true;
                    }
                }
            }
        }
    }
    A.min=minf;
}

//________________________________________________________________________________________________________________________________


bool present(class way A, int f)
//Fonction pour trouver un chemin de faces
{
    for(unsigned int i=A.min; i<A.tree.size(); i=i+5)
        if(A.tree[i+3]==f)
            return false;
    return true;
}

//__________________________________________________________________________________________________________________________________


bool Mesh::composante_connexe(int e, vector<bool>& block)
{
    vector<bool> facette;
    facette.resize(face.size(),false);
    vector<bool> facette2;
    facette2.resize(face.size(),false);

    int f=edge[e].f[1];//triangle voisin à e numéro 1
    int f0=edge[e].f[0];//triangle voisin à e numéro 0
	
    class way A;
    int temp=findFaceEdge(f, e);
    A.tree.push_back(face[f].v[(temp+1)%3]);
    A.tree.push_back(face[f].v[(temp+2)%3]);
    A.tree.push_back(face[f].v[temp]);
    A.tree.push_back(f);
    A.tree.push_back(-1);
    A.min=0;

    facette[f]=true;

    if(block[face[f].e[(temp+1)%3]] || block[face[f].e[(temp+2)%3]])
    {
        bool b=true;
        while(A.tree.size()>A.min && b)
        {
            addBranch(A, facette, block);
            b=present(A, f0);

        }
        if(!b)
        {
            return true;
        }
    }
    return false;
}

//________________________________________________________________________________________________________________________________


int Mesh::edgeExist(int v1, int v2)
//Retourne l'arête reliant v1 et v2 si elle existe. Sinon retroune -1
{
    unsigned int count=0;
    int f=vertex[v1].f;
    int f0=f;
    do
    {
        short k=findFaceVertex(f,v1);
        int e=face[f].e[(k+1)%3];
        int v=face[f].v[(k+2)%3];
        if(v==v2)
            return e;
        f=edge[e].f[1-findEdgeFace(e,f)];
        count++;
    } while(f!=f0 && count<edge.size());
    if(count==edge.size())
    {
        for(unsigned int i=0; i<edge.size(); i++)
        {
            if((edge[i].v[0]==v1 && edge[i].v[1]==v2) || (edge[i].v[0]==v2 && edge[i].v[1]==v1))
            {
                return i;
            }
        }
    }
    return -1;
}

//________________________________________________________________________________________________________________________________

int Mesh::edgeExist(int v1, int v2, vector<int>& result)
//Enregistre les indices des arêtes reliant v1 et v2 dans le tableau result si elles existent. Sinon retroune -1
{
    result.clear();
    unsigned int count=0;
    int f=vertex[v1].f;
    int f0=f;
    do
    {
        short k=findFaceVertex(f,v1);
        int e=face[f].e[(k+1)%3];
        int v=face[f].v[(k+2)%3];
        if(v==v2)
            result.push_back(e);
        f=edge[e].f[1-findEdgeFace(e,f)];
        count++;
    } while(f!=f0 && count<edge.size());

    if(count==edge.size())
    {
        result.clear();
        for(unsigned int i=0; i<edge.size(); i++)
        {
            if((edge[i].v[0]==v1 && edge[i].v[1]==v2) || (edge[i].v[0]==v2 && edge[i].v[1]==v1))
            {
                result.push_back(i);
            }
        }
    }
    return result.size();
}

//________________________________________________________________________________________________________________________________


bool Mesh::edgeExistNoBlock(int v1, int v2, vector<int>& result, const vector<bool>& block)
//Enregistre les indices des arêtes non bloquées reliant v1 et v2 dans le tableau result si elles existent. Sinon retroune -1
{
    result.clear();
    unsigned int count=0;
    int f=vertex[v1].f;
    int f0=f;
    do
    {
        short k=findFaceVertex(f,v1);
        int e=face[f].e[(k+1)%3];
        int v=face[f].v[(k+2)%3];
        if(v==v2)
        {
            if(block[e])
				result.push_back(e);
        }
        f=edge[e].f[1-findEdgeFace(e,f)];
        count++;
    } while(f!=f0 && count<edge.size());
    if(count==edge.size())
    {
        result.clear();
        for(unsigned int i=0; i<edge.size(); i++)
        {
            if((edge[i].v[0]==v1 && edge[i].v[1]==v2) || (edge[i].v[0]==v2 && edge[i].v[1]==v1) && block[i])
                result.push_back(i);
        }

        count=0;
        f=vertex[v1].f;
        f0=f;
        do
        {
            short k=findFaceVertex(f,v1);
            int e=face[f].e[(k+1)%3];
            f=edge[e].f[1-findEdgeFace(e,f)];
            count++;
        } while(f!=f0 && count<20);

    }
    return !result.empty();
}

//______________________________________________________________________________________________________________________________


bool Mesh::problem(const int& f0, const int& f1)
//Arete que l'on interdit de basculer
{
    if((face[f0].v[0]==face[f1].v[0] && face[f0].v[1]==face[f1].v[1] && face[f0].v[2]==face[f1].v[2]) || (face[f0].v[0]==face[f1].v[0] && face[f0].v[1]==face[f1].v[2] && face[f0].v[2]==face[f1].v[1]) || (face[f0].v[0]==face[f1].v[1] && face[f0].v[1]==face[f1].v[2] && face[f0].v[2]==face[f1].v[0]) || (face[f0].v[0]==face[f1].v[1] && face[f0].v[1]==face[f1].v[0] && face[f0].v[2]==face[f1].v[2]) || (face[f0].v[0]==face[f1].v[2] && face[f0].v[1]==face[f1].v[0] && face[f0].v[2]==face[f1].v[1]) || (face[f0].v[0]==face[f1].v[2] && face[f0].v[1]==face[f1].v[1] && face[f0].v[2]==face[f1].v[0]))
        return true;
    return false; 
}


//________________________________________________________________________________________________________________________________

void Mesh::addBranchStar(way& A, const int& v0, const vector<bool>& block, vector<bool>& facette, const bool& noStar)
//Fonction pour trouver un chemin de faces
{
	//std::cout << " Add Branch star (" << v0 << ") A.size()=" << A.tree.size() << std::endl;
    if(!noStar)
    {
        unsigned int minf=A.tree.size();
        for(unsigned int i=A.min; i<minf; i=i+5)
        {
            if(A.tree[i+2]!=v0)
            {
                int m=A.tree[i+3];//j est le dernier triangle
                int temp1=findFaceVertex(m,A.tree[i+1]);
                if(temp1>-1)
                {
                    temp1=face[m].e[temp1];//arete
                    if(temp1>-1 && block[temp1])
                    {
                        int temp2=findEdgeFace(temp1, m);
                        if(temp2>-1)
                        {
                            int temp3=edge[temp1].f[(temp2+1)%2];// facette en question
                            if(temp3>-1 && facette[temp3])
                            {
                                A.tree.push_back(A.tree[i]);
                                A.tree.push_back(A.tree[i+2]);
                                A.tree.push_back(face[temp3].v[findFaceEdge(temp3,temp1)]);
                                A.tree.push_back(temp3);
                                A.tree.push_back(i);
                                facette[temp3]=false;
                            }
                        }
                    }
                }

                temp1=findFaceVertex(m,A.tree[i]);
                if(temp1>-1)
                {

                    temp1=face[m].e[temp1];//arete
                    if(temp1>-1 && block[temp1])
                    {
                        int temp2=findEdgeFace(temp1, m);
                        if(temp2>-1)
                        {
                            int temp3=edge[temp1].f[(temp2+1)%2];// facette en question
                            if(temp3>-1 && facette[temp3])
                            {
                                A.tree.push_back(A.tree[i+1]);
                                A.tree.push_back(A.tree[i+2]);
                                A.tree.push_back(face[temp3].v[findFaceEdge(temp3,temp1)]);
                                A.tree.push_back(temp3);
                                A.tree.push_back(i);
                                facette[temp3]=false;
                            }
                        }
                    }
                }
            }
        }
        A.min=minf;
    }
    else
    {
        unsigned int minf=A.tree.size();
        for(unsigned int i=A.min; i<minf; i=i+5)
        {
            if(A.tree[i+2]!=v0)
            {
                int m=A.tree[i+3];//j est le dernier triangle
                int temp1=findFaceVertex(m,A.tree[i+1]);
                if(temp1>-1)
                {
                    temp1=face[m].e[temp1];//arete
                    if(temp1>-1 && block[temp1])
                    {
                        int temp2=findEdgeFace(temp1, m);
                        if(temp2>-1)
                        {
                            int temp3=edge[temp1].f[(temp2+1)%2];// facette en question
                            bool b=true; //FLAG
                            int r=i;
                            while(r>-1)
                            {
                                if(A.tree[r+3]==temp3)
                                {
                                    b=false;
                                }
                                r=A.tree[r+4];
                            }
                            if(temp3>-1 && b && facette[temp3])
                            {
                                A.tree.push_back(A.tree[i]);
                                A.tree.push_back(A.tree[i+2]);
                                A.tree.push_back(face[temp3].v[findFaceEdge(temp3,temp1)]);
                                A.tree.push_back(temp3);
                                A.tree.push_back(i);
                                facette[temp3]=false;
                            }
                        }
                    }
                }

                temp1=findFaceVertex(m,A.tree[i]);
                if(temp1>-1)
                {
                    temp1=face[m].e[temp1];//arete
                    if(temp1>-1 && block[temp1])
                    {
                        int temp2=findEdgeFace(temp1, m);
                        if(temp2>-1)
                        {
                            int temp3=edge[temp1].f[(temp2+1)%2];// facette en question
                            bool b=true; //FLAG
                            int r=i;
                            while(r>-1)
                            {
                                if(A.tree[r+3]==temp3)
                                {
                                    b=false;
                                }
                                r=A.tree[r+4];
                            }
                            if(temp3>-1 && b && facette[temp3])
                            {
                                A.tree.push_back(A.tree[i+1]);
                                A.tree.push_back(A.tree[i+2]);
                                A.tree.push_back(face[temp3].v[findFaceEdge(temp3,temp1)]);
                                A.tree.push_back(temp3);
                                A.tree.push_back(i);
                                facette[temp3]=false;
                            }
                        }
                    }
                }
            }
        }
        A.min=minf;
    }
}

//________________________________________________________________________________________________________________________________

void Mesh::initWayStar(class way& A, int v1, vector<bool>& facette)
//Fonction pour trouver un chemin de faces
{
    unsigned int count=0;
    int f=vertex[v1].f;
    int f0=f;
    do
    {
        short temp=findFaceVertex(f,v1);

        A.tree.push_back(face[f].v[(temp+1)%3]);
        A.tree.push_back(face[f].v[temp]);
        A.tree.push_back(face[f].v[(temp+2)%3]);
        A.tree.push_back(f);
        facette[f]=false;
        A.tree.push_back(-1);

        int e=face[f].e[(temp+1)%3];
        f=edge[e].f[1-findEdgeFace(e,f)];

        count++;
    } while(f!=f0 && count<edge.size());
    A.min=0;
}

//________________________________________________________________________________________________________________________________


int Mesh::BuildEdgeStar(int v0, int v1, vector<bool>& block, const bool& c, const bool& noStar)
//Relier deux sommets v0 et v1
{
    vector<int> seq;
    bool check;
    class way A;
    vector<bool> facette;
    facette.resize(edge.size(), true);
    A.tree.reserve(face.size());
    int temp=findFaceVertex(vertex[v0].f, v0);
    int e =face[vertex[v0].f].e[(temp+1)%3];
    A.min=0;
    int f0;
    if(block[face[edge[e].f[0]].e[0]] || block[face[edge[e].f[0]].e[1]] || block[face[edge[e].f[0]].e[2]])
    {
        f0=edge[e].f[0];
    }
    else
    {
        f0=edge[e].f[1];
    }
    int tem;
    if(face[f0].v[0]!=edge[e].v[0] && face[f0].v[0]!=edge[e].v[1])
        tem=0;
    else
    {
        if(face[f0].v[1]!=edge[e].v[0] && face[f0].v[1]!=edge[e].v[1])
            tem=1;
        else
            tem=2;
    }
    A.tree.push_back(face[f0].v[(tem+1)%3]);
    A.tree.push_back(face[f0].v[(tem+2)%3]);
    A.tree.push_back(face[f0].v[tem]);
    A.tree.push_back(f0);
    facette[f0]=false;
    A.tree.push_back(-1);

    int f=-1;
    if(f>-1)
    {
        f=over(A, v1, v0, f+5);
    }
	//std::cout << "ici3" << std::endl;
    while(f<0 && A.tree.size()>A.min)
    {
        addBranchStar(A, v1, block, facette, 0);
        f=over(A, v1, v0, A.min);
    }
    if(f>-1)
    {
        vector<int> afaire;
        int ind=f;
        int e;
        vector<int> hist;
        while(A.tree[ind]!=v0 && A.tree[ind+1]!=v0 && A.tree[ind+2]!=v0 && ind>-1)
        {
            e=face[A.tree[ind+3]].e[findFaceVertex(A.tree[ind+3], A.tree[ind+2])];
            afaire.push_back(e);
            for(unsigned int o=0; o<2; o++)
            {
                hist.push_back(edge[e].f[o]);
                for(unsigned int u=0; u<3; u++)
                {
                    int fre=face[edge[e].f[o]].v[u];
                    hist.push_back(fre);
                    hist.push_back(vertex[fre].f);
                    fre=face[edge[e].f[o]].e[u];
                    hist.push_back(fre);
                    for(unsigned int y=0; y<2; y++)
                    {
                        hist.push_back(edge[fre].v[y]);
                        hist.push_back(edge[fre].f[y]);
                    }
                }
            }
            ind=A.tree[ind+4];
        }
        while(!afaire.empty() && (edge[afaire[afaire.size()-1]].v[0]==v0 || edge[afaire[afaire.size()-1]].v[1]==v0))
        {
            afaire.pop_back();
        }
        bool pr=false;
        for(int x=0; x<afaire.size(); x++)
        {
            if(problem(edge[afaire[x]].f[0], edge[afaire[x]].f[1]))
                pr=true;
            Flip(afaire[x]);
        }

        if(!afaire.empty() && ((edge[afaire[afaire.size()-1]].v[0]==v0 && edge[afaire[afaire.size()-1]].v[1]==v1) || (edge[afaire[afaire.size()-1]].v[1]==v0 && edge[afaire[afaire.size()-1]].v[0]==v1)))
        {
            for(unsigned int x=0; x<afaire.size(); x++)
                comb.push_back(afaire[x]);
            return afaire[afaire.size()-1];
        }
        else
        {
            FlipR(hist);
            hist.clear();
        }
    }
    return -1;
}

//________________________________________________________________________________________________________________________________

int Mesh::BuildEdgeStarFromEdge(const int& v0, const int& v1, const int& e, vector<bool>& block, const int& tr, bool c)
//Relier un sommet v0 à une extrémité de l'arête e
{
	//std::cout << "BuildEdgeStarFromEdge - (" << v0 << " , " << v1 << ") from e=" << e << std::endl;
    vector<int> seq;
    int max=33554431;
    class way A;
    A.tree.reserve(face.size());
    vector<bool> facette;
    facette.resize(face.size(), true);
    A.min=0;
    int f0;
    if(block[face[edge[e].f[0]].e[0]] || block[face[edge[e].f[0]].e[1]] || block[face[edge[e].f[0]].e[2]])
    {
        f0=edge[e].f[0];
    }
    else
    {
        f0=edge[e].f[1];
    }

    int tem;
    if(face[f0].v[0]!=edge[e].v[0] && face[f0].v[0]!=edge[e].v[1])
        tem=0;
    else
    {
        if(face[f0].v[1]!=edge[e].v[0] && face[f0].v[1]!=edge[e].v[1])
            tem=1;
        else
            tem=2;
    }

    A.tree.push_back(face[f0].v[(tem+1)%3]);
    A.tree.push_back(face[f0].v[(tem+2)%3]);
    A.tree.push_back(face[f0].v[tem]);
    A.tree.push_back(f0);
    facette[f0]=false;
    A.tree.push_back(-1);

    int f=-2;
    while(f!=-1)
    {
        if(f>-1)
        {
            f=over(A, v1, v0, f+5);
        }
        while(f<0 && A.tree.size()>A.min && (A.tree.size()<max))
        {
            addBranchStar(A, v1, block, facette, 1);
            f=over(A, v1, v0, A.min);
        }
        if(f>-1)
        {
            int ind=f;
            int e;
            vector<int> hist;
            while(A.tree[ind+4]>-1)
            {
                e=face[A.tree[ind+3]].e[findFaceVertex(A.tree[ind+3], A.tree[ind+2])];
                seq.push_back(e);
                for(unsigned int o=0; o<2; o++)
                {
                    hist.push_back(edge[e].f[o]);
                    for(unsigned int u=0; u<3; u++)
                    {
                        int fre=face[edge[e].f[o]].v[u];
                        hist.push_back(fre);
                        hist.push_back(vertex[fre].f);
                        fre=face[edge[e].f[o]].e[u];
                        hist.push_back(fre);
                        for(unsigned int y=0; y<2; y++)
                        {
                            hist.push_back(edge[fre].v[y]);
                            hist.push_back(edge[fre].f[y]);
                        }
                    }
                }
                ind=A.tree[ind+4];
            }
            while(!seq.empty() && (edge[seq[seq.size()-1]].v[0]==v0 || edge[seq[seq.size()-1]].v[1]==v0) /*(face[edge[seq[seq.size()-1]].f[0]].v[0]==v0  || face[edge[seq[seq.size()-1]].f[0]].v[1]==v0 || face[edge[seq[seq.size()-1]].f[0]].v[2]==v0) && (face[edge[seq[seq.size()-1]].f[1]].v[0]==v0  || face[edge[seq[seq.size()-1]].f[1]].v[1]==v0 || face[edge[seq[seq.size()-1]].f[1]].v[2]==v0)*/)
                seq.pop_back();

            bool pr=false;
            for(int x=0; x<seq.size(); x++)
            {
                if(problem(edge[seq[x]].f[0], edge[seq[x]].f[1]))
                    pr=true;
                Flip(seq[x]);
            }

            if(!seq.empty() && ((edge[seq[seq.size()-1]].v[0]==v0 && edge[seq[seq.size()-1]].v[1]==v1) || (edge[seq[seq.size()-1]].v[1]==v0 && edge[seq[seq.size()-1]].v[0]==v1)) && composante_connexe(seq[seq.size()-1], block))
            {
                for(unsigned int x=0; x<seq.size(); x++)
                    comb.push_back(seq[x]);
                return seq[seq.size()-1];
            }
            FlipR(hist);
            hist.clear();
            seq.clear();
            if(c)
                f=-1;
        }
    }
    if(c)
    {
        c=false;
        return BuildEdgeStarFromEdge(v0, v1, e, block, tr, c);
    }
	//std::cout << "BuildEdgeStarFromEdge - echec" << std::endl;
    return -1;
}

//________________________________________________________________________________________________________________________________

int Mesh::BuildEdgeStarFromEdgeWithClose(const int& v0, const int& v1, const int& e1, vector<bool>& block, int t)
//Construction d'une face vers un sommet incident à aucune arête bloquée
{
	//std::cout << " [BuildEdgeStarFromEdgeWithClose]" << std::endl;
    int max=500005;
    class way A;
    A.tree.reserve(face.size());
    vector<bool> facette;
    facette.resize(face.size(), true);

    A.min=0;
    int f0;
    if(block[face[edge[e1].f[0]].e[0]] || block[face[edge[e1].f[0]].e[1]] || block[face[edge[e1].f[0]].e[2]])
    {
        f0=edge[e1].f[0];
    }
    else
    {
        f0=edge[e1].f[1];
    }

    int tem;
    if(face[f0].v[0]!=edge[e1].v[0] && face[f0].v[0]!=edge[e1].v[1])
        tem=0;
    else
    {
        if(face[f0].v[1]!=edge[e1].v[0] && face[f0].v[1]!=edge[e1].v[1])
            tem=1;
        else
            tem=2;
    }
    A.tree.push_back(face[f0].v[(tem+1)%3]);
    A.tree.push_back(face[f0].v[(tem+2)%3]);
    A.tree.push_back(face[f0].v[tem]);
    A.tree.push_back(f0);
    facette[f0]=false;
    A.tree.push_back(-1);

//std::cout << " 1" << std::endl;
    int f=-1;
    while(f<0 && A.tree.size()>A.min && A.tree.size()<max)
    {
        addBranchStar(A, v1, block, facette, 0);
        f=over(A, v1, v0, A.min);
		//std::cout << " x";
    }
    if(f>-1)
    {
		//std::cout << " chamin trouve avec A.tree.size()=" << A.tree.size() << std::endl;
        vector<int> afaire;
        int ind=f;
        int e;
        vector<int> hist;
        while(A.tree[ind+4]>-1)
        {
            e=face[A.tree[ind+3]].e[findFaceVertex(A.tree[ind+3], A.tree[ind+2])];
            afaire.push_back(e);
            for(unsigned int o=0; o<2; o++)
            {
                hist.push_back(edge[e].f[o]);
                for(unsigned int u=0; u<3; u++)
                {
                    int fre=face[edge[e].f[o]].v[u];
                    hist.push_back(fre);
                    hist.push_back(vertex[fre].f);
                    fre=face[edge[e].f[o]].e[u];
                    hist.push_back(fre);
                    for(unsigned int y=0; y<2; y++)
                    {
                        hist.push_back(edge[fre].v[y]);
                        hist.push_back(edge[fre].f[y]);
                    }
                }
            }
            ind=A.tree[ind+4];
        }
		//std::cout << " f=" << f << std::endl;
        while(!afaire.empty() && (edge[afaire[afaire.size()-1]].v[0]==v0 || edge[afaire[afaire.size()-1]].v[1]==v0))
        {
            afaire.pop_back();
        }

        for(int x=0; x<afaire.size(); x++)
        {
            Flip(afaire[x]);
        }
//std::cout << " 2" << std::endl;
        if(!afaire.empty() && ((edge[afaire[afaire.size()-1]].v[0]==v0 && edge[afaire[afaire.size()-1]].v[1]==v1) || (edge[afaire[afaire.size()-1]].v[1]==v0 && edge[afaire[afaire.size()-1]].v[0]==v1)))
        {
            block[afaire[afaire.size()-1]]=false;
            for(unsigned int x=0; x<afaire.size(); x++)
                comb.push_back(afaire[x]);
//std::cout << " AVant Close" << std::endl;
            int r=CloseTriangle(afaire[afaire.size()-1], e1, block);
            if(r>-1)
            {
                block[r]=false;
                if((face[edge[r].f[0]].e[0]==e1 || face[edge[r].f[0]].e[1]==e1 || face[edge[r].f[0]].e[2]==e1) &&
                   (face[edge[r].f[0]].e[0]==afaire[afaire.size()-1] || face[edge[r].f[0]].e[1]==afaire[afaire.size()-1] || face[edge[r].f[0]].e[2]==afaire[afaire.size()-1]))
                        return edge[r].f[0];
                else
                {
                    if((face[edge[r].f[1]].e[0]==e1 || face[edge[r].f[1]].e[1]==e1 || face[edge[r].f[1]].e[2]==e1) &&
                       (face[edge[r].f[1]].e[0]==afaire[afaire.size()-1] || face[edge[r].f[1]].e[1]==afaire[afaire.size()-1] || face[edge[r].f[1]].e[2]==afaire[afaire.size()-1]))
                        return edge[r].f[1];
                    else
                    {
						//std::cout << " echec" << std::endl;
                        return -1;
                    }
                }
            }
            for(unsigned int x=0; x<afaire.size(); x++)
                comb.pop_back();
            block[afaire[afaire.size()-1]]=true;
        }
        FlipR(hist);
        hist.clear();
    }
    else
    {
        for(int y=0; y<A.tree.size(); y=y+5)
        {
            for(int v=0; v<3; v++)
            {
                if(block[face[A.tree[y+3]].e[v]])
                {
                    int ert=face[A.tree[y+3]].e[v];
                    if(edge[ert].f[0]==A.tree[y+3])
                    {
                        bool bn=true;
                        for(int g=0; g<A.tree.size(); g=g+5)
                        {
                            if(edge[ert].f[1]==A.tree[g+3])
                                bn=false;
                        }
                    }
                    else
                    {
                        if(edge[ert].f[1]==A.tree[y+3])
                        {
                            bool bn=true;
                            for(int g=0; g<A.tree.size(); g=g+5)
                            {
                                if(edge[ert].f[0]==A.tree[g+3])
                                    bn=false;
                            }
                        }
                    }
                }
            }
        }
    }
    return -1;
}

//_________________________________________________________________________________________________________________________________


int Mesh::DifferentFaceEdge(int f, int e1, int e2)
{
    if(face[f].e[0]!=e1 && face[f].e[0]!=e2)
        return face[f].e[0];
    else
    {
        if(face[f].e[1]!=e1 && face[f].e[1]!=e2)
            return face[f].e[1];
        else
        {
            if(face[f].e[2]!=e1 && face[f].e[2]!=e2)
                return face[f].e[2];
            else
            {
                return face[f].e[2];
            }
        }
    }
}

//_________________________________________________________________________________________________________________________________

bool Mesh::CheckOrientation(int f, const int& t)
{
    int i;
    if(face[f].v[0]==vertexR[faceT[t].v[0]])
        i=0;
    else
    {
        if(face[f].v[0]==vertexR[faceT[t].v[1]])
            i=1;
        else
        {
            if(face[f].v[0]==vertexR[faceT[t].v[2]])
                i=2;
            else
                return false;
        }
    }
    if(face[f].v[1]==vertexR[faceT[t].v[(i+1)%3]] && face[f].v[2]==vertexR[faceT[t].v[(i+2)%3]])
    {
        return true;
    }
    return false;
}

//_________________________________________________________________________________________________________________________________

int Mesh::EdgeCorresponding(int f, int center, int u)
{
    if((edge[face[f].e[0]].v[0]==center || edge[face[f].e[0]].v[1]==center) && face[f].e[0]!=u)
        return face[f].e[0];
    else
    {
        if((edge[face[f].e[1]].v[0]==center || edge[face[f].e[1]].v[1]==center) && face[f].e[1]!=u)
            return face[f].e[1];
        else
        {
            if((edge[face[f].e[2]].v[0]==center || edge[face[f].e[2]].v[1]==center) && face[f].e[2]!=u)
                return face[f].e[2];
            else
                return face[f].e[2];
        }
    }
}

//_________________________________________________________________________________________________________________________________

int Mesh::troisieme(const int& e1, const int& e2, const vector<bool>& block)
{
    int f;
    if(edge[e1].f[0]==edge[e2].f[0] || edge[e1].f[0]==edge[e2].f[1])
        f=edge[e1].f[0];
    else
    {
        if(edge[e1].f[1]==edge[e2].f[0] || edge[e1].f[1]==edge[e2].f[1])
            f=edge[e1].f[1];
        else
        {
            return -1;
        }
    }
    if((face[f].e[0]==e1 && face[f].e[1]==e2) || (face[f].e[0]==e2 && face[f].e[1]==e1))
    {
        if(block[face[f].e[2]])
            return face[f].e[2];
    }
    else
    {
        if((face[f].e[0]==e1 && face[f].e[2]==e2) || (face[f].e[0]==e2 && face[f].e[2]==e1))
        {
            if(block[face[f].e[1]])
                return face[f].e[1];
        }
        else
        {
            if((face[f].e[2]==e1 && face[f].e[1]==e2) || (face[f].e[2]==e2 && face[f].e[1]==e1))
            {
                if(block[face[f].e[0]])
                    return face[f].e[0];
            }
        }
    }
    return -1;
}

//_________________________________________________________________________________________________________________________________

int Mesh::CloseTriangle(const int& ed, const int& u, const vector<bool>& block)
//Construction d'une face à partir de deux arêtes bloquées
{
    for(int i=0; i<2; i++)
        if(findFaceEdge(edge[ed].f[i], u)>-1)
            for(int y=0; y<3; y++)
                if(face[edge[ed].f[i]].e[y]!=ed && face[edge[ed].f[i]].e[y]!=u)
                    if(block[face[edge[ed].f[i]].e[y]])
                        return face[edge[ed].f[i]].e[y];
	
    vector<int> hist;
    vector<int> afaire;
    vector<int> mem;
    int center;
    bool check;
    int dfg=0;
    bool bbb=true;
    int tempu, temped;
    if(edge[u].v[0]==edge[ed].v[0])
    {
        tempu=0;
        temped=0;
        center=edge[u].v[0];
    }
    else
    {
        if(edge[u].v[0]==edge[ed].v[1])
        {
            tempu=0;
            temped=1;
            center=edge[u].v[0];
        }
        else
        {
            if(edge[u].v[1]==edge[ed].v[0])
            {
                tempu=1;
                temped=0;
                center=edge[u].v[1];
            }
            else
            {
                if(edge[u].v[1]==edge[ed].v[1])
                {
                    tempu=1;
                    temped=1;
                    center=edge[u].v[1];
                }
            }
        }
    }
    int f0=edge[u].f[0];
    bool continu=true;
    int temp=u;
    mem.push_back(u);
    int fg=f0;
    while(continu && temp!=ed)
    {
        temp=EdgeCorresponding(fg, center, temp);
		mem.push_back(temp);
        if(temp!=ed)
        {
            if(!block[temp])
            {
                continu=false;
            }
            else
            {
                afaire.push_back(temp);
                if(edge[temp].f[0]==fg)
                {
                    fg=edge[temp].f[1];
                }
                else
                {
                    if(edge[temp].f[1]==fg)
                    {
                        fg=edge[temp].f[0];
                    }
                    else
                    {
                        fg=edge[temp].f[0];
                    }
                }
            }
        }
    }
    if(continu)
    {
        dfg++;
        for(unsigned int i=0; i< afaire.size(); i++)
        {
            for(unsigned int o=0; o<2; o++)
            {
                hist.push_back(edge[afaire[i]].f[o]);
                for(unsigned int ui=0; ui<3; ui++)
                {
                    int fre=face[edge[afaire[i]].f[o]].v[ui];
                    hist.push_back(fre);
                    hist.push_back(vertex[fre].f);
                    fre=face[edge[afaire[i]].f[o]].e[ui];
                    hist.push_back(fre);
                    for(unsigned int y=0; y<2; y++)
                    {
                        hist.push_back(edge[fre].v[y]);
                        hist.push_back(edge[fre].f[y]);
                    }
                }
            }
        }
		
		
        int x=0;
        check=true;
        bool ghj=false;
        vector<int> seq;
        while(x<afaire.size())
        {
            while(x<afaire.size() && !problem(edge[afaire[x]].f[0], edge[afaire[x]].f[1]))
            {
                Flip(afaire[x]);
                seq.push_back(afaire[x]);
                x++;
            }
			
            if(x<afaire.size())
            {
                ghj=true;
                int y=x;
                while(y<afaire.size() && problem(edge[afaire[y]].f[0], edge[afaire[y]].f[1]))
                    y++;
                if(y==afaire.size())
                {
                    int y=x+1;
                    bool cont=true;
                    int tempc;
                    while(y<=afaire.size() && cont)
                    {
                        if(y<afaire.size())
                            tempc=troisieme(afaire[y-1], afaire[y], block);
                        else
                            tempc=troisieme(afaire[y-1], ed, block);
                        if(tempc!=-1)
                        {
                            Flip(tempc);
                            seq.push_back(tempc);
                            for(int h=y-1; h>x-1; h--)
                            {
                                if(!problem(edge[afaire[h]].f[0], edge[afaire[h]].f[1]))
                                {
                                    Flip(afaire[h]);
                                    seq.push_back(afaire[h]);
                                }
                            }
                            cont=false;
                        }
                        else
                        {
                            y++;
                        }
                    }
					
                    if(y>afaire.size())
                    {
                        check=false;
                    }
                }
                else
                {
                    for(int h=y; h>x-1; h--)
                    {
                        Flip(afaire[h]);
                        seq.push_back(afaire[h]);
                    }
                }
				
                afaire.clear();
                f0=edge[u].f[0];
                continu=true;
                temp=u;
                fg=f0;
                x=0;
                while(continu && temp!=ed)
                {
                    temp=EdgeCorresponding(fg, center, temp);
                    if(temp!=ed)
                    {
                        if(!block[temp])
                        {
                            continu=false;
                        }
                        else
                        {
                            afaire.push_back(temp);
                            if(edge[temp].f[0]==fg)
                            {
                                fg=edge[temp].f[1];
                            }
                            else
                            {
                                if(edge[temp].f[1]==fg)
                                    fg=edge[temp].f[0];
                                else
                                    fg=edge[temp].f[0];
                            }
                        }
                    }
                }
            }
        }
		
        int cas;
        if(seq.empty())
        {
            cas=DifferentFaceEdge(f0, u, ed);
        }
        else
            cas=seq[seq.size()-1];
        if(check && block[cas] && ((edge[cas].v[0]==edge[u].v[(tempu+1)%2] && edge[cas].v[1]==edge[ed].v[(temped+1)%2]) || (edge[cas].v[1]==edge[u].v[(tempu+1)%2] && edge[cas].v[0]==edge[ed].v[(temped+1)%2])))
        {
            for(unsigned int o=0; o<seq.size(); o++)
            {
                comb.push_back(seq[o]);
            }
            return cas;
        }
        FlipR(hist);
        bbb=false;
    }
	
	mem.clear();
    afaire.clear();
    hist.clear();
    f0=edge[u].f[1];
    continu=true;
    temp=u;
    fg=f0;
    mem.push_back(temp);
    while(continu && temp!=ed)
    {
        temp=EdgeCorresponding(fg, center, temp);
        mem.push_back(temp);
		
        if(temp!=ed)
        {
            if(!block[temp])
            {
                continu=false;
            }
            else
            {
                afaire.push_back(temp);
                if(edge[temp].f[0]==fg)
                {
                    fg=edge[temp].f[1];
                }
                else
                {
                    if(edge[temp].f[1]==fg)
                        fg=edge[temp].f[0];
                    else
                    {
                        fg=edge[temp].f[0];
                    }
                }
            }
        }
    }
    if(continu)
    {
        dfg++;
        vector<int> mem;
		
        for(unsigned int i=0; i< afaire.size(); i++)
        {
            mem.push_back(afaire[i]);
            for(unsigned int o=0; o<2; o++)
            {
                hist.push_back(edge[afaire[i]].f[o]);
                for(unsigned int ui=0; ui<3; ui++)
                {
                    int fre=face[edge[afaire[i]].f[o]].v[ui];
                    hist.push_back(fre);
                    hist.push_back(vertex[fre].f);
                    fre=face[edge[afaire[i]].f[o]].e[ui];
                    hist.push_back(fre);
                    for(unsigned int y=0; y<2; y++)
                    {
                        hist.push_back(edge[fre].v[y]);
                        hist.push_back(edge[fre].f[y]);
                    }
                }
            }
        }
		
        int x=0;
        check=true;
        bool ghj=false;
        vector<int> seq;
        while(x<afaire.size())
        {
            while(x<afaire.size() && !problem(edge[afaire[x]].f[0], edge[afaire[x]].f[1]))
            {
                Flip(afaire[x]);
                seq.push_back(afaire[x]);
                x++;
            }
			
            if(x<afaire.size())
            {
                ghj=true;
                int y=x;
                while(y<afaire.size() && problem(edge[afaire[y]].f[0], edge[afaire[y]].f[1]))
                    y++;
                if(y==afaire.size())
                {
                    int y=x+1;
                    bool cont=true;
                    int tempc;
                    while(y<=afaire.size() && cont)
                    {
                        if(y<afaire.size())
                            tempc=troisieme(afaire[y-1], afaire[y], block);
                        else
                            tempc=troisieme(afaire[y-1], ed, block);
                        if(tempc!=-1)
                        {
                            Flip(tempc);
                            seq.push_back(tempc);
                            for(int h=y-1; h>x-1; h--)
                            {
                                if(!problem(edge[afaire[h]].f[0], edge[afaire[h]].f[1]))
                                {
                                    Flip(afaire[h]);
                                    seq.push_back(afaire[h]);
                                }
                            }
                            cont=false;
                        }
                        else
                        {
                            y++;
                        }
                    }
					
                    if(y>afaire.size())
                    {
                        check=false;
                    }
                }
                else
                {
                    for(int h=y; h>x-1; h--)
                    {
                        Flip(afaire[h]);
                        seq.push_back(afaire[h]);
                    }
                }
				
                afaire.clear();
                f0=edge[u].f[1];
                continu=true;
                temp=u;
                fg=f0;
                x=0;
                while(continu && temp!=ed)
                {
                    temp=EdgeCorresponding(fg, center, temp);
                    if(temp!=ed)
                    {
                        if(!block[temp])
                        {
                            continu=false;
                        }
                        else
                        {
                            afaire.push_back(temp);
                            if(edge[temp].f[0]==fg)
                            {
                                fg=edge[temp].f[1];
                            }
                            else
                            {
                                if(edge[temp].f[1]==fg)
                                    fg=edge[temp].f[0];
                                else
                                    fg=edge[temp].f[0];
                            }
                        }
                    }
                }
            }
        }
		
		
        int cas;
        if(afaire.empty())
            cas=DifferentFaceEdge(fg, u, ed);
        else
            cas=seq[seq.size()-1];
        if(check && block[cas] && ((edge[cas].v[0]==edge[u].v[(tempu+1)%2] && edge[cas].v[1]==edge[ed].v[(temped+1)%2]) || (edge[cas].v[1]==edge[u].v[(tempu+1)%2] && edge[cas].v[0]==edge[ed].v[(temped+1)%2])))
        {
            for(unsigned int o=0; o<seq.size(); o++)
            {
                comb.push_back(seq[o]);
            }
            return cas;
        }
        FlipR(hist);
        bbb=false;
    }
    return -1;
}

//_________________________________________________________________________________________________________________________________

bool Mesh::buildTriangle(const int& t, vector<bool>& block, vector<int>& FaceOk, vector<int>& vertexR)
//Construction d'une face par fermeture d'une face ou par une face avec un sommet incident à aucune arête bloquée
{
	//std::cout << " [ buildTriangle]" << std::endl;
    int temp;
    int e[3];

    for(int s=0; s<3; s++)
    {
        if(edgeT[faceT[t].e[s]].f[0]==t)
            temp=1;
        else
            temp=0;

        int k=FaceOk[edgeT[faceT[t].e[s]].f[temp]];
        if(k!=-1)
        {
            if((edge[face[k].e[0]].v[0]==vertexR[edgeT[faceT[t].e[s]].v[0]] && edge[face[k].e[0]].v[1]==vertexR[edgeT[faceT[t].e[s]].v[1]]) ||
               (edge[face[k].e[0]].v[0]==vertexR[edgeT[faceT[t].e[s]].v[1]] && edge[face[k].e[0]].v[1]==vertexR[edgeT[faceT[t].e[s]].v[0]]))
            {
                e[s]=face[k].e[0];
            }
            else
            {
                if((edge[face[k].e[1]].v[0]==vertexR[edgeT[faceT[t].e[s]].v[0]] && edge[face[k].e[1]].v[1]==vertexR[edgeT[faceT[t].e[s]].v[1]]) ||
                   (edge[face[k].e[1]].v[0]==vertexR[edgeT[faceT[t].e[s]].v[1]] && edge[face[k].e[1]].v[1]==vertexR[edgeT[faceT[t].e[s]].v[0]]))
                    e[s]=face[k].e[1];
                else
                    e[s]=face[k].e[2];
            }
        }
        else
            e[s]=-1;
    }
//std::cout << " 1" << std::endl;
    if(e[0]==-1)
    {
        int temp=e[0];
        e[0]=e[2];
        e[2]=temp;
    }
    if(e[1]==-1)
    {
        int temp=e[1];
        e[1]=e[2];
        e[2]=temp;
    }
    else
    {
        int temp=e[1];
        e[1]=e[0];
        e[0]=temp;
    }

//std::cout << " 2" << std::endl;
    if(e[1]==-1)
    {
		//std::cout << " A" << std::endl;
        vector<int> result;
        int newvertex=-1;
        if(vertexR[faceT[t].v[0]]!=edge[e[0]].v[0] && vertexR[faceT[t].v[0]]!=edge[e[0]].v[1])
        {
            newvertex=chooseVertex(vertexR, faceT[t].v[0], e[0], block);
        }
        else
        {
            if(vertexR[faceT[t].v[1]]!=edge[e[0]].v[0] && vertexR[faceT[t].v[1]]!=edge[e[0]].v[1])
            {
                newvertex=chooseVertex(vertexR, faceT[t].v[1], e[0], block);
            }
            else
            {
                if(vertexR[faceT[t].v[2]]!=edge[e[0]].v[0] && vertexR[faceT[t].v[2]]!=edge[e[0]].v[1])
                {
                    newvertex=chooseVertex(vertexR, faceT[t].v[2], e[0], block);
                }
            }
        }

        if(true)
        {

			//Cas d'une création de facette vers un sommet non utilisé
            int temp=findFaceVertex(edge[e[0]].f[0], vertexR[faceT[t].v[0]]);
            if(temp>-1)
            {
                if(face[edge[e[0]].f[0]].v[temp]==vertexR[faceT[t].v[0]] && face[edge[e[0]].f[0]].v[(temp+1)%3]==vertexR[faceT[t].v[1]] && face[edge[e[0]].f[0]].v[(temp+2)%3]==vertexR[faceT[t].v[2]])
                {
                    if(block[face[edge[e[0]].f[0]].e[0]] || block[face[edge[e[0]].f[0]].e[1]] || block[face[edge[e[0]].f[0]].e[2]])
                    {
                        for(int g=0; g<3; g++)
                        {
                            block[face[edge[e[0]].f[0]].e[g]]=false;
                        }
                        FaceOk[t]=edge[e[0]].f[0];
                        return true;
                    }
                }
            }

            temp=findFaceVertex(edge[e[0]].f[1], vertexR[faceT[t].v[0]]);
            if(temp>-1)
            {
                if(face[edge[e[0]].f[1]].v[temp]==vertexR[faceT[t].v[0]] && face[edge[e[0]].f[1]].v[(temp+1)%3]==vertexR[faceT[t].v[1]] && face[edge[e[0]].f[1]].v[(temp+2)%3]==vertexR[faceT[t].v[2]])
                {
                    if(block[face[edge[e[0]].f[1]].e[0]] || block[face[edge[e[0]].f[1]].e[1]] || block[face[edge[e[0]].f[1]].e[2]])
                    {
                        for(int g=0; g<3; g++)
                        {
                            block[face[edge[e[0]].f[1]].e[g]]=false;
                        }
                        FaceOk[t]=edge[e[0]].f[1];
                        return true;
                    }
                }
            }

            for(int bene=0; bene<2; bene++)
            {
                if(edgeExistNoBlock(edge[e[0]].v[bene], newvertex, result, block))
                {
                    for(int s=0; s<result.size(); s++)
                    {
                        block[result[s]]=false;
						//std::cout << " close1" << std::endl;
                        int r=CloseTriangle(result[s], e[0], block);
                        if(r>-1)
                        {
                            block[r]=false;

                            if((face[edge[r].f[0]].e[0]==e[0] || face[edge[r].f[0]].e[1]==e[0] || face[edge[r].f[0]].e[2]==e[0]) &&
                               (face[edge[r].f[0]].e[0]==result[s] || face[edge[r].f[0]].e[1]==result[s] || face[edge[r].f[0]].e[2]==result[s]))
                                 FaceOk[t]=edge[r].f[0];
                            else
                            {
                                if((face[edge[r].f[1]].e[0]==e[0] || face[edge[r].f[1]].e[1]==e[0] || face[edge[r].f[1]].e[2]==e[0]) &&
                                   (face[edge[r].f[1]].e[0]==result[s] || face[edge[r].f[1]].e[1]==result[s] || face[edge[r].f[1]].e[2]==result[s]))
                                    FaceOk[t]=edge[r].f[1];
                                else
                                {
                                    return false;
                                }
                            }
                            return true;
                        }
                        block[result[s]]=true;
                    }
                }
            }

            if(findFaceVertex(edge[e[0]].f[0],newvertex)==-1 && findFaceVertex(edge[e[0]].f[1],newvertex)==-1)
            {
				//std::cout << " Build1" << std::endl;
                int y=BuildEdgeStarFromEdgeWithClose(edge[e[0]].v[0], newvertex, e[0], block, t);
                if(y>-1)
                {
                    FaceOk[t]=y;
                    return true;
                }
            }
        }
    }
    else
    {
		//std::cout << " 5" << std::endl;
        if(e[2]==-1)
        {
			//std::cout << " B" << std::endl;
			// Cas de fermeture d'une facette
			//std::cout << " Close2" << std::endl;
            int u=CloseTriangle(e[0], e[1], block);
            if(u>-1)
            {
                block[u]=false;
                if((face[edge[u].f[0]].e[0]==e[0] || face[edge[u].f[0]].e[1]==e[0] || face[edge[u].f[0]].e[2]==e[0]) &&
                   (face[edge[u].f[0]].e[0]==e[1] || face[edge[u].f[0]].e[1]==e[1] || face[edge[u].f[0]].e[2]==e[1]))
                    FaceOk[t]=edge[u].f[0];
                else
                {
                    if((face[edge[u].f[1]].e[0]==e[0] || face[edge[u].f[1]].e[1]==e[0] || face[edge[u].f[1]].e[2]==e[0]) &&
                       (face[edge[u].f[1]].e[0]==e[1] || face[edge[u].f[1]].e[1]==e[1] || face[edge[u].f[1]].e[2]==e[1]))
                        FaceOk[t]=edge[u].f[1];
                    else
                    {
                        return false;
                    }
                }
                return true;
            }
        }
    }
    return false;
}

//_______________________________________________________________________________________________________________________________________________________________

int Mesh::CloseTriangleWithOrientation(const int& ed, const int& u, const int& t, const vector<bool>& block)
//Fermeture d'une face avec vérification de l'orientation
{
    vector<int> hist;
    vector<int> afaire;
    int center;
    int tempu, temped;
    if(edge[u].v[0]==edge[ed].v[0])
    {
        tempu=0;
        temped=0;
        center=edge[u].v[0];
    }
    else
    {
        if(edge[u].v[0]==edge[ed].v[1])
        {
            tempu=0;
            temped=1;
            center=edge[u].v[0];
        }
        else
        {
            if(edge[u].v[1]==edge[ed].v[0])
            {
                tempu=1;
                temped=0;
                center=edge[u].v[1];
            }
            else
            {
                if(edge[u].v[1]==edge[ed].v[1])
                {
                    tempu=1;
                    temped=1;
                    center=edge[u].v[1];
                }
            }
        }
    }
    int f0=edge[u].f[0];
    bool continu=true;
    int temp=u;
    int fg=f0;
    while(continu && temp!=ed)
    {
        temp=EdgeCorresponding(fg, center, temp);
        if(temp!=ed)
        {
            if(!block[temp])
            {
                continu=false;
            }
            else
            {
                afaire.push_back(temp);
                if(edge[temp].f[0]==fg)
                {
                    fg=edge[temp].f[1];
                }
                else
                {
                    if(edge[temp].f[1]==fg)
                    {
                        fg=edge[temp].f[0];
                    }
                    else
                    {
                        fg=edge[temp].f[0];
                    }
                }
            }
        }
    }
    if(continu)
    {
        for(unsigned int i=0; i< afaire.size(); i++)
        {
            for(unsigned int o=0; o<2; o++)
            {
                hist.push_back(edge[afaire[i]].f[o]);
                for(unsigned int ui=0; ui<3; ui++)
                {
                    int fre=face[edge[afaire[i]].f[o]].v[ui];
                    hist.push_back(fre);
                    hist.push_back(vertex[fre].f);
                    fre=face[edge[afaire[i]].f[o]].e[ui];
                    hist.push_back(fre);
                    for(unsigned int y=0; y<2; y++)
                    {
                        hist.push_back(edge[fre].v[y]);
                        hist.push_back(edge[fre].f[y]);
                    }
                }
            }
        }
        bool check=true;
        for(unsigned int x=0; x<afaire.size(); x++)
        {
            if(problem(edge[afaire[x]].f[0], edge[afaire[x]].f[1]))
            {
                int bene=x+1;
                while(bene<afaire.size() && problem(edge[afaire[bene]].f[0], edge[afaire[bene]].f[1]))
                    bene++;
                if(bene==afaire.size())
                {
                    check=false;
                    x=afaire.size()-1;
                }
                else
                {
                    int libe=x;
                    int temp;
                    while(libe<bene)
                    {
                        temp=afaire[libe];
                        afaire[libe]=afaire[bene];
                        afaire[bene]=temp;
                        libe++;
                        bene--;
                    }
                }
            }
            Flip(afaire[x]);
        }
        int cas;
        if(afaire.empty())
        {
            cas=DifferentFaceEdge(f0, u, ed);
        }
        else
            cas=afaire[afaire.size()-1];
        if(check && CheckOrientation(fg, t) && ((edge[cas].v[0]==edge[u].v[(tempu+1)%2] && edge[cas].v[1]==edge[ed].v[(temped+1)%2]) || (edge[cas].v[1]==edge[u].v[(tempu+1)%2] && edge[cas].v[0]==edge[ed].v[(temped+1)%2])))
        {
            for(unsigned int o=0; o<afaire.size(); o++)
                comb.push_back(afaire[o]);
            return cas;
        }
        FlipR(hist);
    }
    else
    {
        FlipR(hist);
    }

    afaire.clear();
    hist.clear();
    f0=edge[u].f[1];
    continu=true;
    temp=u;
    fg=f0;

    while(continu && temp!=ed)
    {
        temp=EdgeCorresponding(fg, center, temp);
        if(temp!=ed)
        {
            if(!block[temp])
            {
                continu=false;
            }
            else
            {
                afaire.push_back(temp);
                if(edge[temp].f[0]==fg)
                {
                    fg=edge[temp].f[1];
                }
                else
                {
                    if(edge[temp].f[1]==fg)
                    {
                        fg=edge[temp].f[0];
                    }
                    else
                    {
                        fg=edge[temp].f[0];
                    }
                }
            }
        }
    }

    if(continu)
    {
        for(unsigned int i=0; i< afaire.size(); i++)
        {
            for(unsigned int o=0; o<2; o++)
            {
                hist.push_back(edge[afaire[i]].f[o]);
                for(unsigned int ui=0; ui<3; ui++)
                {
                    int fre=face[edge[afaire[i]].f[o]].v[ui];
                    hist.push_back(fre);
                    hist.push_back(vertex[fre].f);
                    fre=face[edge[afaire[i]].f[o]].e[ui];
                    hist.push_back(fre);
                    for(unsigned int y=0; y<2; y++)
                    {
                        hist.push_back(edge[fre].v[y]);
                        hist.push_back(edge[fre].f[y]);
                    }
                }
            }
        }

        bool check=true;
        for(unsigned int x=0; x<afaire.size(); x++)
        {
            if(problem(edge[afaire[x]].f[0], edge[afaire[x]].f[1]))
            {
                int bene=x+1;
                while(bene<afaire.size() && problem(edge[afaire[bene]].f[0], edge[afaire[bene]].f[1]))
                    bene++;
                if(bene==afaire.size())
                {
                    check=false;
                    x=afaire.size()-1;
                }
                else
                {
                    int libe=x;
                    int temp;
                    while(libe<bene)
                    {
                        temp=afaire[libe];
                        afaire[libe]=afaire[bene];
                        afaire[bene]=temp;
                        libe++;
                        bene--;
                    }
                }
            }
            Flip(afaire[x]);
        }

        int cas;
        if(afaire.empty())
            cas=DifferentFaceEdge(fg, u, ed);
        else
        {
            cas=afaire[afaire.size()-1];
        }
        if(check && CheckOrientation(fg, t) && ((edge[cas].v[0]==edge[u].v[(tempu+1)%2] && edge[cas].v[1]==edge[ed].v[(temped+1)%2]) || (edge[cas].v[1]==edge[u].v[(tempu+1)%2] && edge[cas].v[0]==edge[ed].v[(temped+1)%2])))
        {
            for(unsigned int o=0; o<afaire.size(); o++)
                comb.push_back(afaire[o]);
            return cas;
        }
        FlipR(hist);
    }
    else
    {
        FlipR(hist);
    }
    return -1;
}

//_________________________________________________________________________________________________________________________________


int Mesh::init(const int& t, vector<bool>& block, const int& color)
//Construction de la première face de la croissance de région
{
	//std::cout << "init -> " << t << std::endl;
    vector<int> result;
    int e1=edgeExist(vertexR[faceT[t].v[0]], vertexR[faceT[t].v[1]], result);
    if(e1>0)
    {
		//std::cout << "	1" << std::endl;
        for(int a1=0; a1<result.size(); a1++)
        {
            int temp=findFaceVertex(edge[result[a1]].f[0], vertexR[faceT[t].v[0]]);
            if(temp>-1)
            {
				//std::cout << "	1.1" << std::endl;
                if(face[edge[result[a1]].f[0]].v[temp]==vertexR[faceT[t].v[0]] && face[edge[result[a1]].f[0]].v[(temp+1)%3]==vertexR[faceT[t].v[1]] && face[edge[result[a1]].f[0]].v[(temp+2)%3]==vertexR[faceT[t].v[2]])
                {
                    for(int g=0; g<3; g++)
                    {
                        block[face[edge[result[a1]].f[0]].e[g]]=false;
                    }
                    return edge[result[a1]].f[0];
                }
            }
            temp=findFaceVertex(edge[result[a1]].f[1], vertexR[face[t].v[0]]);
            if(temp>-1)
            {
//std::cout << "	1.2" << std::endl;
                if(face[edge[result[a1]].f[1]].v[temp]==vertexR[faceT[t].v[0]] && face[edge[result[a1]].f[1]].v[(temp+1)%3]==vertexR[faceT[t].v[1]] && face[edge[result[a1]].f[1]].v[(temp+2)%3]==vertexR[faceT[t].v[2]])
                {
                    for(int g=0; g<3; g++)
                    {
                        block[face[edge[result[a1]].f[1]].e[g]]=false;
                    }
                    return edge[result[a1]].f[1];
                }
            }
        }
        e1=result[0];
        block[e1]=false;
    }
    else
    {
		//std::cout << "	2" << std::endl;
        e1=BuildEdgeStar(vertexR[faceT[t].v[0]], vertexR[faceT[t].v[1]], block, false, false);
        if(e1>-1)
        {
			//std::cout << "	2.1" << std::endl;
            block[e1]=false;
        }
        else
		{
			//std::cout << "	-1" << std::endl;
            return -1;
		}
    }
    int e2=edgeExist(vertexR[faceT[t].v[1]], vertexR[faceT[t].v[2]]);
    if(e2>-1)
    {
        block[e2]=false;
//std::cout << "	3" << std::endl;
        int r=CloseTriangleWithOrientation(e1, e2, t, block);
        if(r>-1)
        {
            block[r]=false;
//std::cout << "	3.1" << std::endl;
            if((face[edge[r].f[0]].e[0]==e1 || face[edge[r].f[0]].e[1]==e1 || face[edge[r].f[0]].e[2]==e1) &&
               (face[edge[r].f[0]].e[0]==e2 || face[edge[r].f[0]].e[1]==e2 || face[edge[r].f[0]].e[2]==e2))
			{
				//std::cout << "	3.11" << std::endl;
                return edge[r].f[0];
			}
            else
            {
                if((face[edge[r].f[1]].e[0]==e1 || face[edge[r].f[1]].e[1]==e1 || face[edge[r].f[1]].e[2]==e1) &&
                   (face[edge[r].f[1]].e[0]==e2 || face[edge[r].f[1]].e[1]==e2 || face[edge[r].f[1]].e[2]==e2))
				{
					//std::cout << "	3.12" << std::endl;
                    return edge[r].f[1];
				}
                else
                {
					//std::cout << "	-1" << std::endl;
                    return -1;
                }
            }
        }
        block[e2]=true;
    }
    else
    {
		//std::cout << "	4" << std::endl;
        int h=BuildEdgeStar(vertexR[faceT[t].v[1]], vertexR[faceT[t].v[2]], block, false, false);
        if(h>-1)
        {
			//std::cout << "	4.1" << std::endl;
            block[h]=false;
            int r=CloseTriangleWithOrientation(e1, h, t, block);
            if(r>-1)
            {
				//std::cout << "	4.2" << std::endl;
                block[r]=false;
                if((face[edge[r].f[0]].e[0]==e1 || face[edge[r].f[0]].e[1]==e1 || face[edge[r].f[0]].e[2]==e1) &&
                   (face[edge[r].f[0]].e[0]==h || face[edge[r].f[0]].e[1]==h || face[edge[r].f[0]].e[2]==h))
                    return edge[r].f[0];
                else
                {
                    if((face[edge[r].f[1]].e[0]==e1 || face[edge[r].f[1]].e[1]==e1 || face[edge[r].f[1]].e[2]==e1) &&
                       (face[edge[r].f[1]].e[0]==h || face[edge[r].f[1]].e[1]==h || face[edge[r].f[1]].e[2]==h))
                        return edge[r].f[1];
                    else
                    {
						//std::cout << "	-1" << std::endl;
                        return -1;
                    }
                }
            }
            block[h]=true;
        }
    }
    block[e1]=true;
    return -1;
}

//_____________________________________________________________________________________________________________________________________________


bool Mesh::NoConstructible(const int& t, vector<int>& FaceOk, vector<bool>& block)
{
    unsigned int count=0;
    for(int i=0; i<3; i++)
    {
        if(edgeT[faceT[t].e[i]].f[0]==t)
        {

            if(edgeT[faceT[t].e[i]].f[1]!=-1 && FaceOk[edgeT[faceT[t].e[i]].f[1]]!=-1)
            {
                count++;
            }
        }
        else
        {
            if(edgeT[faceT[t].e[i]].f[0]!=-1 && FaceOk[edgeT[faceT[t].e[i]].f[0]]!=-1)
            {
                count++;
            }
        }
    }
    if(count==2)
        return false;

    count=0;
    for(int i=0; i<3; i++)
    {
        bool nb=true;
        int v1=faceT[t].v[i];
        count=0;
        int f=vertexT[v1].f;
        int f0=f;
        do
        {
            if(FaceOk[f]==-1)
            {
                int e;
                if(faceT[f].v[0]==v1)
                {
                    e=faceT[f].e[1];
                }
                else
                {
                    if(faceT[f].v[1]==v1)
                    {
                        e=faceT[f].e[2];
                    }
                    else
                    {
                        if(faceT[f].v[2]==v1)
                        {
                            e=faceT[f].e[0];
                        }
                        else
                        {
                            return false;
                        }
                    }
                }
                if(edgeT[e].f[0]==f)
                    f=edgeT[e].f[1];
                else
                    f=edgeT[e].f[0];
                count++;
            }
            else
                nb=false;
        } while(nb && f!=f0 && count<edge.size());

        if(nb)
        {
            return false;
        }
    }
    return true;
}

//___________________________________________________________________________________________________

bool Mesh::buildTriangleWithModifyGenus(const int& t, vector<bool>& block, vector<int>& FaceOk, vector<int>& vertexR) //Cas d'une création de facette vers un sommet déjà utilisé
//Construction d'une face à partir d'une unique arête bloquée et trois sommets incidents à des arêtes bloquées
{
    int temp;
    int e[3];
	
    for(int s=0; s<3; s++)
    {
        if(edgeT[faceT[t].e[s]].f[0]==t)
            temp=1;
        else
            temp=0;
		
        int k=FaceOk[edgeT[faceT[t].e[s]].f[temp]];
        if(k!=-1)
        {
			
            if((edge[face[k].e[0]].v[0]==vertexR[edgeT[faceT[t].e[s]].v[0]] && edge[face[k].e[0]].v[1]==vertexR[edgeT[faceT[t].e[s]].v[1]]) ||
               (edge[face[k].e[0]].v[0]==vertexR[edgeT[faceT[t].e[s]].v[1]] && edge[face[k].e[0]].v[1]==vertexR[edgeT[faceT[t].e[s]].v[0]]))
            {
                e[s]=face[k].e[0];
            }
            else
            {
                if((edge[face[k].e[1]].v[0]==vertexR[edgeT[faceT[t].e[s]].v[0]] && edge[face[k].e[1]].v[1]==vertexR[edgeT[faceT[t].e[s]].v[1]]) ||
                   (edge[face[k].e[1]].v[0]==vertexR[edgeT[faceT[t].e[s]].v[1]] && edge[face[k].e[1]].v[1]==vertexR[edgeT[faceT[t].e[s]].v[0]]))
                    e[s]=face[k].e[1];
                else
                    e[s]=face[k].e[2];
            }
        }
        else
            e[s]=-1;
    }
	
    if(e[0]==-1)
    {
        int temp=e[0];
        e[0]=e[2];
        e[2]=temp;
    }
    if(e[1]==-1)
    {
        int temp=e[1];
        e[1]=e[2];
        e[2]=temp;
    }
    else
    {
        int temp=e[1];
        e[1]=e[0];
        e[0]=temp;
    }
	
    if(e[1]==-1)
    {
        vector<int> result;
        int newvertex=-1;
        if(vertexR[faceT[t].v[0]]!=edge[e[0]].v[0] && vertexR[faceT[t].v[0]]!=edge[e[0]].v[1])
        {
            newvertex=vertexR[faceT[t].v[0]];
        }
        else
        {
            if(vertexR[faceT[t].v[1]]!=edge[e[0]].v[0] && vertexR[faceT[t].v[1]]!=edge[e[0]].v[1])
            {
                newvertex=vertexR[faceT[t].v[1]];
            }
            else
            {
                if(vertexR[faceT[t].v[2]]!=edge[e[0]].v[0] && vertexR[faceT[t].v[2]]!=edge[e[0]].v[1])
                {
                    newvertex=vertexR[faceT[t].v[2]];
                }
            }
        }
		
        for(int y=0; y<2; y++)
        {
            result.clear();
            if(edgeExistNoBlock(edge[e[0]].v[y], newvertex, result, block))
            {
                int s=0;
                while(s<result.size())
                {
                    if(composante_connexe(result[s], block))
                    {
                        block[result[s]]=false;
                        int r=CloseTriangle(result[s], e[0], block);
                        if(r>-1)
                        {
                            block[r]=false;
							
                            if((face[edge[r].f[0]].e[0]==e[0] || face[edge[r].f[0]].e[1]==e[0] || face[edge[r].f[0]].e[2]==e[0]) &&
                               (face[edge[r].f[0]].e[0]==result[s] || face[edge[r].f[0]].e[1]==result[s] || face[edge[r].f[0]].e[2]==result[s]))
                                FaceOk[t]=edge[r].f[0];
                            else
                            {
                                if((face[edge[r].f[1]].e[0]==e[0] || face[edge[r].f[1]].e[1]==e[0] || face[edge[r].f[1]].e[2]==e[0]) &&
                                   (face[edge[r].f[1]].e[0]==result[s] || face[edge[r].f[1]].e[1]==result[s] || face[edge[r].f[1]].e[2]==result[s]))
                                    FaceOk[t]=edge[r].f[1];
                                else
                                {
                                    std::cout << " GROS SOUCIS [cas 3] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
                                    return false;
                                }
                            }
                            return true;
                        }
                        else
                        {
                            block[result[s]]=true;
                        }
                    }
                    s++;
                }
            }
        }
		
        if(findFaceVertex(edge[e[0]].f[0],newvertex)==-1 && findFaceVertex(edge[e[0]].f[1],newvertex)==-1)
        {
            int e1=BuildEdgeStarFromEdge(edge[e[0]].v[0], newvertex, e[0], block, t, true);
            if(e1>-1)
            {
                block[e1]=false;
                int r=CloseTriangle(e1, e[0], block);
                if(r>-1)
                {
                    block[r]=false;
                    if((face[edge[r].f[0]].e[0]==e[0] || face[edge[r].f[0]].e[1]==e[0] || face[edge[r].f[0]].e[2]==e[0]) &&
                       (face[edge[r].f[0]].e[0]==e1 || face[edge[r].f[0]].e[1]==e1 || face[edge[r].f[0]].e[2]==e1))
                        FaceOk[t]=edge[r].f[0];
                    else
                    {
                        if((face[edge[r].f[1]].e[0]==e[0] || face[edge[r].f[1]].e[1]==e[0] || face[edge[r].f[1]].e[2]==e[0]) &&
                           (face[edge[r].f[1]].e[0]==e1 || face[edge[r].f[1]].e[1]==e1 || face[edge[r].f[1]].e[2]==e1))
                            FaceOk[t]=edge[r].f[1];
                        else
                        {
                            return false;
                        }
                    }
                    return true;
                }
                else
                {
                    block[e1]=true;
                }
            }
        }
    }
    return false;
}

//_________________________________________________________________________________________________________________________________


void Mesh::buildEdgeWithCondition(bool simply)
//Algo principal pour déterminer une séquence de bascules d'arêtes entre les deux triangulations
{
    if(!edge.empty() && edgeT.size()==edge.size())
    {
		comb.clear();
		reverse=false;

		std::deque<int> remains;

		vector<bool> block;
		block.resize(edge.size(), true);

		vector<int> FaceOk;
		FaceOk.resize(face.size(), -1);
		
		int tt=0;
        int t=0;
		while(tt==0 && t<face.size())
		{
			int k=init(t, block, 1);
			if(k!=-1)
			{
				FaceOk[t]=k;
				tt++;
				for(int i=0; i<3; i++)
				{
					if(edgeT[faceT[t].e[i]].f[0]==t)
						remains.push_back(edgeT[faceT[t].e[i]].f[1]);
					else
						remains.push_back(edgeT[faceT[t].e[i]].f[0]);
				}
			}
		}
		
		int coun=0;
		std::cout << "	[ " << coun << "%";
		coun=coun+10;
		int i=0;
		int genus=0;
		int last_unconstructible=-1;

		while(!remains.empty() && tt!=t)
		{
			//std::cout << "	A" << std::endl;
			t=tt;
			while(!remains.empty())
			{
				//std::cout << "								facette : " << remains[0] << "  [ " <<  faceT.size()-tt << "  ] " << std::endl;
				if(FaceOk[remains[0]]!=-1 || NoConstructible(remains[0], FaceOk, block))
				{
					if(FaceOk[remains[0]]==-1)
						last_unconstructible=remains[0];
					remains.pop_front();
				}
				else
				{
					if(buildTriangle(remains[0], block, FaceOk, vertexR))
					{
						for(int i=0; i<3; i++)
						{
							if(FaceOk[edgeT[faceT[remains[0]].e[i]].f[0]]==-1)
								remains.push_back(edgeT[faceT[remains[0]].e[i]].f[0]);
							if(FaceOk[edgeT[faceT[remains[0]].e[i]].f[1]]==-1)
								remains.push_back(edgeT[faceT[remains[0]].e[i]].f[1]);
						}
						remains.pop_front();
						tt++;
						
						
						if(tt%((face.size()/10)+2)==0)
						{
							std::cout << " - " << coun << "%";
							coun=coun+10;
						}
						
					}
					else
					{
						int y=remains[0];
						remains.pop_front();
						remains.push_back(y);
					}
				}
			}
			if(tt<faceT.size()-1)
			{
				//std::cout << std::endl << " * Cas 3 *" << std::endl << std::endl;
				
				int imax=i;
				i=(i+1)%faceT.size();
				bool ok=true;
				if(last_unconstructible==-1 || FaceOk[last_unconstructible]!=-1)
				{
					while(i!=imax && (FaceOk[i]!=-1 || (FaceOk[edge[face[i].e[0]].f[0]]==-1 && FaceOk[edge[face[i].e[0]].f[1]]==-1 && FaceOk[edge[face[i].e[1]].f[0]]==-1 && FaceOk[edge[face[i].e[1]].f[1]]==-1 && FaceOk[edge[face[i].e[2]].f[0]]==-1 && FaceOk[edge[face[i].e[2]].f[1]]==-1)))
						i=(i+1)%faceT.size();
					if(i==imax)
						ok=false;
					else
					{
						last_unconstructible=i;
						i=(i+1)%faceT.size();
					}
				}
				//std::cout << " i -> " << i << std::endl;
				while(ok && !buildTriangleWithModifyGenus(last_unconstructible, block, FaceOk, vertexR))
				{
					while(i!=(i+1)%faceT.size() && (FaceOk[i]!=-1 || (FaceOk[edge[face[i].e[0]].f[0]]==-1 && FaceOk[edge[face[i].e[0]].f[1]]==-1 && FaceOk[edge[face[i].e[1]].f[0]]==-1 && FaceOk[edge[face[i].e[1]].f[1]]==-1 && FaceOk[edge[face[i].e[2]].f[0]]==-1 && FaceOk[edge[face[i].e[2]].f[1]]==-1)))

						i=(i+1)%faceT.size();
					if(i==imax)
					{
						ok=false;
					}
					else
					{
						last_unconstructible=i;
						i=(i+1)%faceT.size();
					}
					//std::cout << " i -> " << i << std::endl;
				}
				if(ok)
				{
					tt++;
					genus++;
					for(int w=0; w<3; w++)
					{
						if(FaceOk[edgeT[faceT[last_unconstructible].e[w]].f[0]]==-1)
						{
							remains.push_back(edgeT[faceT[last_unconstructible].e[w]].f[0]);
						}
						if(FaceOk[edgeT[faceT[last_unconstructible].e[w]].f[1]]==-1)
						{
							remains.push_back(edgeT[faceT[last_unconstructible].e[w]].f[1]);
						}
					}
				}
			}
		}
		
		if(faceT.size()<tt+3)
		{
			std::cout << " - 100% ]"<< std::endl;
			std::cout << "	Genus : " << genus/2 << std::endl;
			std::cout << "	Number of remaining triangles : " << faceT.size()-tt << std::endl;
			std::cout << "	Number of edge flips : " << comb.size() << std::endl;
		}
		else
		{
			std::cout << " - failure ]"<< std::endl;
			std::cout << "	Number of remaining triangles : " << faceT.size()-tt << std::endl;
			std::cout << "	Number of edge flips : " << comb.size() << std::endl;
		}
	}
}

//________________________________________________________________________

int Mesh::chooseVertex(vector<int>& vertexR, const int& v, const int & e, const vector<bool>& block)
{
    if(vertexR[v]==-1)
    {
        vertexR[v]=v;
        return v;
    }
    else
        return vertexR[v];
}

//________________________________________________________________________________________________________________________________________________________________

void Mesh::generateFile()
//Genère un fichier .txt contenant une liste des arêtes basculées
{
    std::ofstream temp;
    temp.open ("Sequence_Of_Flips.txt");
    temp << "#Flips: " << comb.size() << std::endl;
    for(int i=comb.size()-1; i>-1; i--)
    {
        Flip(comb[i]);
    }
    for(int i=0; i<comb.size(); i++)
    {
        temp << edge[comb[i]].v[0] << " " << edge[comb[i]].v[1] << std::endl;
        Flip(comb[i]);
    }
    temp.close();
}

//________________________________________________________________________________________________________________________________________________________________


/*/*******************/
/*/* Simplification **/
/*/*******************/

//________________________________________________________________________________________________________________________________________________________________

int Mesh::simplifysequence()
//Algo principale pour réduire une séquence de bascules d'arêtes
{
    int tailleAvant=comb.size();
    gain=0;
    simplifyLastFlip(gain,  0);
    moy=(100*gain)/tailleAvant;
    if(tailleAvant>0)
    {
        std::cout << "	Number of edge flips after reduction : " << comb.size() << "  (gain : " << gain << "  i.e  " << moy << "%)" << std::endl;
    }
}

//________________________________________________________________________________________________________________________________________________________________


bool Mesh::commute(const int& a, const int& b)
// return true if a commute with b, if ot false
{
    if(a==-1 || b==-1)
        return true;
    int A[2], B[2];
    A[0]=edge[a].v[0];
    A[1]=edge[a].v[1];
    B[0]=edge[b].v[0];
    B[1]=edge[b].v[1];
    Flip(a);
    Flip(b);
    Flip(a);
    Flip(b);
    if((edge[a].v[0]==A[0] || edge[a].v[1]==A[0]) && (edge[a].v[0]==A[1] || edge[a].v[1]==A[1]) && (edge[b].v[0]==B[0] || edge[b].v[1]==B[0]) && (edge[b].v[0]==B[1] || edge[b].v[1]==B[1]))
    {
        return true;
    }
    Flip(b);
    Flip(a);
    Flip(b);
    Flip(a);
    return false;
}

//_________________________________________________________________________________________________________________________________________________________________


bool Mesh::simplificationNoCommutativity(const int& a, const int& b)
//Retroune Vrai si l'on peut simplifier une séquence par non commutativité
{
    int A[2], B[2];
    A[0]=edge[a].v[0];
    A[1]=edge[a].v[1];
    B[0]=edge[b].v[0];
    B[1]=edge[b].v[1];
        Flip(a);
        Flip(b);
        Flip(a);
        Flip(b);
        Flip(a);
    if((edge[b].v[0]==A[0] || edge[b].v[1]==A[0]) && (edge[b].v[0]==A[1] || edge[b].v[1]==A[1]) && (edge[a].v[0]==B[0] || edge[a].v[1]==B[0]) && (edge[a].v[0]==B[1] || edge[a].v[1]==B[1]))
    {
        Flip(a);
        Flip(b);
        Flip(a);
        Flip(b);
        Flip(a);
        return true;
    }
    Flip(a);
    Flip(b);
    Flip(a);
    Flip(b);
    Flip(a);
    return false;
}

//-------------------------------------------------------------------------------------------------------------------------------

int Mesh::Simplify(int i, int & first)
{
    int mmax=-1;
    if(i>0 && i<comb.size() && comb[i]!=-1)
    {
        int temp=comb[i];
        int other=i-1;
        while(other>-1 && comb[other]!=temp)
            other--;

        if(other>-1)// there a second flip i in comb
        {
            int j=i;
            int tamp;


            if(first>j-1)
            {
                for(int y=first-1; y>j-2; y--)
                    if(comb[y]!=-1)
                        Flip(comb[y]);
            }
            if(first<j-1)
            {
                for(int y=first; y<j-1; y++)
                    if(comb[y]!=-1)
                        Flip(comb[y]);
            }

            while(j>other+1 && commute(comb[j], comb[j-1]))
            {
                tamp=comb[j];
                comb[j]=comb[j-1];
                comb[j-1]=tamp;
                j--;
                if(comb[j-1]!=-1)
                    Flip(comb[j-1]);
            }
            first=j-1;

            if(j-1>other)
            {
                for(int u=first-1; u>other-1; u--)
                    if(comb[u]!=-1)
                        Flip(comb[u]);
                while(other+1<j && commute(comb[other+1], comb[other]))
                {
                    tamp=comb[other];
                    comb[other]=comb[other+1];
                    comb[other+1]=tamp;
                    other++;
                    if(comb[other-1]!=-1)
                        Flip(comb[other-1]);
                }

                first=other;
            }

            if(j==other+1)
            {
                comb[j]=-1;
                comb[j-1]=-1;
                for(unsigned int u=j+1; u<comb.size(); u++)
                    comb[u-2]=comb[u];
                comb.pop_back();
                comb.pop_back();

            }
            else
            {
                if(j==other+2 && simplificationNoCommutativity(comb[other], comb[j-1]))
                {
                    tamp=comb[j-2];
                    comb[j-2]=comb[j-1];
                    comb[j-1]=tamp;

                    int k=1;
                    for(unsigned int u=j; u<comb.size()-1; u++)
                    {
                        if(comb[u+1]==-1)
                        {
                            k++;
                        }
                        else
                        {
                            if(comb[u+k]==comb[j-1])
                                comb[u]=comb[j-2];
                            else
                            {
                                if(comb[u+k]==comb[j-2])
                                    comb[u]=comb[j-1];
                                else
                                    comb[u]=comb[u+k];
                            }
                        }
                    }
                    for(int b=0; b<k; b++)
                    {
                        comb.pop_back();
                    }
                }
                else
                {
                    return -1;
                }
            }
            mmax=first;
        }
    }
    return mmax;
}


//_______________________________________________________________________________________________________________________________________________________________________________________________________


void Mesh::simplifyLastFlip(int& gain, int i)
{
    int ret=i;
    int first=comb.size();
    int out, taille;
    while(ret<comb.size())
    {
        taille=comb.size();
        out=Simplify(ret, first);
        gain=gain+taille-comb.size();
        if(out==-1)
            ret++;
        else
        {
            ret=out;
        }
    }
    int k=0;
    for(unsigned int u=0; u<comb.size(); u++)
    {
        if(comb[u]==-1)
        {
             k=k+1;
        }
        else
        {
             if(k>0)
             {
                comb[u-k]=comb[u];
             }
         }
    }
    for(int b=0; b<k; b++)
    {
        comb.pop_back();
    }

    for(unsigned int s=first; s<comb.size(); s++)
            Flip(comb[s]);
}

//______________________________________________________________________________________________________________

void Mesh::SimplifySeq()
//Second algorithme indépendant pour réduire une séquence de bascules d'arêtes
{
    gain=0;
    int t=comb.size();
    int b=0;
    int t1=0;
    vector<vector<int> > seq;

    while(comb.size()-t1>5)
    {
        b++;
        t1=comb.size();
        seq.clear();
        seq.resize(comb.size());
		for(int i=0; i<comb.size(); i++)
        {
            int y=comb.size()-1-i;
            Flip(comb[y]);
            seq[i].resize(4,-1);
            seq[i][0]=comb[i];
        }

        int pos=0;
        SimplifySubSeq(seq, pos, 0, comb.size(), -1, -1);
        MovePos(seq, pos, comb.size());
        comb.clear();

        for(int s=0; s<seq.size(); s++)
        {
            if(seq[s][0]!=-1)
            {
                comb.push_back(seq[s][0]);
            }
        }
    }
    gain=t-comb.size();
    int pos=100;
    if(gain!=t && t!=0)
    {
        pos=1;
        while(100*gain>pos*t)
        {
            pos++;
        }
    }
    moy=pos;
	std::cout << "	* Number of edge flips after reduction : " << comb.size() << "  (gain : " << gain << "  i.e  " << moy << "%)" << std::endl;
}

//___________________________________________________________________________________________________________________________


void Mesh::SimplifySubSeq(vector<vector<int> >& seq, int& pos, const int& start, const int& end, const int& col1, const int& col2)
{
    if(start<end)
    {
        for(int i=start; i<end; i++)
        {
            if(seq[i][0]==col1 || seq[i][1]==col1 || seq[i][2]==col1 || seq[i][0]==col2 || seq[i][1]==col2 || seq[i][2]==col2)
            {
                if(FillOneColumn(seq, pos, i))
                    i--;
            }
        }
    }
}


//___________________________________________________________________________________________________________________________

void Mesh::EmptyColumn(vector<vector<int> >& seq, const int& i)
{
    for(int o=0; o<4; o++)
        seq[i][o]=-1;
}


//___________________________________________________________________________________________________________________________

bool Mesh::FillOneColumn(vector<vector<int> >& seq, int& pos, const int& col) //gauche à droite!!!
{
    if(seq[col][0]!=-1)
    {
        int other=-1;
        int i=col+1;
        while(other==-1 && i<comb.size())
        {
            if(seq[i][0]==seq[col][0])
                other=i;
            i++;
        }
        if(other==-1)
        {
            for(i=1; i<3; i++)
                seq[col][i]=-1;
        }
        else
        {
            MovePos(seq, pos, col);
            bool commuter=true;

            i=col+1;
            while(commuter && i<other)
            {
                commuter=false;
                if(seq[i][0]==-1)
                {
                    pos++;
                    commuter=true;
                    i++;
                }
                else
                {
                    if(edge[seq[col][0]].f[0]!=edge[seq[i][0]].f[0]
								&& edge[seq[col][0]].f[1]!=edge[seq[i][0]].f[0]
                                && edge[seq[col][0]].f[0]!=edge[seq[i][0]].f[1]
                                && edge[seq[col][0]].f[1]!=edge[seq[i][0]].f[1])
                    {
                        Flip(seq[i][0]);
                        pos++;
                        commuter=true;
                        i++;
                    }
                }
            }
            Flip(seq[col][0]);
            pos++;
            if(i==other)
            {
                seq[col][1]=-1;
                seq[col][2]=-1;

                Flip(seq[i][0]);
                pos++;
                int temp1=seq[col][3];
                if(seq[other][3]!=-1 && seq[other][3]<temp1)
                    temp1=seq[other][3];
                int temp2=seq[col][0];
                EmptyColumn(seq, other);
                EmptyColumn(seq, col);
                if(seq[col][3]!=-1)
                    SimplifySubSeq(seq, pos, temp1, other, temp2, temp2);
                return false;
            }
            else
            {
                seq[col][1]=seq[i][0];
                if(seq[i][3]==-1 || seq[i][3]>col)
                {
                    seq[i][3]=col;
                }
                if(seq[other][3]==-1 || seq[other][3]>col)
                {
                    seq[other][3]=col;
                }
            }

            MovePos(seq, pos, other-1);
            commuter=true;

            int j=other-1;
            while(commuter && i>col)
            {
                commuter=false;
                if(seq[j][0]==-1 || (edge[seq[other][0]].f[0]!=edge[seq[j][0]].f[0]
                                        && edge[seq[other][0]].f[1]!=edge[seq[j][0]].f[0]
										&& edge[seq[other][0]].f[0]!=edge[seq[j][0]].f[1]
                                        && edge[seq[other][0]].f[1]!=edge[seq[j][0]].f[1]))
                {
                    j--;
                    if(seq[j][0]!=-1)
                        Flip(seq[j][0]);
                    pos--;
                    commuter=true;
                }
            }
            if(j==col)
            {
                seq[col][1]=-1;
                seq[col][2]=-1;
            }
            else
            {
                seq[col][2]=seq[j][0];
                if(seq[j][3]==-1 || seq[j][3]>col)
                    seq[j][3]=col;
                if(j==i)
                {
                    return simplificationWithNoCommutativity(seq, pos, col, j, other);
                }
            }
        }
    }
    return false;
}

//___________________________________________________________________________________________________________________________


bool Mesh::simplificationWithNoCommutativity(vector<vector<int> >& seq, int& pos, const int& col, int& i, int& other)
{
    MovePos(seq, pos, col);

    int j=col;
    int temp;
    while(j<i-1)
    {
        for(int y=0; y<4; y++)
        {
            temp=seq[j][y];
            seq[j][y]=seq[j+1][y];
            seq[j+1][y]=temp;
        }
        if(seq[j][0]!=-1)
            Flip(seq[j][0]);
        pos++;
        j++;
    }

    j=other-1;
    while(j>i)
    {
        for(int y=0; y<4; y++)
        {
            temp=seq[j][y];
            seq[j][y]=seq[j+1][y];
            seq[j+1][y]=temp;
        }
        j--;
    }
    if(simplificationNoCommutativity(seq[i-1][0], seq[i][0]))
    {
        int min=seq[i-1][3];
        if(min==-1 || min>col)
            min=col;
        EmptyColumn(seq, i-1);
        int a=seq[i+1][0];
        int b=seq[i][0];
        for(int k=i+2; k<comb.size(); k++)
        {
            if(seq[k][0]==a)
            {
                seq[k][0]=b;
            }
            else
            {
                if(seq[k][0]==b)
                {
                    seq[k][0]=a;
                }
            }
        }
        if(min>seq[i+1][3] && seq[i+1][3]!=-1)
            min=seq[i+1][3];
        if(min>seq[i][3] && seq[i][3]!=-1)
            min=seq[i][3];

        if(min!=-1)
            SimplifySubSeq(seq, pos, min, i+1, a, b);
        return true;
    }
}

//_________________________________________________________________________________


void Mesh::MovePos(const vector<vector<int> >& seq, int& pos, const int& souhaite) //addapte la séquence de flips
{
    if(souhaite!=pos)
    {
        if(souhaite<pos)
        {
            for(int i=pos-1; i>souhaite-1; i--)
            {
                if(seq[i][0]!=-1)
                {
                    if(edge[seq[i][0]].v[0]==edge[seq[i][0]].v[1])
                        std::cout << "Forbiden Edge Flip -> " << i;
                    Flip(seq[i][0]);
                }
                pos--;
            }
        }
        else
        {
            for(int i=pos; i<souhaite; i++)
            {
                if(seq[i][0]!=-1)
                {
                    if(edge[seq[i][0]].v[0]==edge[seq[i][0]].v[1])
                        std::cout << "Forbiden Edge Flip -> " << i;
                    Flip(seq[i][0]);
                }
                pos++;
            }
        }
    }
}

//________________________________________________________________________


/********************************************************************/
/* Generer une sequence de bits à partir d'une séquence de bascules */
/********************************************************************/


void Mesh::CodageSimultaneousFlip2()
//Codage compact d'une séquence de bascules en une séquence de bits
{
	interval=100;
	
    vector<bool> cod;
    int stage=0;

    vector<int> now;
    now.resize(edge.size(), 0);
    vector<int> next;
    next.resize(edge.size(), -1);

    vector<int> last;
    last.resize(edge.size(), 0);

    if(comb.empty())
        cod.push_back(false);
    else
    {
        for(int i=comb.size()-1; i>-1; i--)
            Flip(comb[i]);

        vector<bool> Flips;
        Flips.resize(comb.size(),true);

        int flipped=0;

        vector<bool> Fait;

        while(flipped<comb.size())
        {
            stage++;
            cod.push_back(true);

            Fait.clear();
            Fait.resize(edge.size(), false);

            for(int it=0; it<comb.size(); it++)
            {
                if(Flips[it])
                {
                    if(!Fait[face[edge[comb[it]].f[0]].e[0]] && !Fait[face[edge[comb[it]].f[0]].e[1]] && !Fait[face[edge[comb[it]].f[0]].e[2]] && !Fait[face[edge[comb[it]].f[1]].e[0]] && !Fait[face[edge[comb[it]].f[1]].e[1]] && !Fait[face[edge[comb[it]].f[1]].e[2]])
                    {
                        Flip(comb[it]);
                        now[comb[it]]=1;
                        next[comb[it]]=-1;
                        for(int ff=0; ff<2; ff++)
                        {
                            for(int ee=0; ee<3; ee++)
                            {
                                if(face[edge[comb[it]].f[ff]].e[ee]!=comb[it])
                                {
                                    if(face[edge[comb[it]].f[ff]].e[ee]>comb[it])
                                    {
                                        now[face[edge[comb[it]].f[ff]].e[ee]]=-1;
                                    }
                                    next[face[edge[comb[it]].f[ff]].e[ee]]=0;
                                }
                            }
                        }

                        Flips[it]=false;
                        flipped++;
                    }
                    Fait[comb[it]]=true;
                }
            }

            for(int h=0; h<now.size(); h++)
            {
                if(now[h]==1)
				{
					//std::cout << " * Flip -> " << h << std::endl;
                    cod.push_back(true);
				}
                else
                {
                    if(now[h]==0)
                    {
						cod.push_back(false);
                    }
                }
            }
            last=now;
            now=next;
            next.clear();
            next.resize(edge.size(), -1);
        }
        cod.push_back(false);
    }
    double nF=0, nT=0;
	
	
	vtkArithmeticCoderBase *Coder = vtkArithmeticCoderBase::New();
	Coder->OpenFile("Sequence.flip", 1);
	Coder->StartCoding();
	
	qsmodel model;
	model.initqsmodel(2, 15, interval, NULL, 1);
		for(int j=0; j<cod.size(); j++)
		{
			
			if(cod[j])
			{
				Coder->Encode(1, &model);
				nT=nT+1;
			}
			if(!cod[j])
			{
				Coder->Encode(0, &model);
				nF=nF+1;
			}
		}
	
	
	int size = Coder->StopCoding();
	cout << "	" << size << " bytes written -> " << (double)8*size/vertex.size() << " bpv" << std::endl;
	Coder->CloseFile();
}

//________________________________________________________________________

void Mesh::ReadCodageSimultaneousFlip2()
//Lecture de la séquence de bits compact
{
	interval=100;
	
	vector<int> now;
	now.resize(edge.size(), 0);

	vector<int> next;
	next.resize(edge.size(), -1);

	vector<int> last;
	last.resize(edge.size(), -1);
		
	vtkArithmeticCoderBase *Coder = vtkArithmeticCoderBase::New();
	Coder->OpenFile("Sequence.flip", 0);
	Coder->StartDecoding();

	qsmodel model;
	model.initqsmodel(2, 15, interval, NULL, 0);
		
	int x;
	x=Coder->Decode(&model);
	//std::cout << "First -> " << x << std::endl;
	int stage=1;
	while(x==1)
	{
		for(int i=0; i<edge.size(); i++)
		{
			int f0=edge[i].f[0];
			int f1=edge[i].f[1];
			if(now[i]!=-1)
			{
				x=Coder->Decode(&model);
				//std::cout << x << std::endl;
				if(x==1)
				{
					//std::cout << "Flip -> " << i << std::endl;
					Flip(i);
					next[i]=-1;
					for(int ff=0; ff<2; ff++)
					{
						for(int ee=0; ee<3; ee++)
						{
							if(face[edge[i].f[ff]].e[ee]!=i)
							{
								if(face[edge[i].f[ff]].e[ee]>i)
								{
									now[face[edge[i].f[ff]].e[ee]]=-1;
								}
								next[face[edge[i].f[ff]].e[ee]]=0;
							}
						}
					}//
				}
			}
		}
		now=next;
		next.clear();
		next.resize(edge.size(), -1);
		
		x=Coder->Decode(&model);
		stage++;
	}
	Coder->StopDecoding();
}
