#include "Maillage.hpp"

#include <string>
#include <fstream>
#include <ostream>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <math.h>
#include "Triangle.hpp"

using namespace std;

Point Maillage::operator[](int i) const
{
    if (i < 0 || i >= vPoint_.size()) {
        cout << "Indice incorrect pour le nombre de points" << endl;
        exit(1);
    }
    return(vPoint_[i]);
}

Triangle Maillage::operator()(int i) const
{
    if (i < 0 || i >= vTriangle_.size()) {
        cout << "Indice incorrect pour le nombre de triangles" << endl;
        exit(1);
    }
    return(vTriangle_[i]);
}

void Maillage::loadGMSH(const std::string & filename)
{
	vPoint_.resize(0);
	vTriangle_.resize(0);
	vTriangleAroundPoint_.resize(0);
	vPointAroundTriangle_.resize(0);
	char p_space = 32;


	ifstream infile;
	infile.open(filename.c_str(), ifstream::in);
	/*
	if !(infile) {
        cout << "Le fichier n'a pas ete trouve.";
        infile.close();
        exit(1);
	}*/
    string line;
    //while(!infile.eof()) {}
    getline(infile, line);
    while (line != "$Nodes")
    {
        getline(infile, line);
    }
    getline(infile, line);

    vector<string> line_split = split(line, p_space);
    int numEntityBlocks = atoi(line_split[0].c_str());
    int numNodes = atoi(line_split[1].c_str());
    int minNodeTag = atoi(line_split[2].c_str());
    int maxNodeTag = atoi(line_split[3].c_str());
    for (int i=0; i<numEntityBlocks; ++i) {
        getline(infile, line);
        vector<string> node_info = split(line, p_space);
        int entityDim = atoi(node_info[0].c_str());
        int entityTag = atoi(node_info[1].c_str());
        int parametric = atoi(node_info[2].c_str());
        int numNodesInBlock = atoi(node_info[3].c_str());

        if (entityDim != 2 && false) {
            for (int j=0; j<2*numNodesInBlock; ++j) {
                getline(infile, line);
            }
        }
        else{
            int prev_size = vPoint_.size();
            vPoint_.resize(prev_size + numNodesInBlock);
            for (int j=0; j<numNodesInBlock; ++j) {
                getline(infile, line);
            }
            for (int j=0; j<numNodesInBlock; ++j) {
                getline(infile, line);
                vector<string> node_coor = split(line, p_space);
                Point P(atof(node_coor[0].c_str()), atof(node_coor[1].c_str()), atof(node_coor[2].c_str()));
                vPoint_[prev_size + j] = P;
            }
        }
    }

	vTriangleAroundPoint_.resize(gNbPt());

    getline(infile, line);
    if (line != "$EndNodes") {
        cout << "Issue when reading .gmsh file" << endl;
    }
    getline(infile, line);
    getline(infile, line);
    line_split = split(line, p_space);
    numEntityBlocks = atoi(line_split[0].c_str());
    int numElements = atoi(line_split[1].c_str());
    int minElementTag = atoi(line_split[2].c_str());
    int maxElementTag = atoi(line_split[3].c_str());
    for (int i=0; i<numEntityBlocks; ++i) {
        getline(infile, line);
        vector<string> element_info = split(line, p_space);
        int entityDim = atoi(element_info[0].c_str());
        int entityTag = atoi(element_info[1].c_str());
        int elementType = atoi(element_info[2].c_str());
        int numElementsInBlock = atoi(element_info[3].c_str());

        if (elementType != 2) {
            for (int j=0; j<numElementsInBlock; ++j) {
                getline(infile, line);
            }
        }
        else{
            int prev_size = vTriangle_.size();
            vTriangle_.resize(prev_size + numElementsInBlock);
            vPointAroundTriangle_.resize(prev_size + numElementsInBlock);
            for (int j=0; j<numElementsInBlock; ++j) {
                getline(infile, line);
                vector<string> element_coor = split(line, p_space);
                Triangle T(vPoint_[atoi(element_coor[1].c_str())-1], vPoint_[atoi(element_coor[2].c_str())-1], vPoint_[atoi(element_coor[3].c_str())-1]);
                vTriangle_[prev_size + j] = T;

                vPointAroundTriangle_[j].resize(3);
                vPointAroundTriangle_[j][0] = atoi(element_coor[1].c_str())-1;
                vPointAroundTriangle_[j][1] = atoi(element_coor[2].c_str())-1;
                vPointAroundTriangle_[j][2] = atoi(element_coor[3].c_str())-1;

                vTriangleAroundPoint_[atoi(element_coor[1].c_str())-1].push_back(j);
                vTriangleAroundPoint_[atoi(element_coor[2].c_str())-1].push_back(j);
                vTriangleAroundPoint_[atoi(element_coor[3].c_str())-1].push_back(j);
            }

        }
    }

	infile.close();
}

bool Maillage::crossTriangle(const Point & P)
{
    bool crossT(false);

    for (int i=0; i<gNbTri(); ++i) {
        if (vTriangle_[i].isIntersected(P)) {
            crossT = true;
        }
    }
    return(crossT);
}

bool Maillage::crossTriangle(const Point & P, const Point & Q)
{
    bool crossT(false);

    for (int i=0; i<gNbTri(); ++i) {
        if (vTriangle_[i].isIntersected(P, Q)) {
            crossT = true;
        }
    }
    return(crossT);
}

bool Maillage::crossTriangleObs(const Point & P_obs, const Point & P_pixel)
{
    bool crossT(false);

    for (int i=0; i<gNbTri(); ++i) {
        if (vTriangle_[i].isIntersectedObs(P_obs, P_pixel)) {
            crossT = true;
        }
    }
    return(crossT);
}

bool Maillage::crossOtherTriangle(int idTriangle, const Point & P, const Point & Q)
{
    bool crossT(false);

    for (int i=0; i<gNbTri(); ++i) {
        if (vTriangle_[i].isIntersected(P, Q) && i != idTriangle) {
            crossT = true;
        }
    }
    return(crossT);
}

Triangle Maillage::closestIntersectedTriangle(const Point & P)
{
    Point P_intersection;
    double norm_intersection;
    Triangle T_min;
    double norm_min(numeric_limits<double>::infinity());

    for (int i=0; i<gNbTri(); ++i) {
        if (vTriangle_[i].isIntersected(P)) {
            P_intersection = vTriangle_[i].ptIntersection(P);
            norm_intersection = norm_euc(P_intersection-P);
            if (norm_intersection < norm_min) {
                norm_min = norm_intersection;
                T_min = vTriangle_[i];
            }
        }
    }
    return(T_min);
}

Triangle Maillage::closestIntersectedTriangle(const Point & P, const Point & Q)
{
    Point P_intersection;
    double norm_intersection;
    Triangle T_min;
    double norm_min(numeric_limits<double>::infinity());

    for (int i=0; i<gNbTri(); ++i) {
        if (vTriangle_[i].isIntersected(P, Q)) {
            P_intersection = vTriangle_[i].ptIntersection(P, Q);
            norm_intersection = norm_euc(P_intersection-P);
            if (norm_intersection < norm_min) {
                norm_min = norm_intersection;
                T_min = vTriangle_[i];
            }
        }
    }
    return(T_min);
}

int Maillage::closestIntersectedTriangleId(const Point & P)
{
    Point P_intersection;
    double norm_intersection;
    int id_T_min(-1);
    double norm_min(numeric_limits<double>::infinity());

    for (int i=0; i<gNbTri(); ++i) {
        if (vTriangle_[i].isIntersected(P)) {
            P_intersection = vTriangle_[i].ptIntersection(P);
            norm_intersection = norm_euc(P_intersection-P);
            if (norm_intersection < norm_min) {
                norm_min = norm_intersection;
                id_T_min = i;
            }
        }
    }
    return(id_T_min);
}

int Maillage::closestIntersectedTriangleId(const Point & P, const Point & Q)
{
    Point P_intersection;
    double norm_intersection;
    int id_T_min(-1);
    double norm_min(numeric_limits<double>::infinity());

    for (int i=0; i<gNbTri(); ++i) {
        if (vTriangle_[i].isIntersected(P, Q)) {
            P_intersection = vTriangle_[i].ptIntersection(P, Q);
            norm_intersection = norm_euc(P_intersection-P);
            if (norm_intersection < norm_min) {
                norm_min = norm_intersection;
                id_T_min = i;
            }
        }
    }
    return(id_T_min);
}

int Maillage::closestIntersectedTriangleIdObs(const Point & P_obs, const Point & P_pixel)
{
    Point P_intersection;
    double norm_intersection;
    int id_T_min(-1);
    double norm_min(numeric_limits<double>::infinity());

    for (int i=0; i<gNbTri(); ++i) {
        if (vTriangle_[i].isIntersectedObs(P_obs, P_pixel)) {
            P_intersection = vTriangle_[i].ptIntersectionObs(P_obs, P_pixel);
            norm_intersection = norm_euc(P_intersection-P_obs);
            if (norm_intersection < norm_min) {
                norm_min = norm_intersection;
                id_T_min = i;
            }
        }
    }
    return(id_T_min);
}

std::vector<int> Maillage::getPointAroundTriangle(int i)
{
    return(vPointAroundTriangle_[i]);
}

std::vector<int> Maillage::getTriangleAroundPoint(int i)
{
    return(vTriangleAroundPoint_[i]);
}

Point Maillage::vectNormPoint(int idxPoint)
{
    Point normPoint;
    for (int i=0; i<vTriangleAroundPoint_[idxPoint].size(); ++i) {
        normPoint += vTriangle_[vTriangleAroundPoint_[idxPoint][i]].vectNormal() / vTriangleAroundPoint_[idxPoint].size();
    }
    return(normPoint);
}

vector<double> Maillage::angleOmbrage(int id_P, const Point & P_lum, const Point & P_pixel)
{
    vector<double> vAng(2);
    Point vectNorm = vectNormPoint(id_P);
    vAng[0] = acos((vectNorm | (P_lum - vPoint_[id_P])) / norm_euc(P_lum - vPoint_[id_P])); //theta
    vAng[1] = acos((vectNorm | (P_pixel - vPoint_[id_P]))/norm_euc(P_pixel - vPoint_[id_P])) - vAng[0]; //alpha
    return(vAng);
}

/*
function [Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nomfile)
fid=fopen(nomfile,'r');
if fid <=0,
   msg=['Le fichier de maillage : ' nomfile ' n''a pas ?t? trouv?'];
   error(msg);
end

while ~strcmp(fgetl(fid),'$Nodes'), end
Nbpt = str2num(fgetl(fid));
Coorneu = zeros(Nbpt,2);
Refneu = zeros(Nbpt,1);
Nbaretes = 0;
Numaretes = [];
Refaretes = [];
RefneuBis = zeros(Nbpt,1);
for i=1:Nbpt
tmp= str2num(fgetl(fid));
Coorneu(i,:) = tmp(2:3);
end
while ~strcmp(fgetl(fid),'$Elements'), end
Nbtri = str2num(fgetl(fid));
tmp= str2num(fgetl(fid));test = tmp(2);
while ( (test~=15) && (test~=1) && (test~=2)), tmp= str2num(fgetl(fid)); test = tmp(2); Nbtri = Nbtri-1;end
% on traite les noeuds des coins
while ((test~=1) && (test~=2))
     Refneu(tmp(end))=tmp(end-2);
     tmp= str2num(fgetl(fid)); test = tmp(2); Nbtri = Nbtri-1;
 end
% on traite les noeuds du bord
k=1;
while (test~=2)
    RefneuBis(tmp(6:7)')=tmp(4);
    Numaretes= [Numaretes;tmp(6:7)];
    Refaretes= [Refaretes;tmp(4)];
    tmp= str2num(fgetl(fid));
    test = tmp(2); Nbaretes=Nbaretes+1; Nbtri=Nbtri-1;k=k+1;
end
Refneu(find(Refneu==0))=RefneuBis(find(Refneu==0));
Numtri = zeros(Nbtri,3);
Reftri = zeros(Nbtri,1);
% on traite les references des triangles
for i=1:Nbtri
    Numtri(i,:) = tmp(end-2:end);
    Reftri(i)=tmp(end-3);
    tmp= str2num(fgetl(fid));
end
fclose(fid);
*/
