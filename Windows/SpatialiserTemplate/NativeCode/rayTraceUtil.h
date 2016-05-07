#ifndef rayTraceUtil_h
#define rayTraceUtil_h
#include <vector>
#include <deque>
#include <string>
#include <sstream>
#include <utility>


#if UNITY_WIN
#define ABA_API __declspec(dllexport)
#else
#define ABA_API
#endif

typedef void(*DebugCallback) (const char *str);
DebugCallback gDebugCallback;

extern "C" ABA_API void RegisterDebugCallback (DebugCallback callback)
{
    if (callback)
    {
        gDebugCallback = callback;
    }
}

void DebugInUnity(std::string message)
{
    if (gDebugCallback)
    {
        gDebugCallback(message.c_str());
    }
}

const float kEpsilon = 1e-8;


class Vector3 {
public:
    float X, Y, Z;
    inline Vector3(){};
    inline Vector3(float n_x, float n_y, float n_z)
    { X = n_x; Y = n_y; Z = n_z;}
    inline void setPos(float n_x, float n_y, float n_z)
    {X = n_x; Y = n_y; Z = n_z;}
    inline Vector3 operator + (const Vector3& A) const
    {return Vector3(X + A.X, Y + A.Y, Z + A.Z); }
    
    inline Vector3 operator * (const float& A) const
    {return Vector3(X * A, Y * A, Z * A); }
    inline Vector3 operator - (const Vector3& A) const
    {return Vector3(X - A.X, Y - A.Y, Z - A.Z); }
    inline float Dot( const Vector3& A ) const
    { return A.X*X + A.Y*Y + A.Z*Z; }
    inline Vector3 cross(const Vector3 &A) const
    {return Vector3(Y * A.Z - Z * A.Y,Z * A.X - X * A.Z,X * A.Y - Y * A.X);}
    inline bool operator<(const Vector3 &A) const {
        return (X < A[0] && Y < A[1] && Z < A[2]);
    }
    inline float operator [] (const int& A) const{
        switch(A){
        case 0:
        return X;
        break;
        case 1:
        return Y;
        break;
        case 2:
        return Z;
        break;
        default:
        return NULL;
        }
    }
};


class Tri {
public:
    Vector3 P1,P2,P3,faceNorm;
    int isListener;
    inline Tri(){};
    inline Tri(Vector3 n_p1,Vector3 n_p2,Vector3 n_p3,Vector3 n_fn,int n_ls)
    {P1 = n_p1; P2 = n_p2; P3 = n_p3; faceNorm = n_fn;isListener = n_ls;}
    inline void setVerts(Vector3 n_p1,Vector3 n_p2,Vector3 n_p3)
    {P1 = n_p1; P2 = n_p2; P3 = n_p3;}
    inline void setFaceNorm(Vector3 n_fn)
    {faceNorm = n_fn;}
};


class Ray {
public:
    Vector3 origin,direction,invDirection;
    float pathLength = 0;
    int sign[3];
    inline Ray(){};
    inline Ray(Vector3 n_o, Vector3 n_d)
    {origin = n_o;
        direction = n_d;
        invDirection = Vector3(1/n_d.X, 1/n_d.Y, 1/n_d.Z);
        sign[0] = (invDirection.X < 0);
        sign[1] = (invDirection.Y < 0);
        sign[2] = (invDirection.Z < 0);
    }
    inline void setOrigin(Vector3 n_o)
    {origin = n_o;}
    inline void setDirection(Vector3 n_d)
    {direction = n_d;
        invDirection = Vector3(1/n_d.X, 1/n_d.Y, 1/n_d.Z);
    }
    inline void addLength(float add)
    {pathLength += add;}
	/*inline float getAngle(Vector3 normal)
	{
		float angle = acosf(direction.Dot(normal)/sqrtf(direction.lenSq * normal.lenSq));
		return angle;
	}*/
	inline Vector3 getNewDirection(Vector3 normal)
	{
		Vector3 newDir = direction - normal * 2.0f * direction.Dot(normal);
		return newDir;
	}
    inline bool testIntersect(Tri *tri,float &t,float &u, float&v){
       
        Vector3 v0v1 = tri->P2 - tri->P1;
        Vector3 v0v2 = tri->P3 - tri->P1;
        Vector3 pvec = direction.cross(v0v2);
        float det = v0v1.Dot(pvec);
        
        if(fabs(det) < kEpsilon) return false;
        float invDet = 1/det;
        
        Vector3 tvec = origin - tri->P1;
        u = tvec.Dot(pvec) * invDet;
        if(u < 0 || u > 1) return false;
    
        Vector3 qvec = tvec.cross(v0v1);
        v = direction.Dot(qvec) * invDet;
        if(v < 0 || u + v > 1) return false;
        t = v0v2.Dot(qvec)*invDet;
        return true;
    }
    
};

class Bounds {
public:
    Vector3 parameters[2];
    inline Bounds() {}
    inline Bounds(Vector3 n_max,Vector3 n_min){
        assert(n_min < n_max);
        parameters[0] = n_min;
        parameters[1] = n_max;
    }
    inline void setBB(Vector3 n_max,Vector3 n_min){
        assert(n_min < n_max);
        parameters[0] = n_min;
        parameters[1] = n_max;
    }
    inline bool testIntersect(Ray *r, float t0, float t1) {
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
       
        tmin = (parameters[r->sign[0]].X - r->origin.X) * r->invDirection.X;
        tmax = (parameters[1-r->sign[0]].X - r->origin.X) * r->invDirection.X;
        tymin = (parameters[r->sign[1]].Y - r->origin.Y) * r->invDirection.Y;
        tymax = (parameters[1-r->sign[1]].Y - r->origin.Y) * r->invDirection.Y;
       
        if ( (tmin > tymax) || (tymin > tmax) )
            return false;
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;
        tzmin = (parameters[r->sign[2]].Z - r->origin.Z) * r->invDirection.Z;
        tzmax = (parameters[1-r->sign[2]].Z - r->origin.Z) * r->invDirection.Z;
        if ( (tmin > tzmax) || (tzmin > tmax) )
            return false;
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;
        return ( (tmin < t1) && (tmax > t0) );
    }
};

class Node {
private:
    Bounds boundingBox;
    int depth;
    std::vector<Tri> nodeTriangles;
    bool isLeaf = false;
public:
    Node *leftChild,*rightChild;
    inline Node() {}
    inline Node(int parentDepth, int maxDepth, std::deque<float> *boundingBoxArray,std::deque<float> *triangleArray,std::deque<int> *leafSizeArray,std::deque<int> *idList) {
        this->depth = parentDepth+1;
        this->boundingBox = Bounds(Vector3(boundingBoxArray->at(0),boundingBoxArray->at(1),boundingBoxArray->at(2)),Vector3(boundingBoxArray->at(3),boundingBoxArray->at(4),boundingBoxArray->at(5)));
        boundingBoxArray->erase(boundingBoxArray->begin(), boundingBoxArray->begin()+6);
        if(this->depth > maxDepth) {
            this->isLeaf = true;
            int sizeOfLeaf = leafSizeArray->front();
            leafSizeArray->pop_front();
            nodeTriangles.resize(sizeOfLeaf);
            for (int i=0; i < sizeOfLeaf; i++) {
                Vector3 P1 = Vector3(triangleArray->at(0),triangleArray->at(1),triangleArray->at(2));
                Vector3 P2 = Vector3(triangleArray->at(3),triangleArray->at(4),triangleArray->at(5));
                Vector3 P3 = Vector3(triangleArray->at(6),triangleArray->at(7),triangleArray->at(8));
                Vector3 faceNorm = Vector3(triangleArray->at(9),triangleArray->at(10),triangleArray->at(11));
                triangleArray->erase(triangleArray->begin(), triangleArray->begin()+12);
                nodeTriangles[i] = *new Tri(P1,P2,P3,faceNorm,idList->front());
                idList->pop_front();
            }
        }else{
            leftChild = new Node(depth,maxDepth,boundingBoxArray,triangleArray,leafSizeArray,idList);
            rightChild = new Node(depth,maxDepth,boundingBoxArray,triangleArray,leafSizeArray,idList);
        }
    }
    inline bool intersects(Ray *testRay){
        return boundingBox.testIntersect(testRay, 0, INFINITY);
    }
    inline bool getIsLeaf() {
        return isLeaf;
    }
    inline std::vector<Tri>* getTri() {
        return &nodeTriangles;
    }
    inline int numTri() {
        return nodeTriangles.size();
    }
};

class GeomeTree {
    public:
    Node masterNode;
    inline GeomeTree() {};
    inline GeomeTree(int maxDepth,std::deque<float> *boundingBoxArray,std::deque<float> *triangleArray,std::deque<int> *leafSizeArray,std::deque<int> *idList) {
        masterNode = *new Node(0,maxDepth,boundingBoxArray,triangleArray,leafSizeArray,idList);
    };
    
    inline std::deque<Tri> getCandidates(Ray *testRay) {
        std::deque<Tri> candidates;
        collision(&masterNode, testRay, &candidates);
        return candidates;
    }
    
    inline void collision(Node *currentNode, Ray *testRay, std::deque<Tri> *candidates) {
        if(currentNode->intersects(testRay)) {
            if (currentNode->getIsLeaf()) {
                if(currentNode->numTri() != 0) {
                    for(int i = 0; i < currentNode->getTri()->size();i++){
                        candidates->push_back(currentNode->getTri()->at(i));
                    }
                }
            }else {
                collision(currentNode->leftChild, testRay, candidates);
                collision(currentNode->rightChild, testRay, candidates);
        
        }
    }
}
};

    class raySphere {
    public:
        std::vector<Ray> rays;
        inline raySphere(){};
        inline raySphere(int nRays) {
            float rnd1 = rand() % 2;
            float rnd2 = rand() % 2;
            for(int i = 0; i < nRays; i++) {
                for (int j = 0; j < nRays; j++) {
            float l = 2 * sqrt(((i + rnd1) / nRays) - pow(((i + rnd1) / nRays), 2.0))*cos(2 * M_PI*((j + rnd2) / nRays));
            float m = 2 * sqrt(((i + rnd1) / nRays) - pow(((i + rnd1) / nRays), 2.0))*sin(2 * M_PI*((j + rnd2) / nRays));
            float n = 1 - 2 * ((i + rnd1) / nRays);
            Vector3 directionVector = Vector3(l, m, n);
            Vector3 positionVector = directionVector;;
            rays.push_back(Ray(positionVector,directionVector));
                }
            }
        }
        inline std::vector<Ray> getRayList(Vector3 sourcePos){
            std::vector<Ray> translatedRays;
            for(int i = 0; i< rays.size(); i++){
                Ray tempRay = Ray(rays[i].origin + sourcePos, rays[i].direction);
                translatedRays.push_back(tempRay);
            }
            return translatedRays;
        }
};



#endif /* rayTraceUtil_h */
