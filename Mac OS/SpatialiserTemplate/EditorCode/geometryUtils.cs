using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;


// Disable warning, is a byproduct of using classes recursively as a datastructure (the leaves childeren never get "used")
// The private field Node.leftChild is assigned but never used (lines 72 and 73)
#pragma warning disable 0414 

public class triangle {
	// X,Y,Z of 3 triangle corners and the face normal 
	public Vector3 P1,P2,P3,facenorm;
	public int isListener;
	public triangle() {}
	public triangle(Vector3 n_p1,Vector3 n_p2,Vector3 n_p3, Vector3 n_fn,int n_ls){
		P1 = n_p1;
		P2 = n_p2;
		P3 = n_p3;
		facenorm = n_fn;
		isListener = n_ls;
	}
	public void setData(Vector3 n_p1,Vector3 n_p2,Vector3 n_p3,Vector3 n_fn,int n_ls) {
		P1 = n_p1;
		P2 = n_p2;
		P3 = n_p3;
		facenorm = n_fn;
		isListener = n_ls;
	}
}
	
public class GeomeTree {
	// dim - How many dimensions to partition the triangles in
	// bucketSize - Maximum number of triangles in a leaf
	private int dim,bucketSize,numTri;
	private int leafCount = 0;
	private List<Node> nodeList = new List<Node>();
	public Node masterNode;

	public GeomeTree(int inputDim,int inputBucket,triangle[] inputTriangles) {
		dim = inputDim;
		bucketSize = inputBucket;
		numTri = inputTriangles.Length;
		// Create the master node, additonal nodes are recursively created under this one.
		// Leaf count and nodeList are references as the recursive algorithm makes it
		// difficult to store these accurately and safely (use of static is possible, but opens 
		// up more issues than it is worth)
		masterNode = new Node (inputTriangles, dim, bucketSize,0,ref nodeList,ref leafCount,0);
	}

	public int getLeafCount() {
		return leafCount;
	}

	public List<Node> getNodeList() {
		return nodeList;
	}

	public int numTriangles() {
		return numTri;
	}
}

public class Node {
	private int depth = 0;
	private Bounds boundingBox;
	private Node leftChild, rightChild;
	private triangle[] nodeTriangles = null;
	private bool isLeaf = false;

	public Node(triangle[] inputTriangles, int dim, int bucketSize, int axis,ref List<Node> nodeList, ref int leafCount,int parentDepth){
		depth = parentDepth + 1;
		// Resize the bounding box of this node so that it includes the triangles 
		// BB's are initialied with no size.
		resizeBB (inputTriangles);
		// If the number of triangles passed to this node is less than the maximum
		// mark this node as a leaf and assign the triangles to it
		if (inputTriangles.Length <= bucketSize) {
			leafCount++;
			nodeTriangles = inputTriangles;
			isLeaf = true;
			nodeList.Add (this);
		// If the number of triangles is greater, split the triangles along
			// the next dimension and create left/right nodes(recursion).
		}  else {
			nodeList.Add (this);
			triangle[] leftTris, rightTris;
			splitBB (inputTriangles,axis++, out leftTris, out rightTris);				
			leftChild = new Node (leftTris, dim, bucketSize,axis % dim, ref nodeList, ref leafCount,depth);
			rightChild = new Node (rightTris, dim, bucketSize,axis % dim, ref nodeList, ref leafCount,depth);
		}
	}

	private void resizeBB(triangle[] inputTriangles) {
		float avgX = 0;
		float avgY = 0;
		float avgZ = 0;
		for (int i = 0; i < inputTriangles.Length; i++) {
			avgX += inputTriangles [i].P1.x;
			avgX += inputTriangles [i].P2.x;
			avgX += inputTriangles [i].P3.x;
			avgY += inputTriangles [i].P1.y;
			avgY += inputTriangles [i].P2.y;
			avgY += inputTriangles [i].P3.y;
			avgZ += inputTriangles [i].P1.z;
			avgZ += inputTriangles [i].P2.z;
			avgZ += inputTriangles [i].P3.z;
		}
		avgX /= (inputTriangles.Length * 3);
		avgY /= (inputTriangles.Length * 3);
		avgZ /= (inputTriangles.Length * 3);
		boundingBox.center = new Vector3 (avgX, avgY, avgZ);
		for (int i = 0; i < inputTriangles.Length; i++) {
			boundingBox.Encapsulate (inputTriangles [i].P1);
			boundingBox.Encapsulate (inputTriangles [i].P2);
			boundingBox.Encapsulate (inputTriangles [i].P3);
		}
	}

	private void splitBB(triangle[] inputTriangles,int axis, out triangle[] leftTris, out triangle[] rightTris) {
		triangle[] sortedTriangles = sortAlongDim (inputTriangles, axis);
		leftTris = new triangle[inputTriangles.Length / 2];
		rightTris = new triangle[inputTriangles.Length - leftTris.Length];
		Array.Copy(sortedTriangles, leftTris, leftTris.Length);
		Array.Copy(sortedTriangles, leftTris.Length, rightTris, 0, rightTris.Length);
	}

	private triangle[] sortAlongDim(triangle[] inputTriangles, int axis) {
		int len = inputTriangles.Length;
		float[] sortArray = new float[len];
		int[] sortOrder = new int[len];
		triangle[] sortedTriangles = new triangle[len];
		for(int i = 0; i < len; i++) {
			float max = inputTriangles [i].P1 [axis];
			if(inputTriangles [i].P2 [axis] > max) {
				max = inputTriangles [i].P2 [axis];
			}
			if(inputTriangles [i].P3 [axis] > max) {
				max = inputTriangles [i].P3 [axis];
			}
			sortOrder [i] = i;
			sortArray [i] =max;
		}
		Array.Sort (sortArray,sortOrder);
		for(int i = 0; i < len; i++) {
			sortedTriangles [i] = inputTriangles[sortOrder [i]];
		}
		return sortedTriangles;
	}
		
	public bool getIsLeaf() {
		return isLeaf;
	}

	public triangle[] getTriangles() {
		return nodeTriangles;
	}

	public Bounds getBB() {
		return boundingBox;
	}
}

