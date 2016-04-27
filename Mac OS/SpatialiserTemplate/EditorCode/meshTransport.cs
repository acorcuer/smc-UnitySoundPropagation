using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
public class meshTransport : MonoBehaviour {

	public string includeTag = "geo";
	public int numTriPerLeaf = 10;
	public bool rescanFlag = false;
	public bool showCombinedMesh = false;
		
	[DllImport("AudioPluginSpatializerTemplate", CallingConvention = CallingConvention.StdCall)]
	private static extern void marshalGeomeTree (int numNodes,int numTri, int depth,int bbl, float[] boundingBoxes,int tl,float[] triangles, int lsl,int[] leafSizes);
	  
	void Start () {
		GeomeTree KDTree = calcTree ();
		sendTree (KDTree);
	}

	void Update() {
		if (rescanFlag) {
			GeomeTree KDTree = calcTree ();
			sendTree (KDTree);
			rescanFlag = false;
		}
		transform.GetComponent<MeshRenderer> ().enabled =  showCombinedMesh;
	}

	GeomeTree calcTree() {
		GameObject[] includedObjects = GameObject.FindGameObjectsWithTag (includeTag);
		Mesh combMesh = combineMeshes(includedObjects);
		triangle[] triArray = extractTris (combMesh);
		GeomeTree theTree = new GeomeTree (3,numTriPerLeaf, triArray);
		return theTree;
	}

	Mesh combineMeshes(GameObject[] includedObjects) {
		CombineInstance[] combine = new CombineInstance[includedObjects.Length];
		for (int i = 0; i < includedObjects.Length; i++) {
			MeshFilter thisMeshfilt = includedObjects [i].GetComponent<MeshFilter> ();
			combine [i].mesh = thisMeshfilt.sharedMesh;
			combine[i].transform = thisMeshfilt.transform.localToWorldMatrix;
		}
		transform.GetComponent<MeshFilter>().mesh = new Mesh();
		transform.GetComponent<MeshFilter>().mesh.CombineMeshes(combine);
		transform.GetComponent<MeshRenderer> ().material = new Material (Shader.Find ("Diffuse"));
		return transform.GetComponent<MeshFilter> ().mesh;
	}

	triangle[] extractTris(Mesh inputMesh) {
		int numTri = inputMesh.triangles.Length/3;
		int [] triIdx = inputMesh.triangles;
		Vector3[] vertVecs = inputMesh.vertices;
		Vector3[] normVecs = inputMesh.normals;
		triangle[] outputArray = new triangle[numTri];
		for (int i = 0; i < numTri; i++) {
			// Vectors in the mesh are stored on local scale (relating them to each other) in order to 
			// get meaningful coordinates transform.TransformPoint is used to convert them to world scale.
			// It basically looks at the transform of all parent objects and adds accordingly.
			Vector3 P1 = transform.TransformPoint(vertVecs [triIdx[i*3]]);
			Vector3 P2 = transform.TransformPoint(vertVecs [triIdx[(i*3)+1]]);
			Vector3 P3 = transform.TransformPoint(vertVecs [triIdx[(i*3)+2]]);
			Vector3 faceNorm = transform.TransformPoint((normVecs [triIdx[i*3]] + normVecs [triIdx[(i*3)+1]] + normVecs [triIdx[(i*3)+2]])/3);
			outputArray[i] = new triangle(P1,P2,P3,faceNorm);
		}
		return outputArray;
	}

	void sendTree(GeomeTree inputTree) {
		Node[] nodeList = inputTree.getNodeList ().ToArray();
		int numTris = inputTree.numTriangles ();
		int numLeaves = inputTree.getLeafCount ();
		int numNodes = (numLeaves * 2) - 1;
		int depth = (int)Mathf.Sqrt ((numNodes + 1));
		// init arrays to send to CPP, the binary tree has to be flattend into arrays
		// for easier marshalling
		float[] boundingBoxList = new float[numNodes*6];
		int[] leafSizeList = new int[numLeaves];
		int leafSizeIdx = 0;
		float[] triangleList = new float[numTris * 12];
		int triListIdx = 0;
		for(int i = 0; i < nodeList.Length; i++) {
			Bounds tempBounds = nodeList [i].getBB ();
			Vector3 tempMax = tempBounds.max;
			Vector3 tempMin = tempBounds.min;
			boundingBoxList [i * 6] 	  = tempMax.x;
			boundingBoxList [(i * 6) + 1] = tempMax.y;
			boundingBoxList [(i * 6) + 2] = tempMax.z;
			boundingBoxList [(i * 6) + 3] = tempMin.x;
			boundingBoxList [(i * 6) + 4] = tempMin.y;
			boundingBoxList [(i * 6) + 5] = tempMin.z;
			if(nodeList[i].getIsLeaf()) {
				triangle[] nodeTriangles = nodeList [i].getTriangles();
				leafSizeList [leafSizeIdx++] = nodeTriangles.Length;
				for (int j = 0; j < nodeTriangles.Length; j++) {
					triangleList [triListIdx*12] = nodeTriangles [j].P1.x;
					triangleList [(triListIdx*12)+1] = nodeTriangles [j].P1.y;
					triangleList [(triListIdx*12)+2] = nodeTriangles [j].P1.z;
					triangleList [(triListIdx*12)+3] = nodeTriangles [j].P2.x;
					triangleList [(triListIdx*12)+4] = nodeTriangles [j].P2.y;
					triangleList [(triListIdx*12)+5] = nodeTriangles [j].P2.z;
					triangleList [(triListIdx*12)+6] = nodeTriangles [j].P3.x;
					triangleList [(triListIdx*12)+7] = nodeTriangles [j].P3.y;
					triangleList [(triListIdx*12)+8] = nodeTriangles [j].P3.z;
					triangleList [(triListIdx*12)+9] = nodeTriangles [j].facenorm.x;
					triangleList [(triListIdx*12)+10] = nodeTriangles [j].facenorm.y;
					triangleList [(triListIdx*12)+11] = nodeTriangles [j].facenorm.z;
					triListIdx++;
				}
			}
		}
		print (numTris.ToString ());
		marshalGeomeTree (numNodes,numTris, depth, boundingBoxList.Length,boundingBoxList,triangleList.Length, triangleList,leafSizeList.Length, leafSizeList); 
	}
}
	