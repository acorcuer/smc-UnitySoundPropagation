using UnityEngine;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
public class meshTransport : MonoBehaviour {

	public string includeTag = "geo";
	public int numRays = 1000;
	public int maxPathLength = 100;
	public int maxNumReflecs = 75;
	public float absCoeff = 0.5f;
	public int numTriPerLeaf = 10;
	public bool debugEnable = false;
	public bool rescanFlag = false;
	public bool showWireFrame = false;
	public bool showSurfaceNormals = false;
	public bool showBoundingHeirarchy = false;
	public bool getRays = false;
	public bool showRays = false;
	private triangle[] debugTri;
	private Node[] debugNode;
	private Vector3[] rayOrigins, rayDirections;
	private float[] rayLengths;

	[DllImport("AudioPluginSpatializerTemplate")]
	private static extern void marshalGeomeTree (int numNodes,int numTri, int depth,int bbl, float[] boundingBoxes,int tl,float[] triangles, int lsl,int[] leafSizes,int tidl,int[] triangleIds,int tml,float[] triangleMatList);

	[DllImport("AudioPluginSpatializerTemplate")]
	private static extern void debugToggle (bool state);

	[DllImport("AudioPluginSpatializerTemplate")]
	private static extern void setTraceParam (int numberOfRays,int maxLen,int maxReflec,float absorbtionCoeff);

	[DllImport("AudioPluginSpatializerTemplate")]
	private static extern void getRayData (out int length, out IntPtr array);

	private float[] marshalRayData() {
		int arraySize;
		IntPtr arrayPtr;
		getRayData (out arraySize,out arrayPtr);
		float[] theArray =  new float[arraySize];
		Marshal.Copy(arrayPtr, theArray, 0, arraySize);
		Marshal.FreeCoTaskMem (arrayPtr);
		return theArray;
	}

	void Start () {
		debugToggle (debugEnable);
		setTraceParam (Mathf.FloorToInt (Mathf.Sqrt (numRays)),maxPathLength, maxNumReflecs,absCoeff);
		GeomeTree KDTree = calcTree ();
		sendTree (KDTree);
	}

	void Update() {
		if (rescanFlag) {
			GeomeTree KDTree = calcTree ();
			sendTree (KDTree);
			rescanFlag = false;
		}
		if (getRays) {
			float[] temp = marshalRayData ();
			rayOrigins = new Vector3[temp.Length / 7];
			rayDirections = new Vector3[temp.Length / 7];
			rayLengths = new float[temp.Length / 7];
			for (int i = 0; i < rayOrigins.Length; i++) {
				rayOrigins[i] = new Vector3(temp[i*7],temp[(i*7)+1],temp[(i*7)+2]);
				rayDirections[i] = new Vector3(temp[(i*7)+3],temp[(i*7)+4],temp[(i*7)+5]);
				rayLengths [i] = temp[(i*7)+6];
			}
			getRays = false;
		}
	}

	GeomeTree calcTree() {
		GameObject[] includedObjects = GameObject.FindGameObjectsWithTag (includeTag);
		Mesh combMesh = combineMeshes(includedObjects);
		triangle[] triArray = extractTris (combMesh);
		debugTri = triArray;
		GeomeTree theTree = new GeomeTree (3,numTriPerLeaf, triArray);
		return theTree;
	}

	Mesh combineMeshes(GameObject[] includedObjects) {
		CombineInstance[] combine = new CombineInstance[includedObjects.Length];

		for (int i = 0; i < includedObjects.Length; i++) {
			float coeff;
			Mesh tempMesh;
			if (includedObjects [i].name.Equals ("Right") || includedObjects [i].name.Equals ("Left")) {
				
			 tempMesh = Mesh.Instantiate (includedObjects [i].GetComponent<MeshFilter> ().sharedMesh)as Mesh;
			} else {
				tempMesh = includedObjects [i].GetComponent<MeshFilter> ().sharedMesh;
			}

			coeff = includedObjects [i].GetComponent<raytraceMaterial> ().absorbitonCoeff;
			Color[] colors = new Color[tempMesh.vertices.Length];
				
			for (int j = 0; j < colors.Length; j++) {
				if(includedObjects[i].name.Equals("Right")) {
					colors [j] = new Color (0,1,0,1);
				}else {
					if(includedObjects[i].name.Equals("Left")) {
						colors [j] = new Color (0, 0, 1,1);
					} else {

					colors [j] = new Color (coeff, 0, 0,1);
					}
				}
			}
			tempMesh.colors = colors;
			combine [i].mesh = tempMesh;
			combine [i].transform = includedObjects [i].GetComponent<MeshFilter> ().transform.localToWorldMatrix;
			}
		transform.GetComponent<MeshFilter>().mesh = new Mesh();
		transform.GetComponent<MeshFilter>().mesh.CombineMeshes(combine);
		transform.GetComponent<MeshRenderer> ().enabled = false;
		return transform.GetComponent<MeshFilter> ().mesh;

	}

	triangle[] extractTris(Mesh inputMesh) {
		int numTri = inputMesh.triangles.Length/3;
		int [] triIdx = inputMesh.triangles;
		Vector3[] vertVecs = inputMesh.vertices;
		triangle[] outputArray = new triangle[numTri];
		for (int i = 0; i < numTri; i++) {
			// Vectors in the mesh are stored on local scale (relating them to each other) in order to 
			// get meaningful coordinates transform.TransformPoint is used to convert them to world scale.
			// It basically looks at the transform of all parent objects and adds accordingly.
			Vector3 P1 = transform.TransformPoint(vertVecs [triIdx[i*3]]);
			Vector3 P2 = transform.TransformPoint(vertVecs [triIdx[(i*3)+1]]);
			Vector3 P3 = transform.TransformPoint(vertVecs [triIdx[(i*3)+2]]);
			Vector3 faceNorm = Vector3.Cross((P2 - P1),(P3 - P1));
			faceNorm.Normalize ();
			int listenerTag;

			if (inputMesh.colors [triIdx [i * 3]].g == 1) {
				listenerTag = 1;
			} else {
				if (inputMesh.colors [triIdx [i * 3]].b == 1) {
					listenerTag = 2;

				} else {
					listenerTag = 0;
				}
			}

			outputArray[i] = new triangle(P1,P2,P3,faceNorm,listenerTag,inputMesh.colors [triIdx[i*3]].r);
		}
		return outputArray;
	}

	void sendTree(GeomeTree inputTree) {
		Node[] nodeList = inputTree.getNodeList ().ToArray();
		debugNode = nodeList;
		int numTris = inputTree.numTriangles ();
		int numLeaves = inputTree.getLeafCount ();
		int numNodes = nodeList.Length;
		int depth = (int)(Mathf.Log (numNodes) / Mathf.Log (2.0f));  
		// init arrays to send to CPP, the binary tree has to be flattend into arrays
		// for easier marshalling
		float[] boundingBoxList = new float[numNodes*6];
		int[] leafSizeList = new int[numLeaves];
		int leafSizeIdx = 0;
		float[] triangleList = new float[numTris * 12];
		int[] triangleIdList = new int[numTris];
		float[] triangleMatList = new float[numTris];
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

					triangleIdList [triListIdx] = nodeTriangles [j].objectType;
					print (triangleIdList [triListIdx].ToString());
					triangleMatList [triListIdx] = nodeTriangles [j].absorbtionCoeff;
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
 		marshalGeomeTree (numNodes,numTris, depth, boundingBoxList.Length,boundingBoxList,triangleList.Length, triangleList,leafSizeList.Length,leafSizeList,triangleIdList.Length,triangleIdList,triangleMatList.Length,triangleMatList); 
	}

	void OnDrawGizmos() {
		Gizmos.color = Color.green;
		if (showWireFrame) {
			Gizmos.DrawWireMesh (GetComponent<MeshFilter> ().mesh);
		}
		Gizmos.color = Color.red;
		if (showSurfaceNormals) {
			for (int i = 0; i < debugTri.Length; i++) {
				Ray tempRay = new Ray((debugTri [i].P1 +debugTri [i].P2+debugTri [i].P3)/3, debugTri[i].facenorm);
				Gizmos.DrawRay(tempRay);
			}
		}	
		Gizmos.color = Color.blue;
		if (showBoundingHeirarchy) {
			for(int i = 0; i < debugNode.Length; i++) {
				Gizmos.DrawWireCube (debugNode [i].getBB ().center,debugNode [i].getBB ().size);
			}
		}
		Gizmos.color = Color.yellow;
		if (showRays) {
			for (int i = 0; i < rayDirections.Length; i++) {
				Gizmos.DrawLine (rayOrigins [i], rayOrigins [i] + (rayDirections [i] * rayLengths[i]));
			}
		}
	}
	 
}
	