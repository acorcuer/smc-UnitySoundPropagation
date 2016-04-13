using UnityEngine;
using System;
using System.Collections;
using System.Runtime.InteropServices;

public class meshTransport : MonoBehaviour {

	public bool rescanFlag = false;

	[DllImport("AudioPluginSpatializerTemplate", CallingConvention = CallingConvention.StdCall)]
	private static extern void setStruct (int numVertsIn,float[] vertArrayIn,float[] normArrayIn,float[] sizeArrayIn);

	[StructLayout(LayoutKind.Sequential,Pack = 1),Serializable]
	private struct meshData {
		// Fields
		private int numVerts;
		private float[] verts;
		private float[] norms;
		private float[] size;
		// Constructor
		public meshData(Mesh inputMesh) {	
			numVerts = inputMesh.vertexCount;
			Vector3[] vertVec = inputMesh.vertices;
			Vector3[] normVec = inputMesh.normals;
			verts = new float[3*numVerts];
			norms = new float[3*numVerts];
			size = new float[3];
			size[0] = inputMesh.bounds.size.x;
			size[1] = inputMesh.bounds.size.y;
			size[2] = inputMesh.bounds.size.z;
			int count = 0;
			for(int i = 0; i < numVerts; i++) {
				verts[count] = vertVec[i].x;
				verts[count+1] = vertVec[i].y;
				verts[count+2] = vertVec[i].z;
				norms[count] = normVec[i].x;
				norms[count+1] = normVec[i].y;
				norms[count+2] = normVec[i].z;
				count = count+3;
			}
		}
		// Methods
		public int meshLength() {
			return verts.GetLength (1);
		}
		// Call function in plugin that sets C++ struct to values c# struct
		public void transportToPlugin() {
			setStruct (numVerts, verts, norms, size);

		}
}
	void Start () {
		meshData geometry = new meshData(GetComponent<MeshFilter> ().mesh);
		geometry.transportToPlugin ();
	}

	void Update() {
		if (rescanFlag) {
			meshData geometry = new meshData(GetComponent<MeshFilter> ().mesh);
			geometry.transportToPlugin ();
			rescanFlag = false;
		}
	}
}
	