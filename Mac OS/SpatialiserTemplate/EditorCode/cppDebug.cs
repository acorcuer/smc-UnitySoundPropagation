using UnityEngine;
using UnityEditor;
using System;
using System.Runtime.InteropServices;

[InitializeOnLoad]
public class cppDebug : MonoBehaviour {

	private delegate void DebugCallback(string message);

	[DllImport("AudioPluginSpatializerTemplate")]
	private static extern void RegisterDebugCallback(DebugCallback callback);

	static cppDebug()
	{
		RegisterDebugCallback(new DebugCallback(DebugMethod));
	}

	private static void DebugMethod(string message)
	{
		Debug.Log("FromCPP: " + message);
	}
}