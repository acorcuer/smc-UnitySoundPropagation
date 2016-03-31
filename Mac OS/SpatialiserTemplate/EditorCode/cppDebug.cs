using UnityEngine;
using System;
using System.Runtime.InteropServices;

public class cppDebug : MonoBehaviour {

	private delegate void DebugCallback(string message);

	[DllImport("AudioPluginSpatializerTemplate")]
	private static extern void RegisterDebugCallback(DebugCallback callback);

	private void Start()
	{
		RegisterDebugCallback(new DebugCallback(DebugMethod));
	}

	private static void DebugMethod(string message)
	{
		Debug.Log("NAPlugin: " + message);
	}
}