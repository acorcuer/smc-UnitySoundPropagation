// Please note that this will only work on Unity 5.2 or higher.

#include "AudioPluginUtil.h"
#include <string>
#include <vector>
#include <sstream>



extern float hrtfSrcData[];
extern float reverbmixbuffer[];

namespace Spatializer
{
    
    struct meshData{
        int numVerts;
        std::vector<float> verts;
        std::vector<float> norms;
        std::vector<float> size;
    };
    
    meshData geometry = *new meshData;

    //{ Start Interop Functions
    
   // #if UNITY_WIN
    #define ABA_API __declspec(dllexport)
 //   #else
  //  #define ABA_API
  //  #endif
    
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
    
    extern "C" ABA_API void setStruct (int numVertsIn,float vertArrayIn[],float normArrayIn[],float sizeArrayIn[])
    {
        geometry.numVerts = numVertsIn;
        geometry.verts.assign(vertArrayIn, vertArrayIn+numVertsIn*3);
        geometry.norms.assign(normArrayIn, normArrayIn+numVertsIn*3);
        geometry.size.assign(sizeArrayIn, sizeArrayIn+numVertsIn*3);

		std::stringstream sstr;
		sstr << "Recveived mesh with ";
		sstr << geometry.numVerts;
		sstr << " verts";
		std::string s1 = sstr.str();
		DebugInUnity(std::string(s1));
    }
    
    
    // End Interop Functions}
   
    enum
    {
        P_AUDIOSRCATTN,
        P_FIXEDVOLUME,
        P_CUSTOMFALLOFF,
        P_NUM
    };

    const int HRTFLEN = 512;

    const float GAINCORRECTION = 2.0f;

    class HRTFData
    {
        struct CircleCoeffs
        {
            int numangles;
            float* hrtf;
            float* angles;

            void GetHRTF(UnityComplexNumber* h, float angle, float mix)
            {
                int index1 = 0;
                while (index1 < numangles && angles[index1] < angle)
                    index1++;
                if (index1 > 0)
                    index1--;
                int index2 = (index1 + 1) % numangles;
                float* hrtf1 = hrtf + HRTFLEN * 4 * index1;
                float* hrtf2 = hrtf + HRTFLEN * 4 * index2;
                float f = (angle - angles[index1]) / (angles[index2] - angles[index1]);
                for (int n = 0; n < HRTFLEN * 2; n++)
                {
                    h[n].re += (hrtf1[0] + (hrtf2[0] - hrtf1[0]) * f - h[n].re) * mix;
                    h[n].im += (hrtf1[1] + (hrtf2[1] - hrtf1[1]) * f - h[n].im) * mix;
                    hrtf1 += 2;
                    hrtf2 += 2;
                }
            }
        };

    public:
        CircleCoeffs hrtfChannel[2][14];

    public:
        HRTFData()
        {
            float* p = hrtfSrcData;
            for (int c = 0; c < 2; c++)
            {
                for (int e = 0; e < 14; e++)
                {
                    CircleCoeffs& coeffs = hrtfChannel[c][e];
                    coeffs.numangles = (int)(*p++);
                    coeffs.angles = p;
                    p += coeffs.numangles;
                    coeffs.hrtf = new float[coeffs.numangles * HRTFLEN * 4];
                    float* dst = coeffs.hrtf;
                    UnityComplexNumber h[HRTFLEN * 2];
                    for (int a = 0; a < coeffs.numangles; a++)
                    {
                        memset(h, 0, sizeof(h));
                        for (int n = 0; n < HRTFLEN; n++)
                            h[n + HRTFLEN].re = p[n];
                        p += HRTFLEN;
                        FFT::Forward(h, HRTFLEN * 2, false);
                        for (int n = 0; n < HRTFLEN * 2; n++)
                        {
                            *dst++ = h[n].re;
                            *dst++ = h[n].im;
                        }
                    }
                }
            }
        }
    };

    static HRTFData sharedData;

    struct InstanceChannel
    {
        UnityComplexNumber h[HRTFLEN * 2];
        UnityComplexNumber x[HRTFLEN * 2];
        UnityComplexNumber y[HRTFLEN * 2];
        float buffer[HRTFLEN * 2];
    };

	struct Filterbank
	{
		BiquadFilter Octave1[2];
		BiquadFilter Octave2[2];
		BiquadFilter Octave3[2];
		BiquadFilter Octave4[2];
		BiquadFilter Octave5[2];
		BiquadFilter Octave6[2];
		BiquadFilter Octave7[2];
		BiquadFilter Octave8[2];
		BiquadFilter DisplayFilterCoeffs[3];
		float sr;
		Random random;
	};


    struct EffectData
    {
        float p[P_NUM];
        InstanceChannel ch[2];
    };

    inline bool IsHostCompatible(UnityAudioEffectState* state)
    {
        // Somewhat convoluted error checking here because hostapiversion is only supported from SDK version 1.03 (i.e. Unity 5.2) and onwards.
        return
            state->structsize >= sizeof(UnityAudioEffectState) &&
            state->hostapiversion >= UNITY_AUDIO_PLUGIN_API_VERSION;
    }

    int InternalRegisterEffectDefinition(UnityAudioEffectDefinition& definition)
    {

        int numparams = P_NUM;
        definition.paramdefs = new UnityAudioParameterDefinition[numparams];
        RegisterParameter(definition, "AudioSrc Attn", "", 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, P_AUDIOSRCATTN, "AudioSource distance attenuation");
        RegisterParameter(definition, "Fixed Volume", "", 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, P_FIXEDVOLUME, "Fixed volume amount");
        RegisterParameter(definition, "Custom Falloff", "", 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, P_CUSTOMFALLOFF, "Custom volume falloff amount (logarithmic)");
        definition.flags |= UnityAudioEffectDefinitionFlags_IsSpatializer;
        return numparams;
    }

    static UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK DistanceAttenuationCallback(UnityAudioEffectState* state, float distanceIn, float attenuationIn, float* attenuationOut)
    {
        EffectData* data = state->GetEffectData<EffectData>();
        *attenuationOut =
            data->p[P_AUDIOSRCATTN] * attenuationIn +
            data->p[P_FIXEDVOLUME] +
            data->p[P_CUSTOMFALLOFF] * (1.0f / FastMax(1.0f, distanceIn));
        return UNITY_AUDIODSP_OK;
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK CreateCallback(UnityAudioEffectState* state)
    {
        EffectData* effectdata = new EffectData;
        memset(effectdata, 0, sizeof(EffectData));
        state->effectdata = effectdata;
        if (IsHostCompatible(state))
            state->spatializerdata->distanceattenuationcallback = DistanceAttenuationCallback;
        InitParametersFromDefinitions(InternalRegisterEffectDefinition, effectdata->p);
                DebugInUnity(std::string("Spatializer Plugin Initialized"));
        return UNITY_AUDIODSP_OK;
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK ReleaseCallback(UnityAudioEffectState* state)
    {
        EffectData* data = state->GetEffectData<EffectData>();
        delete data;
        return UNITY_AUDIODSP_OK;
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK SetFloatParameterCallback(UnityAudioEffectState* state, int index, float value)
    {
        EffectData* data = state->GetEffectData<EffectData>();
        if (index >= P_NUM)
            return UNITY_AUDIODSP_ERR_UNSUPPORTED;
        data->p[index] = value;
        return UNITY_AUDIODSP_OK;
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK GetFloatParameterCallback(UnityAudioEffectState* state, int index, float* value, char *valuestr)
    {
        EffectData* data = state->GetEffectData<EffectData>();
        if (index >= P_NUM)
            return UNITY_AUDIODSP_ERR_UNSUPPORTED;
        if (value != NULL)
            *value = data->p[index];
        if (valuestr != NULL)
            valuestr[0] = 0;
        return UNITY_AUDIODSP_OK;
    }

    int UNITY_AUDIODSP_CALLBACK GetFloatBufferCallback(UnityAudioEffectState* state, const char* name, float* buffer, int numsamples)
    {
        return UNITY_AUDIODSP_OK;
    }

    static void GetHRTF(int channel, UnityComplexNumber* h, float azimuth, float elevation)
    {
        float e = FastClip(elevation * 0.1f + 4, 0, 12);
        float f = floorf(e);
        int index1 = (int)f;
        if (index1 < 0)
            index1 = 0;
        else if (index1 > 12)
            index1 = 12;
        int index2 = index1 + 1;
        if (index2 > 12)
            index2 = 12;
        sharedData.hrtfChannel[channel][index1].GetHRTF(h, azimuth, 1.0f);
        sharedData.hrtfChannel[channel][index2].GetHRTF(h, azimuth, e - f);
    }

    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK ProcessCallback(UnityAudioEffectState* state, float* inbuffer, float* outbuffer, unsigned int length, int inchannels, int outchannels)
    {
        // Check that I/O formats are right and that the host API supports this feature
        if (inchannels != 2 || outchannels != 2 ||
            !IsHostCompatible(state) || state->spatializerdata == NULL)
        {
            memcpy(outbuffer, inbuffer, length * outchannels * sizeof(float));
            return UNITY_AUDIODSP_OK;
        }

        EffectData* data = state->GetEffectData<EffectData>();

        static const float kRad2Deg = 180.0f / kPI;

        float* m = state->spatializerdata->listenermatrix;
        float* s = state->spatializerdata->sourcematrix;

        // Currently we ignore source orientation and only use the position
        float px = s[12];
        float py = s[13];
        float pz = s[14];

        float dir_x = m[0] * px + m[4] * py + m[8] * pz + m[12];
        float dir_y = m[1] * px + m[5] * py + m[9] * pz + m[13];
        float dir_z = m[2] * px + m[6] * py + m[10] * pz + m[14];

        float azimuth = (fabsf(dir_z) < 0.001f) ? 0.0f : atan2f(dir_x, dir_z);
        if (azimuth < 0.0f)
            azimuth += 2.0f * kPI;
        azimuth = FastClip(azimuth * kRad2Deg, 0.0f, 360.0f);

		/* We compute the distance between the source and the sound
		and the distance decay coefficient*/
		float distance = sqrt((dir_x*dir_x) + (dir_y*dir_y));
		float dis_decay = 1.0f / FastMax(1.0f, distance);

        float elevation = atan2f(dir_y, sqrtf(dir_x * dir_x + dir_z * dir_z) + 0.001f) * kRad2Deg;
        float spatialblend = state->spatializerdata->spatialblend;
        float reverbmix = state->spatializerdata->reverbzonemix;
	
        GetHRTF(0, data->ch[0].h, azimuth, elevation);
        GetHRTF(1, data->ch[1].h, azimuth, elevation);

        // From the FMOD documentation:
        //   A spread angle of 0 makes the stereo sound mono at the point of the 3D emitter.
        //   A spread angle of 90 makes the left part of the stereo sound place itself at 45 degrees to the left and the right part 45 degrees to the right.
        //   A spread angle of 180 makes the left part of the stero sound place itself at 90 degrees to the left and the right part 90 degrees to the right.
        //   A spread angle of 360 makes the stereo sound mono at the opposite speaker location to where the 3D emitter should be located (by moving the left part 180 degrees left and the right part 180 degrees right). So in this case, behind you when the sound should be in front of you!
        // Note that FMOD performs the spreading and panning in one go. We can't do this here due to the way that impulse-based spatialization works, so we perform the spread calculations on the left/right source signals before they enter the convolution processing.
        // That way we can still use it to control how the source signal downmixing takes place.
        float spread = cosf(state->spatializerdata->spread * kPI / 360.0f);
        float spreadmatrix[2] = { 2.0f - spread, spread };

        float* reverb = reverbmixbuffer;
        for (int sampleOffset = 0; sampleOffset < length; sampleOffset += HRTFLEN)
        {
            for (int c = 0; c < 2; c++)
            {
                // stereopan is in the [-1; 1] range, this acts the way fmod does it for stereo
                float stereopan = 1.0f - ((c == 0) ? FastMax(0.0f, state->spatializerdata->stereopan) : FastMax(0.0f, -state->spatializerdata->stereopan));

                InstanceChannel& ch = data->ch[c];

                for (int n = 0; n < HRTFLEN; n++)
                {
                    float left  = inbuffer[n * 2];
                    float right = inbuffer[n * 2 + 1];
                    ch.buffer[n] = ch.buffer[n + HRTFLEN];
                    ch.buffer[n + HRTFLEN] = left * spreadmatrix[c] + right * spreadmatrix[1 - c];
                }

                for (int n = 0; n < HRTFLEN * 2; n++)
                {
                    ch.x[n].re = ch.buffer[n];
                    ch.x[n].im = 0.0f;
                }

                FFT::Forward(ch.x, HRTFLEN * 2, false);

                for (int n = 0; n < HRTFLEN * 2; n++)
                    UnityComplexNumber::Mul<float, float, float>(ch.x[n], ch.h[n], ch.y[n]);

                FFT::Backward(ch.y, HRTFLEN * 2, false);

                for (int n = 0; n < HRTFLEN; n++)
                {
                    float s = inbuffer[n * 2 + c] * stereopan;
					float y = s + (ch.y[n].re * GAINCORRECTION - s) * spatialblend; //* dis_decay;
                    outbuffer[n * 2 + c] = y;
			//		outbuffer[n * 2 + c] = 0;
                    reverb[n * 2 + c] += y * reverbmix;
                }
            }

            inbuffer += HRTFLEN * 2;
            outbuffer += HRTFLEN * 2;
			outbuffer += HRTFLEN * 2;
            reverb += HRTFLEN * 2;
        }

        return UNITY_AUDIODSP_OK;
    }
}
