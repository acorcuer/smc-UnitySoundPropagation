#include "AudioPluginUtil.h"
#include "rayTraceUtil.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <future>
extern float hrtfSrcData[];
extern float reverbmixbuffer[];

namespace Spatializer
{
    // Variables that need to be shared between instances of the plugin
    static int numRays = 20;
    static int maxPathLength = 100;
    static int maxNumReflecs = 75;
    static float absCoeff = 0.5f;
    const static float airAbsorbtion[6] = {0.002, 0.005, 0.005, 0.007, 0.012, 0.057};
    const static float C = 343.2;
    const static int impLength = std::ceil(44100 * (maxPathLength/C));
    static GeomeTree* triangleTree;
    static raySphere sourceSphere = raySphere(numRays);
    std::vector<float> rayOutputData;
    static bool newTree = true;
    static bool treeInit = false;
    static bool enableDebug;
    bool textWritten = true;
    
    enum
    {
        P_AUDIOSRCATTN,
        P_FIXEDVOLUME,
        P_CUSTOMFALLOFF,
        P_NUM
    };
    
    struct filterData {
        BiquadFilter Octave1[2];
        BiquadFilter Octave2[2];
        BiquadFilter Octave3[2];
        BiquadFilter Octave4[2];
        BiquadFilter Octave5[2];
        BiquadFilter Octave6[2];
        BiquadFilter Octave7[2];
        BiquadFilter Octave8[2];
    };
    
    struct EffectData
    {
        
        float p[P_NUM];
        Vector3 prevPositons[2];
        std::vector<Ray> sucessfullRays;
        std::vector<float> convOverlapL;
        std::vector<float> convOverlapR;
        union
        {
            filterData filterBank;
            unsigned char pad[(sizeof(filterData) + 15) & ~15]; // This entire structure must be a multiple of 16 bytes (and and instance 16 byte aligned) for PS3 SPU DMA requirements
        };
        
    };
    
    static void SetupFilterCoeffs(filterData* data, float samplerate)
    {
        data->Octave1->SetupLowpass(63, samplerate, 1);
        data->Octave2->SetupBandpass(125, samplerate, 1);
        data->Octave3->SetupBandpass(250, samplerate, 1);
        data->Octave4->SetupBandpass(500, samplerate, 1);
        data->Octave5->SetupBandpass(1000, samplerate, 1);
        data->Octave6->SetupBandpass(2000, samplerate, 1);
        data->Octave7->SetupBandpass(4000, samplerate, 1);
        data->Octave8->SetupHighpass(8000, samplerate, 1);
    }
    
    
    inline bool IsHostCompatible(UnityAudioEffectState* state)
    {
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
        rayOutputData.clear();
        EffectData* effectdata = new EffectData;
        memset(effectdata, 0, sizeof(EffectData));
        state->effectdata = effectdata;
        if (IsHostCompatible(state))
            state->spatializerdata->distanceattenuationcallback = DistanceAttenuationCallback;
        InitParametersFromDefinitions(InternalRegisterEffectDefinition, effectdata->p);
        SetupFilterCoeffs(&effectdata->filterBank, state->samplerate);
        if(enableDebug){
            DebugInUnity(std::string("Spatialiser plugin loaded sucessfully"));
        }
        
        return UNITY_AUDIODSP_OK;
    }
    
    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK ReleaseCallback(UnityAudioEffectState* state)
    {
        if(enableDebug){
            DebugInUnity(std::string("Spatailiser plugin released:"));
        }
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
    
    extern "C" ABA_API void getRayData(long* len, float **data){
        *len = rayOutputData.size();
        auto size = (*len)*sizeof(float);
        *data = static_cast<float*>(malloc(size));
        memcpy(*data, rayOutputData.data(), size);
    }
    
    extern "C" ABA_API void debugToggle(bool state){
        enableDebug = state;
    }
    
    extern "C" ABA_API void setTraceParam(int numberOfRays,int maxLen,int maxReflec, float absorbtionCoeff){
        absCoeff = absorbtionCoeff;
        numRays  = numberOfRays;
        sourceSphere = raySphere(numRays);
        maxPathLength = maxLen;
        maxNumReflecs = maxReflec;
    }
    
    extern "C" ABA_API void marshalGeomeTree(int numNodes,int numTri, int depth,int bbl,float boundingBoxes[],int tl,float triangles[],int lsl,int leafSizes[],int tidl, int triangleIds[],int tml,float triangleMatList[]) {
        rayOutputData.clear();
        std::deque<float> boundingList(boundingBoxes, boundingBoxes + bbl);
        std::deque<float> triangleList(triangles, triangles + tl);
        std::deque<int> leafSizeList(leafSizes, leafSizes + lsl);
        std::deque<int> triangleIdList(triangleIds, triangleIds + tidl);
        std::deque<float> triangleMats(triangleMatList, triangleMatList + tml);
        
        triangleTree = new GeomeTree(depth,&boundingList,&triangleList,&leafSizeList,&triangleIdList,&triangleMats);
        if(enableDebug){
            std::stringstream sstr;
            sstr << "received and constructed tree with ";
            sstr << numNodes;
            sstr << " nodes";
            sendStringStream(&sstr);
        }
        newTree = true;
        treeInit = true;
    }
    
    void addToDebugList(Ray* inRay,float* length) {
        rayOutputData.push_back(inRay->origin.X);
        rayOutputData.push_back(inRay->origin.Y);
        rayOutputData.push_back(inRay->origin.Z);
        rayOutputData.push_back(inRay->direction.X);
        rayOutputData.push_back(inRay->direction.Y);
        rayOutputData.push_back(inRay->direction.Z);
        if(length != NULL) {
            rayOutputData.push_back(fabsf(*length));
        }else{
            rayOutputData.push_back(5.0f);
        }
        
    }
    
    void shootRays(std::vector<Ray> *inputRayList,std::vector<Ray> *outputRayList){
        // For each ray in the raylist, itterate backwards to avoid indexing problems when removing
        // from the list
        for (int i = inputRayList->size()-1; i >= 0; i--) {
            // Get list of triangles the ray might intersect with according to bounding heirarchy
            std::deque<Tri> candidates = triangleTree->getCandidates(&inputRayList->at(i));
            // If candidate list has been populated test the ray against each, otherwise delete!
            if(candidates.size() > 0) {
                float min = INFINITY;
                int minDistIdx = 0;
                for(int j = 0; j < candidates.size();j++){
                    float t;
                    bool intersectResult = inputRayList->at(i).testIntersect(&candidates[j], &t);
                    if(intersectResult){
                        if(t < min && t > kEpsilon) {
                            min = t;
                            minDistIdx = j;
                        }
                    }
                }
                // if minimum distance is not infinity (therefore no valid intersections)
                // update the ray
                if(min != INFINITY){
                    
                    addToDebugList(&inputRayList->at(i), &min);
                    // Update origin
                    inputRayList->at(i).origin = inputRayList->at(i).origin + (inputRayList->at(i).direction * fabsf(min));
                    // Update number of reflections
                    inputRayList->at(i).numReflecs++;
                    // Update path length
                    inputRayList->at(i).pathLength += fabsf(min);
                    //Update angle
                    inputRayList->at(i).updateDirec(candidates[minDistIdx].faceNorm);
                    // update absorbtion
                    inputRayList->at(i).absorbtion *= (1.0f-candidates[minDistIdx].absorbitonCoeff);
                    if(candidates[minDistIdx].objectType != 0){
                        inputRayList->at(i).listenerTag = candidates[minDistIdx].objectType;
                        addToDebugList(&inputRayList->at(i), &min);
                        outputRayList->push_back(inputRayList->at(i));
                        inputRayList->erase(inputRayList->begin()+i);
                    }else{
                        if(inputRayList->at(i).numReflecs >= maxNumReflecs || inputRayList->at(i).pathLength >= maxPathLength || inputRayList->at(i).absorbtion < 0.01){
                            addToDebugList(&inputRayList->at(i), NULL);
                            inputRayList->erase(inputRayList->begin()+i);
                        }
                    }
                }else {
                    addToDebugList(&inputRayList->at(i), NULL);
                    inputRayList->erase(inputRayList->begin()+i);
                }
            }else {
                addToDebugList(&inputRayList->at(i), NULL);
                inputRayList->erase(inputRayList->begin()+i);
            }
        }
        if(inputRayList->size() > 0) {
            shootRays(inputRayList, outputRayList);
        }
    }
    
    void task() {
        std::this_thread::sleep_for(std::chrono::seconds(10));
        std::stringstream sstr;
        sstr << "after process";
        sendStringStream(&sstr);
    }
    
    void calcImpResponse(float* listenerMatrix,float* sourceMatrix, float octavePower[],EffectData* data,std::vector<float>* left,std::vector<float>* right) {
        
        Vector3 sourcePos = Vector3(sourceMatrix[12], sourceMatrix[13], sourceMatrix[14]);
        Vector3 listenerPos = Vector3(listenerMatrix[12], listenerMatrix[13], sourceMatrix[14]);
        float maxL = -INFINITY;
        float maxR = -INFINITY;
        std::vector<int> idxL;
        std::vector<int> idxR;
        bool reShoot = false;
        
        if(sourcePos != data->prevPositons[0]){
            reShoot = true;
            data->prevPositons[0] = sourcePos;
        }
        if(listenerPos != data->prevPositons[1]){
            reShoot = true;
            data->prevPositons[1] = listenerPos;
        }
        
        if(reShoot || newTree ) {
            data->sucessfullRays.clear();
            std::vector<Ray> sourceRays = sourceSphere.getRayList(sourcePos);
            shootRays(&sourceRays,&data->sucessfullRays);
            newTree = false;
            if(enableDebug){
                std::stringstream sstr;
                sstr << data->sucessfullRays.size();
                sstr << " Sucessfull Rays";
                sendStringStream(&sstr);
                
            }
        }
        
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < data->sucessfullRays.size(); j++) {
                int sampIdx = (int)std::round((data->sucessfullRays[j].pathLength/C)*44100);
                if(data->sucessfullRays[j].listenerTag == 1) {
                    if(i == 0) {
                        idxR.push_back(sampIdx);
                    }
                    if(sampIdx > right->size()) {
                        right->resize(sampIdx+1);
                    }
                    right->at(sampIdx) += octavePower[i]*expf(-airAbsorbtion[i]*data->sucessfullRays[j].pathLength)*data->sucessfullRays[j].absorbtion;
                    if(right->at(sampIdx) > maxR) {
                        maxR = right->at(sampIdx);
                    }
                }else {
                    if(i == 0) {
                        idxL.push_back(sampIdx);
                    }
                    if(sampIdx > left->size()){
                        left->resize(sampIdx+1);
                    }
                    left->at(sampIdx) += octavePower[i]*expf(-airAbsorbtion[i]*data->sucessfullRays[j].pathLength)*data->sucessfullRays[j].absorbtion;
                    
                    if(left->at(sampIdx) > maxL) {
                        maxL = left->at(sampIdx);
                    }
                    
                }
            }
        }
        
        int minL = left->size();
        int minR = right->size();
        // Normalise IR and remove front zeros
        for(int i = 0; i < idxL.size(); i++) {
            int idx =idxL[i];
            left->at(idx) = left->at(idx)/maxL ;
            if(idx < minL) {
                minL = idx;
            }
        }
        for(int i = 0; i < idxR.size(); i++) {
            int idx =idxR[i];
            right->at(idx) = right->at(idx)/maxR;
            if(idx < minR) {
                minR = idx;
            }
        }
        int min = std::min(minL, minR);
        left->erase(left->begin(), left->begin()+min);
        right->erase(right->begin(), right->begin()+min);
        
        if(textWritten){
            std::ofstream arrayDataL("/Users/Alex/IRL.txt");
            std::ofstream arrayDataR("/Users/Alex/IRR.txt");
            
            for(int i = 0;i < left->size();i++){
                arrayDataL << left->at(i);
                arrayDataL<< ",";
            }
            for(int i = 0;i < right->size();i++){
                arrayDataR << right->at(i);
                arrayDataR<< ",";
            }
            textWritten = false;
            if(enableDebug) {
                std::stringstream sstr;
                sstr << " saved IR";
                sendStringStream(&sstr);
            }
        }
    };
    
    UNITY_AUDIODSP_RESULT UNITY_AUDIODSP_CALLBACK ProcessCallback(UnityAudioEffectState* state, float* inbuffer, float* outbuffer, unsigned int length, int inchannels, int outchannels)
    {
        
        //        clock_t start;
        //        double diff;
        //        start = clock();
        
        
        // Check that I/O formats are right and that the host API supports this feature
        if (inchannels != 2 || outchannels != 2 ||
            !IsHostCompatible(state) || state->spatializerdata == NULL)
        {
            memcpy(outbuffer, inbuffer, length * outchannels * sizeof(float));
            return UNITY_AUDIODSP_OK;
        }
        
        EffectData* data = state->GetEffectData<EffectData>();
        float octavePower[6] = {0.0f};
        float mono[length];
        UnityComplexNumber compMono[length/2];
        for (unsigned int n = 0; n < length; n++)
        {
            mono[n] = (inbuffer[n * 2] + inbuffer[n * 2 + 1]) / 2.0f;
            compMono[n].re = mono[n];
            compMono[n].im = 0.0f;
            octavePower[0] += powf(data->filterBank.Octave1[0].Process(mono[n]),2.0f);
            octavePower[1] += powf(data->filterBank.Octave2[0].Process(mono[n]),2.0);
            octavePower[2] += powf(data->filterBank.Octave3[0].Process(mono[n]),2.0f);
            octavePower[3] += powf(data->filterBank.Octave4[0].Process(mono[n]),2.0f);
            octavePower[4] += powf(data->filterBank.Octave5[0].Process(mono[n]),2.0f);
            octavePower[5] += powf(data->filterBank.Octave6[0].Process(mono[n]),2.0f);
        }
        
        for(int i = 0 ;  i < 6; i++) {
            octavePower[i] = sqrtf(octavePower[i]/length);
        }
        
        std::vector<float> impLeft;
        std::vector<float> impRight;
        float* l = state->spatializerdata->listenermatrix;
        float* s = state->spatializerdata->sourcematrix;
        bool impCalc = false;
        if(treeInit){
            calcImpResponse(l,s,octavePower,data,&impLeft,&impRight);
            impCalc = true;
        }
                if(impCalc == true) {
                    
                    
                }
        //        std::stringstream sstr2;
        //        sstr2 << diff;
        //        sendStringStream(&sstr2);
        return UNITY_AUDIODSP_OK;
    }
}
